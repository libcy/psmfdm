#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cufft.h>
#include "psv.h"
#include "const.h"
#include "cuda.h"
#include "plot.h"
#include "malloc.h"

int main(int argc , char *argv[]) {
	//output file name
	char *oname="opsv";
	char *wname="wpsv";
	int i,j;
	int calculate,out,wav;

	if(argc==2){
		calculate=0;
		out=0;
		wav=0;
		for(i=0;i<argv[1][i]!='\0';i++){
			if(argv[1][i]=='o') out=1;
			else if(argv[1][i]=='w') wav=1;
			else if(argv[1][i]=='c') calculate=1;
		}
	}
	else{
		calculate=1;
		out=0;
		wav=1;
	}

	//calculate
	if(calculate){
		FILE *wfile=fopen(wname,"w");
		FILE *ofile=fopen(oname,"w");

		// //dimension

		float **sxx=cudaMat(nx,ny),**sxy=cudaMat(nx,ny),**syy=cudaMat(nx,ny);
		float **den=cudaMat(nx2,ny2),**rig=cudaMat(nx2,ny2),**lam=cudaMat(nx2,ny2);

		float **ux=cudaMat(nx,ny),**uy=cudaMat(nx,ny);
		float **vx=cudaMat(nx,ny),**vy=cudaMat(nx,ny);
		float **up=floatMat(nx,ny),**us=floatMat(nx,ny);

		float **dxux=cudaMat(nx,ny),**dxuy=cudaMat(nx,ny);
		float **dyux=cudaMat(nx,ny),**dyuy=cudaMat(nx,ny);
		float **dxvx=cudaMat(nx,ny),**dxvy=cudaMat(nx,ny);
		float **dyvx=cudaMat(nx,ny),**dyvy=cudaMat(nx,ny);
		float **dxsxx=cudaMat(nx,ny),**dxsxy=cudaMat(nx,ny);
		float **dysxy=cudaMat(nx,ny),**dysyy=cudaMat(nx,ny);

		float **ggg=cudaMat(nx,ny);
		float **dvp=floatMat(2,nx),**dvs=floatMat(2,nx),**dden=floatMat(2,nx);
		float **cxwork=cudaMat(nx,ny),**cywork=cudaMat(nx,ny);
		float **uxall=floatMat(nst,ntskp),**uyall=floatMat(nst,ntskp);
		float *gx=floatVec(nx),*gy=floatVec(ny);

		int *istx=intVec(nst),*isty=intVec(nst),**imap=intMat(nx,ny);

		for(i=0;i<nst;i++){
			istx[i]=i*4+1;
			isty[i]=na+1;
		}

		// velocity structure

		float VPB = 6.9;
		float VSB = 3.9;
		float ROB = 3.1;

		float RRIGB = ROB * VSB*VSB;
		float RLANB = ROB * VPB*VPB - 2.0 * RRIGB;
		float RDENB = ROB;

		for(i=0;i<nx2;i++){
			for(j=0;j<ny2;j++){
				rig[i][j]=RRIGB;
				den[i][j]=RDENB;
				lam[i][j]=RLANB;
			}
		}
		for(i=0;i<nx;i++){
			for(j=0;j<ny;j++){
				imap[i][j]=0;
			}
		}

		for(i=0;i<nst;i++){
			imap[istx[i]][isty[i]]=7;
		}

		// initialize
		int kx=nbegi2(nx);
		int ky=nbegi2(ny);

		for(i=0;i<nx;i++){
			gx[i]=dx*(i+1);
		}
		for(i=0;i<ny;i++){
			gy[i]=dy*(i+1);
		}

		float ftmax=t0+at*2;

		clear(vx,nx,ny,0.0);
		clear(vy,nx,ny,0.0);
		clear(ux,nx,ny,0.0);
		clear(uy,nx,ny,0.0);
		clear(sxx,nx,ny,0.0);
		clear(sxy,nx,ny,0.0);
		clear(syy,nx,ny,0.0);

		// absorbing boundary confition
		float apara=0.015;
		float gg;
		for(i=0;i<nx;i++){
			for(j=0;j<ny;j++){
				if(i+1<nxa){
					gg=exp(-pow(apara*(nxa-i-1),2));
				}
				else if(i+1>(nx-nxa+1)){
					gg=exp(-pow(apara*(i-nx+nxa),2));
				}
				else if(j+1>(ny-nya+1)){
					gg=exp(-pow(apara*(j-ny+nya),2));
				}
				else{
					gg=1.0;
				}
				ggg[i][j]=gg;
			}
		}

		cufftHandle plan;
		cufftComplex *data;
		int dimension[1]={nx};
		cufftPlanMany(&plan,1,dimension,NULL,1,1,NULL,1,1,CUFFT_C2C,ny*2);
		cudaMallocManaged((void**)&data, sizeof(cufftComplex)*nx*ny*2);


		//time step start
		int ntw=0;
		int ntt=0;
		float t;
		clock_t start0;
		float c0=9.0/8.0;
		float c1=1.0/24.0;

		start0=clock();
		for(int it=0;it<ntmax;it++){
			if(it%((int)ntmax/10)==0) printf("%d%%\n",10*it/((int)ntmax/10));
			ntt++;
			t=dt*it;
			ntw++;

			diffxspm(vx,dxvx,vy,dxvy,plan,data,nx,ny,dx);
			cudaFinidyyx<<<2*nx*nbt,ny/nbt>>>(vy,dyvy,vx,dyvx,nx,ny,dx,dy,dt,c0,c1);
			cudaDeviceSynchronize();
			cudaPrep<<<nx*nbt,ny/nbt>>>(sxx,syy,sxy,lam,rig,ggg,dxvx,dxvy,dyvx,dyvy);
			cudaDeviceSynchronize();
			diffxspm(sxy,dxsxy,sxx,dxsxx,plan,data,nx,ny,dx);
			cudaFinidyyx<<<2*nx*nbt,ny/nbt>>>(sxy,dysxy,syy,dysyy,nx,ny,dx,dy,dt,c0,c1);
			cudaDeviceSynchronize();
			cudaCalc<<<nx*nbt,ny/nbt>>>(vx,vy,ux,uy,dxsxx,dxsxy,dysxy,dysyy,ggg,den,t,
				ftmax,rmxx,rmxy,rmyx,rmyy,fxx,fzz,dpxx,dpzz);
			cudaDeviceSynchronize();
			
			if(ntt==nskip){
				int isx,isy,it1;
				for(int ns=0;ns<nst;ns++){
					ntt=0;
					isx=istx[ns]-1;
					isy=isty[ns]-1;
					it1=(it+1)/nskip;

					uxall[ns][it1]=ux[isx][isy];
					uyall[ns][it1]=uy[isx][isy];
				}
			}

			if(ntw==nwrite){
				ntw=0;

				diffxspm(ux,dxux,uy,dxuy,plan,data,nx,ny,dx);
				cudaFinidyyx<<<2*nx*nbt,ny/nbt>>>(uy,dyuy,ux,dyux,nx,ny,dx,dy,dt,c0,c1);
				cudaDeviceSynchronize();

				for(i=0;i<nx;i++){
					for(j=0;j<ny;j++){
						up[i][j]=dxux[i][j]+dyuy[i][j];
						us[i][j]=dxuy[i][j]-dyux[i][j];
					}
				}

				fprintMat(ofile,up,nx,ny);
				fprintMat(ofile,us,nx,ny);
			}
		}

		fprintMat(wfile,uxall,nst,ntskp);
		fprintMat(wfile,uyall,nst,ntskp);

		printf("100%%\n%.2f\n",(double)(clock()-start0)/CLOCKS_PER_SEC);
	}
		
	if(out){
		snapPSV(oname);
	}

	if(wav){
		wavePSV(wname);
	}

	return 0;
}

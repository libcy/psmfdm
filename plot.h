/*
 * psv.h
 *
 *  Created on: Mar 7, 2014
 *      Author: widget
 */

#ifndef PLOT_H_
#define PLOT_H_

#include "const.h"
#include "malloc.h"

void snapPSV(char *filename){
	FILE *file=fopen(filename,"r");

	float **up=floatMat(nx,ny),**us=floatMat(nx,ny);
	float **u=floatMat(nx,ny),**p=floatMat(nx,ny),**s=floatMat(nx,ny);
	int n[5]={0,1,2,3,4};

	FILE **snapshot=(FILE **)malloc(5*sizeof(FILE *));
	*snapshot=fopen("snap1","w");
	*(snapshot+1)=fopen("snap2","w");
	*(snapshot+2)=fopen("snap3","w");
	*(snapshot+3)=fopen("snap4","w");
	*(snapshot+4)=fopen("snap5","w");

	int i,j;
	float pmax,smax,cp,lp,cs,ls,x,y;

	for(int isnap=0;isnap<nsnap;isnap++){
		for(i=0;i<nx;i++){
			for(j=0;j<ny;j++){
				u[i][j]=0;
			}
		}
		fscanfMat(file,up,nx,ny);
		fscanfMat(file,us,nx,ny);

		pmax=0.0;
		smax=0.0;

		for(i=0;i<nx;i++){
			for(j=0;j<ny;j++){
				if(pmax<abs(up[i][j])){
					pmax=abs(up[i][j]);
				}
				if(smax<abs(us[i][j])){
					smax=abs(us[i][j]);
				}
			}
		}
		// printf("Pmax=%f Smax=%f\n",pmax,smax);

		for(i=0;i<nx;i++){
			for(j=0;j<ny;j++){
				cp=pamp;
				lp=0.1*pmax;
				if(abs(up[i][j])>cp&&up[i][j]<0.0){
					up[i][j]=-cp;
				}
				else if(abs(up[i][j])>cp&&up[i][j]>0.0){
					up[i][j]=cp;
				}
				if(abs(us[i][j])<lp){
					up[i][j]=0.0;
				}
			}
		}

		for(i=0;i<nx;i++){
			for(j=0;j<ny;j++){
				cs=samp;
				ls=0.1*smax;
				if(abs(us[i][j])>cs&&us[i][j]<0.0){
					us[i][j]=-cs;
				}
				else if(abs(us[i][j])>cs&&us[i][j]>0.0){
					us[i][j]=cs;
				}
				if(abs(us[i][j])<ls){
					us[i][j]=0.0;
				}
			}
		}

		if(isnap==n[0]||isnap==n[1]||isnap==n[2]||isnap==n[3]||isnap==n[4]){
			for(j=0;j<ny;j++){
				for(i=0;i<nx;i++){
					x=i*dx;
					y=j*dy;
					p[i][j]=up[i][j]/pampall;
					s[i][j]=us[i][j]/sampall;
					// if(up[i][j]>1e-5||us[i][j]>1e-5){
						// printf("%f %f\n", up[i][j],us[i][j]);
					// }
				}
			}
			for(j=0;j<ny;j++){
				for(i=0;i<nx;i++){
					x=i*dx;
					y=j*dy;
					if(abs(s[i][j])>abs(p[i][j])){
						u[i][j]=-abs(s[i][j]);
					}
					else if(abs(p[i][j])>abs(s[i][j])){
						u[i][j]=abs(s[i][j]);
					}
					fprintf(*(snapshot+isnap), "%f %f %f\n", x,y,u[i][j]);
				}
			}
		}
	}

	fclose(file);
	fclose(*(snapshot));
	fclose(*(snapshot+1));
	fclose(*(snapshot+2));
	fclose(*(snapshot+3));
	fclose(*(snapshot+4));
}

void wavePSV(char *filename){
	int ndskip=1;
	float dt2=dt*10,dx2=dx*4;
	float **ux=floatMat(nst,ntskp),**uz=floatMat(nst,ntskp);

	FILE *file=fopen(filename,"r");

	FILE *filex=fopen("ux","w");
	FILE *filez=fopen("uz","w");

	fscanfMat(file,ux,nst,ntskp);
	fscanfMat(file,uz,nst,ntskp);

	int i,j;
	float tm,shift;
	for(i=0;i<nst;i+=nsskip){
		fprintf(filex, ">\n");
		fprintf(filez, ">\n");
		for(j=0;j<ntskp;j+=ndskip){
			tm=j*dt2;
			shift=i*dx2;
			fprintf(filex, "%f %f\n", tm,ux[i][j]*15.0+shift);
			fprintf(filez, "%f %f\n", tm,uz[i][j]*15.0+shift);
		}
	}

	fclose(file);
}

#endif /* PSV_H_ */

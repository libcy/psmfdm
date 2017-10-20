#ifndef CUDA_H_
#define CUDA_H_
__device__ int cudaIntPow(int x,int y){
	int result=1;
	for(int i=0;i<y;i++){
		result*=x;
	}
	return result;
}
__device__ int cudaBitrev(float *a,int l){
	//preparation
	if(l<=0||l>=23) return 30000;
	float itest[20],inc[20];

	int nn=0,nr=0;
	int i=l-1;
	int m=cudaIntPow(2,i);
	int k=2;
	int mr,mn;
	float w;

	if(i-2>=0){
		while(i-2>0){
			i--;
			itest[i-2]=m-k;
			k+=k;
			inc[i-2]=k-itest[i-2];
		}
		mr=m+nr;
		if(nr-nn>0){
			w=a[nn];
			a[nn]=a[nr];
			a[nr]=w;

			mn=m+nn;
			w=a[mn+1];
			a[mn+1]=a[mr+1];
			a[mr+1]=w;
		}
		nn+=2;
		w=a[nn-1];
		a[nn-1]=a[mr];
		a[mr]=w;
		nr=k+nr;
	}
	while(1){
		mr=m+nr;
		if(nr-nn>0){
			w=a[nn];
			a[nn]=a[nr];
			a[nr]=w;

			mn=m+nn;
			w=a[mn+1];
			a[mn+1]=a[mr+1];
			a[mr+1]=w;
		}
		nn+=2;
		w=a[nn-1];
		a[nn-1]=a[mr];
		a[mr]=w;

		i=1;
		if(nn-m>=0) return 0;
		while(nr-itest[i-1]>=0) i++;
		nr=inc[i-1]+nr;
		mr=m+nr;
		if(nr-nn>0){
			w=a[nn];
			a[nn]=a[nr];
			a[nr]=w;

			mn=m+nn;
			w=a[mn+1];
			a[mn+1]=a[mr+1];
			a[mr+1]=w;
		}
		nn+=2;
		w=a[nn-1];
		a[nn-1]=a[mr];
		a[mr]=w;
		nr=k+nr;
	}
}
__device__ int cudaFftr(float *a,int m){
	if(m<0) return 30000;
	if(m==0) return 0;
	int n2=1;
	if(m>1){
		int i,j,k,l;
		double ss,c,s,c0,s0,cn,sn;
		cudaBitrev(a,m);
		int n=cudaIntPow(2,m);
		float f=2.0/n;
		float p,q;
		int n1,n3,n4,l0,l1,l2,l3,n2k,n3k;
		for(l=0;l<n;l+=2){
			p=a[l+1];
			a[l+1]=(a[l]-p)*f;
			a[l]=(a[l]+p)*f;
		}
		n1=1;
		n2=2;
		n3=4;
		n4=8;
		for(i=0;i<m-2;i++){
			for(l=0;l<n;l+=n3){
				l1=l+n2;
				p=a[l];
				a[l]=p+a[l1];
				a[l1]=p-a[l1];
			}
			c=cos(pi/cudaIntPow(2,i+2));
			s=sin(pi/cudaIntPow(2,i+2));
			ss=s+s;
			c0=1.0;
			s0=0.0;
			for(k=1;k<=n1-1;k++){
				n3k=n3-k;
				n2k=n2-k;
				for(j=1;j<=n;j+=n4){
					l0=j+k;
					l2=l0+n3;
					l1=j+n3k;
					l3=l1+n3;
					l0--;l1--;l2--;l3--;
					p=c*a[l2]-s*a[l3];
					q=s*a[l2]+c*a[l3];
					a[l2]=q-a[l1];
					a[l3]=q+a[l1];
					a[l1]=a[l0]-p;
					a[l0]=a[l0]+p;
					l1=l0+n2;
					l3=l1+n3;
					l0=j+n2k-1;
					l2=l0+n3;
					p=a[l2]*s-a[l3]*c;
					q=a[l2]*c+a[l3]*s;
					a[l2]=q-a[l1];
					a[l3]=q+a[l1];
					a[l1]=a[l0]-p;
					a[l0]=a[l0]+p;
				}
				cn=c0-ss*s;
				sn=ss*c+s0;
				c0=c;
				c=cn;
				s0=s;
				s=sn;
			}
			for(j=0;j<n;j+=n4){
				l0=j+n1;
				l2=l0+n3;
				l1=l2-n2;
				l3=l1+n3;
				p=a[l2]*s-a[l3]*c;
				q=a[l2]*c+a[l3]*s;
				a[l2]=q-a[l1];
				a[l3]=q+a[l1];
				a[l1]=a[l0]-p;
				a[l0]=a[l0]+p;
			}
			n1=n2;
			n2=n3;
			n3=n4;
			n4+=n4;
		}
		for(j=1;j<n1;j++){
			p=a[n2+j];
			a[n2+j]=a[n-j];
			a[n-j]=p;
			
		}
	}
	a[0]=(a[0]+a[n2])*0.5;
	a[n2]=a[0]-a[n2];

	return 0;
}
__device__ int cudaFftri(float *a,int m){
	if(m<0) return 30000;
	if(m==0) return 0;
	int n=cudaIntPow(2,m);
	int n4=n;
	int n3=n/2;
	int n2=n3/2;
	int n1=n2/2;
	int i,j,k,l;

	float p=a[0],q;
	a[0]=p+a[n3];
	a[n3]=p-a[n3];
	if(m>1){
		double c,s,c0,s0,cn,sn,ss;
		int n2k,n3k,l0,l1,l2,l3;
		for(j=1;j<n2;j++){
			p=a[n3+j];
			a[n3+j]=a[n-j];
			a[n-j]=p;
		}
		for(i=m-3;i>=-1;i--){
			for(l=0;l<n;l+=n3){
				l1=l+n2;
				p=a[l];
				a[l]=p+a[l1];
				a[l1]=p-a[l1];
			}
			if(n1==0) break;
			c=cos(pi/cudaIntPow(2,i+2));
			s=sin(pi/cudaIntPow(2,i+2));
			ss=s+s;
			c0=1.0;
			s0=0.0;
			for(k=1;k<=n1-1;k++){
				n3k=n3-k;
				n2k=n2-k;
				for(j=1;j<=n;j+=n4){
					l0=j+k;
					l2=l0+n3;
					l1=j+n3k;
					l3=l1+n3;
					l0--;l1--;l2--;l3--;
					p=a[l0]-a[l1];
					q=a[l2]+a[l3];
					a[l0]=a[l0]+a[l1];
					a[l1]=a[l3]-a[l2];
					a[l2]=c*p+s*q;
					a[l3]=c*q-s*p;

					l1=l0+n2;
					l3=l1+n3;
					l0=j+n2k-1;
					l2=l0+n3;
					p=a[l0]-a[l1];
					q=a[l2]+a[l3];
					a[l0]=a[l0]+a[l1];
					a[l1]=a[l3]-a[l2];
					a[l2]=s*p+c*q;
					a[l3]=s*q-c*p;
				}
				cn=c0-ss*s;
				sn=ss*c+s0;
				c0=c;
				c=cn;
				s0=s;
				s=sn;
			}
			for(j=1;j<=n;j+=n4){
				l0=j+n1-1;
				l2=l0+n3;
				l1=l2-n2;
				l3=l1+n3;
				p=a[l0]-a[l1];
				q=a[l2]+a[l3];
				a[l0]=a[l0]+a[l1];
				a[l1]=a[l3]-a[l2];
				a[l2]=p*s+q*c;
				a[l3]=q*s-p*c;
			}
			n4=n3;n3=n2;
			n2=n1;
			n1/=2;
		}
	}
	return cudaBitrev(a,m);
}
__global__ void cudaDiffxspm(float **pa,float **pb,float **ma,float **mb,
	float **pwork,float **mwork,int nx,int ny,int ns,int kx,float dx){
	int i,j;
	float **a,**b,**work;
	if(blockIdx.x<nbt){
		a=ma;b=mb;work=mwork;j=threadIdx.x+blockIdx.x*ny/nbt;
	}
	else{
		a=pa;b=pb;work=pwork;j=threadIdx.x+(blockIdx.x-nbt)*ny/nbt;
	}
	if(j+1<=ns){
		for(i=0;i<nx;i++){
			a[i][j]=0.0;
		}
	}
	else{
		for(i=0;i<nx;i++){
			work[j][i]=a[i][j];
		}

		cudaFftr(work[j],kx);

		float dkx=6.283185/nx/dx;
		float sft=6.283185/nx/2.0;
		int nxd2=nx/2;

		float ws,wc;

		work[j][0]=0.0;
		work[j][nxd2]=0.0;

		for(i=1;i<nxd2;i++){
			if(blockIdx.x<nbt){
				wc=-dkx*i*(work[j][i]*cos(sft*i)-work[j][nxd2+i]*sin(sft*i));
				ws=dkx*i*(work[j][nxd2+i]*cos(sft*i)+work[j][i]*sin(sft*i));
			}
			else{
				wc=-dkx*i*(work[j][i]*cos(sft*i)+work[j][nxd2+i]*sin(sft*i));
				ws=dkx*i*(work[j][nxd2+i]*cos(sft*i)-work[j][i]*sin(sft*i));
			}
			work[j][i]=ws;
			work[j][nxd2+i]=wc;
		}

		cudaFftri(work[j],kx);

		for(i=0;i<nx;i++){
			b[i][j]=work[j][i];
		}
	}
}
__global__ void cudaFinidyyx(float **ya,float **dyya,float **xa,float **dyxa,int nx,int ny,
	float dx,float dy,float dt,float c0,float c1){
	int i=blockIdx.x%(2*nx),j=threadIdx.x+(blockIdx.x-i)/(2*nx)*ny/nbt;
	if(i<nx){
		if(j==0){
			dyya[i][j]=1.0/dy*(c0*ya[i][0]-c1*ya[i][1]);
		}
		else if(j==1){
			dyya[i][j]=1.0/dy*(c0*(ya[i][1]-ya[i][0])-c1*ya[i][1]);
		}
		else if(j==ny-1){
			dyya[i][j]=1.0/dy*(c0*(ya[i][ny-1]-ya[i][ny-2])+c1*ya[i][ny-3]);
		}
		else{
			dyya[i][j]=1.0/dy*(c0*(ya[i][j]-ya[i][j-1])-c1*(ya[i][j+1]-ya[i][j-2]));
		}
	}
	else{
		i-=nx;
		if(j==0){
			dyxa[i][j]=1.0/dy*(c0*(xa[i][1]-xa[i][0])-c1*xa[i][2]);
		}
		else if(j==ny-2){
			dyxa[i][j]=1.0/dy*(c0*(xa[i][ny-1]-xa[i][ny-2])+c1*xa[i][ny-3]);
		}
		else if(j==ny-1){
			dyxa[i][j]=1.0/dy*(c0*(-xa[i][ny-1])+c1*xa[i][ny-2]);
		}
		else{
			dyxa[i][j]=1.0/dy*(c0*(xa[i][j+1]-xa[i][j])-c1*(xa[i][j+2]-xa[i][j-1]));
		}
	}
}
__device__ float cudaDherrman(float a,float x,float x0){
	float a2=2.0*a;
	float t=x-x0;
	float td=(t+a2)/a;
	if(t<=-a2) return 0.0;
	if(t<=-a) return td/(a2*a);
	if(t<=a) return (-td+2.0)/(a2*a);
	if(t<=a2) return (td-4.0)/(a2*a);
	return 0.0;
}
__device__ float cudaHerrman(float a,float x,float x0){
	float a2=2.0*a;
	float t=x-x0;
	float td=(t+a2)/a;
	if(t<=-a2) return 0.0;
	if(t<=-a) return (0.5*td*td)/a2;
	if(t<=a) return (-0.5*td*td+2.0*td-1.0)/a2;
	if(t<=a2) return (0.5*td*td-4.0*td+8.0)/a2;
	return 0.0;
}
__device__ float cudaFxmxz(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at,float xs,float zs){
	float x0=i0*dx+xs;
	float z0=j0*dz+zs;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=cudaHerrman(ax,x,x0);
	float fhz=-cudaDherrman(az,z,z0);
	float fht=cudaHerrman(at,t,t0);
	return fhx*fhz*fht;
}
__device__ float cudaFzmxz(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at,float xs,float zs){
	float x0=i0*dx+xs;
	float z0=j0*dz+zs;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=-cudaDherrman(ax,x,x0);
	float fhz=cudaHerrman(az,z,z0);
	float fht=cudaHerrman(at,t,t0);
	return fhx*fhz*fht;
}
__device__ float cudaFzmzz(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at,float xs,float zs){
	float x0=i0*dx+xs;
	float z0=j0*dz+zs;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=cudaHerrman(ax,x,x0);
	float fhz=-cudaDherrman(az,z,z0);
	float fht=cudaHerrman(at,t,t0);
	return fhx*fhz*fht;
}
__device__ float cudaFxmxx(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at,float xs,float zs){
	float x0=i0*dx+xs;
	float z0=j0*dz+zs;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=-cudaDherrman(ax,x,x0);
	float fhz=cudaHerrman(az,z,z0);
	float fht=cudaHerrman(at,t,t0);
	return fhx*fhz*fht;
}
__device__ float cudaFx(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at){
	float x0=i0*dx;
	float z0=j0*dz;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=cudaHerrman(ax,x,x0);
	float fhz=cudaHerrman(az,z,z0);
	float fht=cudaHerrman(at,t,t0);
	return fhx*fhz*fht;
}
__device__ float cudaFz(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at){
	float x0=i0*dx;
	float z0=j0*dz;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=cudaHerrman(ax,x,x0);
	float fhz=cudaHerrman(az,z,z0);
	float fht=cudaHerrman(at,t,t0);
	return fhx*fhz*fht;
}
__device__ float cudaExforce(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at){
	float x0=i0*dx;
	float z0=j0*dz;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=-cudaDherrman(ax,x,x0);
	float fhz=cudaHerrman(az,z,z0);
	float fht=cudaHerrman(at,t,t0);
	return fhx*fhz*fht;
}
__device__ float cudaEzforce(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at){
	float x0=i0*dx;
	float z0=j0*dz;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=cudaHerrman(ax,x,x0);
	float fhz=-cudaDherrman(az,z,z0);
	float fht=cudaHerrman(at,t,t0);
	return fhx*fhz*fht;
}
__global__ void cudaCalc(float **vx,float **vy,float **ux,float **uy,float **dxsxx,float **dxsxy,float **dysxy,float **dysyy,
	float **ggg,float **den,float t,float ftmax,float rmxx,float rmxy,float rmyx,float rmyy,float fxx,float fzz,float dpxx,float dpzz){
	int i=blockIdx.x%nx,j=threadIdx.x+(blockIdx.x-i)/nx*ny/nbt;

	float gg=ggg[i][j];
	float denvx=den[i*2][j*2];
	float denvy=den[i*2+1][j*2+1];

	float fx1,fy1;
	if(t<ftmax){
		fx1=rmxx*cudaFxmxx(i,j,ii0,jj0,dx,dy,ax,ay,t,t0,at,0.0,0.0)+
			rmxy*cudaFxmxz(i,j,ii0,jj0,dx,dy,ax,ay,t,t0,at,0.0,0.0)+
			fxx*cudaFx(i,j,ii0,jj0,dx,dy,ax,ay,t,t0,at)+
			dpxx*cudaExforce(i,j,ii0,jj0,dx,dy,ax,ay,t,t0,at);
		fy1=rmyx*cudaFzmxz(i,j,ii0,jj0,dx,dy,ax,ay,t,t0,at,-dx/2,-dy/2)+
			rmyy*cudaFzmzz(i,j,ii0,jj0,dx,dy,ax,ay,t,t0,at,-dx/2,-dy/2)+
			fzz*cudaFz(i,j,ii0,jj0,dx,dy,ax,ay,t,t0,at)+
			dpzz*cudaEzforce(i,j,ii0,jj0,dx,dy,ax,ay,t,t0,at);

	}
	else{
		fx1=0.0;
		fy1=0.0;
	}

	float uxt2ij=(dxsxx[i][j]+dysxy[i][j]+fx1)/denvx;
	float uyt2ij=(dxsxy[i][j]+dysyy[i][j]+fy1)/denvy;

	vx[i][j]=vx[i][j]*gg+dt*uxt2ij;
	vy[i][j]=vy[i][j]*gg+dt*uyt2ij;

	ux[i][j]=ux[i][j]*gg+dt*vx[i][j];
	uy[i][j]=uy[i][j]*gg+dt*vy[i][j];
}
__global__ void cudaPrep(float **sxx,float **syy,float **sxy,float **lam,float **rig,float **ggg,
	float **dxvx,float **dxvy,float **dyvx,float **dyvy){
	int i=blockIdx.x%nx,j=threadIdx.x+(blockIdx.x-i)/nx*ny/nbt;
	float ram1=lam[i*2+1][j*2];
	float rig1=rig[i*2+1][j*2];
	float rig2=rig[i*2][j*2+1];
	float gg=ggg[i][j];

	float sxxt1ij=(ram1+2.0*rig1)*dxvx[i][j]+ram1*dyvy[i][j];
	float syyt1ij=(ram1+2.0*rig1)*dyvy[i][j]+ram1*dxvx[i][j];
	float sxyt1ij=rig2*(dxvy[i][j]+dyvx[i][j]);

	sxx[i][j]=sxx[i][j]*gg+dt*sxxt1ij;
	syy[i][j]=syy[i][j]*gg+dt*syyt1ij;
	sxy[i][j]=sxy[i][j]*gg+dt*sxyt1ij;

	if(j==0) syy[i][na]=0.0;
}
__global__ void diffin(float **vx,float **vy,cufftComplex *data,int nx,int ny){
	int i=blockIdx.x%nx,j=threadIdx.x+(blockIdx.x-i)/nx*ny/nbt;
	data[i*ny+j].x=vx[j][i];
	data[i*ny+j+nx*ny].x=vy[j][i];
	data[i*ny+j].y=0;
	data[i*ny+j+nx*ny].y=0;
}
__global__ void diffout(float **dxvx,float **dxvy,cufftComplex *data,int nx,int ny){
	int i=blockIdx.x%nx,j=threadIdx.x+(blockIdx.x-i)/nx*ny/nbt;
	dxvx[j][i]=data[i*ny+j].x/nx;
	dxvy[j][i]=data[i*ny+j+nx*ny].x/nx;
}
__global__ void diffprocess(cufftComplex *data,int nx,int ny,float dkx,float sft){
	float ws,wc;
	int i=blockIdx.x%nx,j=threadIdx.x+(blockIdx.x-i)/nx*ny/nbt;

	if(j<nx/2){
		if(j==0){
			data[i*ny].x=0;
			data[i*ny].y=0;
			data[i*ny+nx/2].x=0;
			data[i*ny+nx/2].y=0;
		}
		else{
			wc=-dkx*j*(data[i*ny+j].x*cos(sft*j)-data[i*ny+j].y*sin(sft*j));
			ws=dkx*j*(-data[i*ny+j].y*cos(sft*j)-data[i*ny+j].x*sin(sft*j));
			data[i*ny+j].x=ws;
			data[i*ny+j].y=-wc;
			data[i*ny+ny-j].x=ws;
			data[i*ny+ny-j].y=wc;
		}
	}
	else{
		j-=nx/2;
		if(j==0){
			data[i*ny+nx*ny].x=0;
			data[i*ny+nx*ny].y=0;
			data[i*ny+nx/2+nx*ny].x=0;
			data[i*ny+nx/2+nx*ny].y=0;
		}
		else{
			wc=-dkx*j*(data[nx*ny+i*ny+j].x*cos(sft*j)+data[nx*ny+i*ny+j].y*sin(sft*j));
			ws=dkx*j*(-data[nx*ny+i*ny+j].y*cos(sft*j)+data[nx*ny+i*ny+j].x*sin(sft*j));
			data[i*ny+j+nx*ny].x=ws;
			data[i*ny+j+nx*ny].y=-wc;
			data[i*ny+ny-j+nx*ny].x=ws;
			data[i*ny+ny-j+nx*ny].y=wc;
		}
	}
}

void diffxspm(float **vx,float **dxvx,float **vy,float **dxvy,cufftHandle plan,
	cufftComplex *data,int nx,int ny,float dx){
	diffin<<<ny*nbt,nx/nbt>>>(vx,vy,data,nx,ny);
	cufftExecC2C(plan, data, data, CUFFT_FORWARD);
	cudaDeviceSynchronize();
	diffprocess<<<ny*nbt,nx/nbt>>>(data,nx,ny,6.283185/nx/dx,6.283185/nx/2.0);
	cufftExecC2C(plan, data, data, CUFFT_INVERSE);
	cudaDeviceSynchronize();
	diffout<<<ny*nbt,nx/nbt>>>(dxvx,dxvy,data,nx,ny);
}
float *cudaVec(int m){
	float *vec;
	cudaMallocManaged((void **)&vec,m*sizeof(float));
	return vec;
}
float **cudaMat(int m,int n){
	float **mat;
	float *data;
	cudaMallocManaged((void **)&data,m*n*sizeof(float));
	cudaMallocManaged((void **)&mat,m*sizeof(float *));
	for(int i=0;i<m;i++){
		mat[i]=data+n*i;
	}
	return mat;
}



#endif /* CUDA_H_ */

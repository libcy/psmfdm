/*
 * psv.h
 *
 *  Created on: Mar 7, 2014
 *      Author: widget
 */

#ifndef PSV_H_
#define PSV_H_
#include "malloc.h"
#include "const.h"

int intPow(int x,int y){
	int result=1;
	for(int i=0;i<y;i++){
		result*=x;
	}
	return result;
}

int nbegi2(int n){
	int nxx=n;
	int kx=0;
	for(int i=0;i<n;i++){
		kx++;
		nxx/=2;
		if(nxx==1) break;
	}
	return kx;
}

void clear(float **a,int nx,int ny,float value){
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			a[i][j]=value;
		}
	}
}

int bitrev(float *a,int l){
	//preparation
	if(l<=0||l>=23) return 30000;
	float itest[20],inc[20];

	int nn=0,nr=0;
	int i=l-1;
	int m=intPow(2,i);
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

int fftr(float *a,int m){
	if(m<0) return 30000;
	if(m==0) return 0;
	int n2=1;
	if(m>1){
		int i,j,k,l;
		double dc[32],ds[32],ss,c,s,c0,s0,cn,sn;
		for(i=0;i<32;i++){
			dc[i]=cos(pi/intPow(2,i+2));
			ds[i]=sin(pi/intPow(2,i+2));
		}
		bitrev(a,m);
		int n=intPow(2,m);
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
			c=dc[i];
			s=ds[i];
			ss=ds[i]+ds[i];
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

int fftri(float *a,int m){
	if(m<0) return 30000;
	if(m==0) return 0;
	int n=intPow(2,m);
	int n4=n;
	int n3=n/2;
	int n2=n3/2;
	int n1=n2/2;
	int i,j,k,l;

	float p=a[0],q;
	a[0]=p+a[n3];
	a[n3]=p-a[n3];
	if(m>1){
		double dc[32],ds[32];
		double c,s,c0,s0,cn,sn,ss;
		int n2k,n3k,l0,l1,l2,l3;
		for(i=0;i<32;i++){
			dc[i]=cos(pi/intPow(2,i+2));
			ds[i]=sin(pi/intPow(2,i+2));
		}
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
			c=dc[i];
			s=ds[i];
			ss=ds[i]+ds[i];
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
	return bitrev(a,m);
}

void diffxsp(float **a,float **b,float *work,int nx,int ny,int ns,int kx,float dx){
	int i,j;
	
	for(j=0;j<ny;j++){
		if(j+1<=ns){
			for(i=0;i<nx;i++){
				a[i][j]=0.0;
			}
		}
		else{
			for(i=0;i<nx;i++){
				work[i]=a[i][j];
			}

			fftr(work,kx);

			float dkx=6.283185/nx/dx;
			float sft=6.283185/nx/2.0;
			int nxd2=nx/2;

			float ws,wc;

			work[0]=0.0;
			work[nxd2]=0.0;

			for(i=1;i<nxd2;i++){
				wc=-dkx*i*(work[i]*cos(sft*i)+work[nxd2+i]*sin(sft*i));
				ws=dkx*i*(work[nxd2+i]*cos(sft*i)-work[i]*sin(sft*i));
				work[i]=ws;
				work[nxd2+i]=wc;
			}

			fftri(work,kx);

			for(i=0;i<nx;i++){
				b[i][j]=work[i];
			}
		}
	}
}

void diffxsm(float **a,float **b,float *work,int nx,int ny,int ns,int kx,float dx){
	int i,j;
	for(j=0;j<ny;j++){
		if(j+1<=ns){
			for(i=0;i<nx;i++){
				a[i][j]=0.0;
			}
		}
		else{
			for(i=0;i<nx;i++){
				work[i]=a[i][j];
			}

			fftr(work,kx);

			float dkx=6.283185/nx/dx;
			float sft=6.283185/nx/2.0;
			int nxd2=nx/2;

			float ws,wc;

			work[0]=0.0;
			work[nxd2]=0.0;

			for(i=1;i<nxd2;i++){
				wc=-dkx*i*(work[i]*cos(sft*i)-work[nxd2+i]*sin(sft*i));
				ws=dkx*i*(work[nxd2+i]*cos(sft*i)+work[i]*sin(sft*i));
				work[i]=ws;
				work[nxd2+i]=wc;
			}

			fftri(work,kx);

			for(i=0;i<nx;i++){
				b[i][j]=work[i];
			}
		}
	}
}

void finidyx(float **a,float **dya,int nx,int ny,float dx,float dy,float dt){
	float c0=9.0/8.0;
	float c1=1.0/24.0;
	int i,j;

	for(i=0;i<nx;i++){
		dya[i][0]=1.0/dy*(c0*(a[i][1]-a[i][0])-c1*a[i][2]);
		dya[i][ny-2]=1.0/dy*(c0*(a[i][ny-1]-a[i][ny-2])+c1*a[i][ny-3]);
		dya[i][ny-1]=1.0/dy*(c0*(-a[i][ny-1])+c1*a[i][ny-2]);
	}

	for(j=1;j<ny-2;j++){
		for(i=0;i<nx;i++){
			dya[i][j]=1.0/dy*(c0*(a[i][j+1]-a[i][j])-c1*(a[i][j+2]-a[i][j-1]));
		}
	}
}

void finidyy(float **a,float **dya,int nx,int ny,float dx,float dy,float dt){
	float c0=9.0/8.0;
	float c1=1.0/24.0;
	int i,j;

	for(i=0;i<nx;i++){
		dya[i][0]=1.0/dy*(c0*a[i][0]-c1*a[i][1]);
		dya[i][1]=1.0/dy*(c0*(a[i][1]-a[i][0])-c1*a[i][1]);
		dya[i][ny-1]=1.0/dy*(c0*(a[i][ny-1]-a[i][ny-2])+c1*a[i][ny-3]);
	}

	for(j=2;j<ny-1;j++){
		for(i=0;i<nx;i++){
			dya[i][j]=1.0/dy*(c0*(a[i][j]-a[i][j-1])-c1*(a[i][j+1]-a[i][j-2]));
		}
	}
}

float dherrman(float a,float x,float x0){
	float a2=2.0*a;
	float t=x-x0;
	float td=(t+a2)/a;
	if(t<=-a2) return 0.0;
	if(t<=-a) return td/(a2*a);
	if(t<=a) return (-td+2.0)/(a2*a);
	if(t<=a2) return (td-4.0)/(a2*a);
	return 0.0;
}

float herrman(float a,float x,float x0){
	float a2=2.0*a;
	float t=x-x0;
	float td=(t+a2)/a;
	if(t<=-a2) return 0.0;
	if(t<=-a) return (0.5*td*td)/a2;
	if(t<=a) return (-0.5*td*td+2.0*td-1.0)/a2;
	if(t<=a2) return (0.5*td*td-4.0*td+8.0)/a2;
	return 0.0;
}

float fxmxz(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at,float xs,float zs){
	float x0=i0*dx+xs;
	float z0=j0*dz+zs;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=herrman(ax,x,x0);
	float fhz=-dherrman(az,z,z0);
	float fht=herrman(at,t,t0);
	return fhx*fhz*fht;
}

float fzmxz(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at,float xs,float zs){
	float x0=i0*dx+xs;
	float z0=j0*dz+zs;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=-dherrman(ax,x,x0);
	float fhz=herrman(az,z,z0);
	float fht=herrman(at,t,t0);
	// printf("***  %f %f %f %f %d %d\n",ax,x,x0,dherrman(ax,x,x0),i,j);
	return fhx*fhz*fht;
}

float fzmzz(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at,float xs,float zs){
	float x0=i0*dx+xs;
	float z0=j0*dz+zs;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=herrman(ax,x,x0);
	float fhz=-dherrman(az,z,z0);
	float fht=herrman(at,t,t0);
	return fhx*fhz*fht;
}

float fxmxx(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at,float xs,float zs){
	float x0=i0*dx+xs;
	float z0=j0*dz+zs;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=-dherrman(ax,x,x0);
	float fhz=herrman(az,z,z0);
	float fht=herrman(at,t,t0);
	return fhx*fhz*fht;
}

float fx(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at){
	float x0=i0*dx;
	float z0=j0*dz;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=herrman(ax,x,x0);
	float fhz=herrman(az,z,z0);
	float fht=herrman(at,t,t0);
	return fhx*fhz*fht;
}

float fz(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at){
	float x0=i0*dx;
	float z0=j0*dz;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=herrman(ax,x,x0);
	float fhz=herrman(az,z,z0);
	float fht=herrman(at,t,t0);
	return fhx*fhz*fht;
}

float exforce(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at){
	float x0=i0*dx;
	float z0=j0*dz;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=-dherrman(ax,x,x0);
	float fhz=herrman(az,z,z0);
	float fht=herrman(at,t,t0);
	return fhx*fhz*fht;
}

float ezforce(int i,int j,int i0,int j0,float dx,float dz,float ax,float az,float t,float t0,float at){
	float x0=i0*dx;
	float z0=j0*dz;
	float x=(i+1)*dx;
	float z=(j+1)*dz;
	float fhx=herrman(ax,x,x0);
	float fhz=-dherrman(az,z,z0);
	float fht=herrman(at,t,t0);
	return fhx*fhz*fht;
}

#endif /* PSV_H_ */

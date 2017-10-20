/*
 * const.h
 *
 *  Created on: Mar 7, 2014
 *      Author: widget
 */

#ifndef CONST_H_
#define CONST_H_

const int nx=2048,ny=1024,nx2=nx*2,ny2=ny*2;
const float dx=0.0342,dy=0.0342,dt=1.0e-3;
const int ntmax=30000,nwrite=150;
const float ax=dx,ay=dy,at=0.1/4.0,t0=at*2;
const int na=0;
const int ii0=292,jj0=292;
const float rmxx=1.0,rmxy=0.0,rmyy=-1.0,rmyx=0.0;
const float fxx=0.0,fyy=0.0,fzz=0.0;
const float dpxx=0.0,dpyy=0.0,dpzz=0.0;
const int nst=nx/4,nsskip=nx/nst;
const int nxa=20,nya=20;
const int nskip=10,ntskp=ntmax/nskip+1;

const int nsnap=60;
const float pamp=0.5,samp=2.2;
const float pampall=0.5,sampall=2.2;

const int nmax=1024;
const double pi=3.1415926535897932384626;

const int nbt=16;

#endif /* CONST_H_ */

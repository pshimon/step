/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include <math.h>
#define ZZERO 1.0e-12
#define EPS 1.0e-6

static inline double newg0(double u,double v,double q,double zz) {
	double vv,z,vz,svq,s,t,ww,qq,vq,ss,tt;
	double g1,g2,g3,g4,g5,L1,L2,L3;
	double vv2,vv1,ts,s2t2,cvts;
	vv=v*v;
	vv2=(vv+1.0);
	vv1=vv2-2.0;
	z=sqrt(zz);
	qq=q*q;
	ww=zz+qq;
	s=sqrt(u*u+ww)-u;
	t=s+2.0*u;
	vz=v*z;
	vq=v*q;
	svq=s+vq;
	g1=2.0*z*vv2*atan2(svq,vz);
	ss=s*s;
	tt=t*t;
	ts=t-s;
	cvts=v*ts;
	s2t2=0.5*(ss+tt)*vv;
	g2=cvts;
	L1=log(svq*svq+vz*vz);
	L2=log(s);
	L3=log(vv2);
	g3=L2*2.0*q ;
	g4=(L3+L2-L1)*cvts;
	g5=L1* vv1*q;
	return g1+g2+g3+g4+g5;
}

static inline double ff0(double wtt,double wbt,double wtb,double wbb,
		double qt,double qb,double Vt,double Vb,
		double zz) {
	double res;
	res=(-newg0(wtt,Vt,qt,zz)+newg0(wbt,Vt,qt,zz))/(Vt*Vt+1.0)
		+(newg0(wtb,Vb,qb,zz)-newg0(wbb,Vb,qb,zz))/(Vb*Vb+1.0);
	return res;
}
double trgint0(double x,double y,double z,
             double x0,double y0,double z0,
             double x1,double y1,double z1,
             double x2,double y2,double z2) {
	double x10,y10,z10,x20,y20,z20;
	double nx,ny,nz,kx,ky,kz,n2x,n2y,n2z;
	double d1,d2,c2,s2,u2,u1,v2,zz,u3,v3,v4;
	double Pt,pt,pb,qt,qb,Pb,ppt,ppb,wtt,wbb,wtb,wbt,Vb,Vt,res/*,norm=0.25/M_PI*/;
	x10=x1-x0;
	y10=y1-y0;
	z10=z1-z0;
	d1=sqrt(x10*x10+y10*y10+z10*z10);
	x20=x2-x0;
	y20=y2-y0;
	z20=z2-z0;
	d2=sqrt(x20*x20+y20*y20+z20*z20);
	nx=x10/d1;
	ny=y10/d1;
	nz=z10/d1;
	n2x=x20/d2;
	n2y=y20/d2;
	n2z=z20/d2;
	c2=nx*n2x+ny*n2y+nz*n2z;
	s2=sqrt(1.0-c2*c2);
	u2=d2*c2;
	u1=d1;
	v2=d2*s2;
	kx=(n2x-c2*nx)/s2;
	ky=(n2y-c2*ny)/s2;
	kz=(n2z-c2*nz)/s2;
	x0-=x;
	y0-=y;
	z0-=z;
	u3=nx*x0+ny*y0+nz*z0;
	v3=kx*x0+ky*y0+kz*z0;	
	zz=x0*x0+y0*y0+z0*z0-u3*u3-v3*v3;
	if(zz<ZZERO) zz=ZZERO;
	v4=v2+v3;
	Pt=(u2-u1)/v2;
	Pb=u2/v2;
	ppt=sqrt(1.0+Pt*Pt);
	ppb=sqrt(1.0+Pb*Pb);
	qt=(u3+u1-Pt*v3)/ppt;
	qb=(u3-Pb*v3)/ppb;
	wtt=ppt*v4+Pt*qt;
	wtb=ppb*v4+Pb*qb;
	wbt=ppt*v3+Pt*qt;
	wbb=ppb*v3+Pb*qb;
	pt=Pt/ppt;
	pb=Pb/ppb;
	Vt=sqrt((1.0+pt)/(1.0-pt));
	Vb=sqrt((1.0+pb)/(1.0-pb));  
	res=ff0(wtt,wbt,wtb,wbb,qt,qb,Vt,Vb,zz);
    return res;

}
//maple checked!!!
inline static dbl psi(dbl s,dbl u,dbl v,dbl w2) {
    dbl vv=(1.0-v*v);
    dbl dd=w2*vv-u*u;
    dbl d=sqrt(dd);
    dbl r=sqrt(s*s+w2);
    dbl b=r+v*s+u;
    dbl a=d*(r*v+s);
    dbl c=dd+u*b;
    return -s*vv+(s*vv-v*u)*log(b)-u*log(r-s)+d*atan2(a,c);
};
/* psi0(s,u,w)= psi(s,u,0,w) */
inline static dbl psi0(dbl s,dbl u,dbl w2) {
    dbl dd=w2-u*u;
    dbl d=sqrt(dd);
    dbl r=sqrt(s*s+w2);
    dbl b=r+u;
    dbl a=d*s;
    dbl c=dd+u*b;
    return -s+s*log(b)-u*log(r-s)+d*atan2(a,c);
};
/* triangle */
inline static dbl pot_patch0(dbl a,dbl p,dbl h,dbl x,dbl y,dbl z2) {
    dbl sR,sL,w2L,w2R,uL,uR,vL,vR;
    dbl a1L,a2L,a1R,a2R;
    dbl aL,bL,aR,bR;
    aL=p/h;aR=(p-a)/h;
    bL=aL*y-x;bR=aR*y-x+a;
    a2R=1.0+aR*aR;
    a1R=sqrt(a2R);
    a2L=1.0+aL*aL;
    a1L=sqrt(a2L);
    sR=aR*bR/a2R-y;
    sL=aL*bL/a2L-y;
    w2R=z2/a2R+bR*bR/(a2R*a2R);
    w2L=z2/a2L+bL*bL/(a2L*a2L);
    vR=aR/a1R;
    vL=aL/a1L;
    uR=bR/(a1R*a2R);
    uL=bL/(a1L*a2L);
    return h*log(a1R/a1L)+(psi(sR+h,uR,vR,w2R)-psi(sR,uR,vR,w2R))*a2R+(psi(sL,uL,vL,w2L)-psi(sL+h,uL,vL,w2L))*a2L;
};
/* parallelogram */
inline static dbl pot_patch1(dbl a,dbl p,dbl h,dbl x,dbl y,dbl z2) {
    dbl sR,sL,w2L,w2R,uL,uR,vL;
    dbl a1L,a2L;
    dbl aL,bL,bR;
    aL=p/h;
    bL=aL*y-x;bR=bL+a;
    a2L=1.0+aL*aL;
    a1L=sqrt(a2L);
    sR=aL*bR/a2L-y;
    sL=aL*bL/a2L-y;
    w2R=z2/a2L+bR*bR/(a2L*a2L);
    w2L=z2/a2L+bL*bL/(a2L*a2L);
    vL=aL/a1L;
    uR=bR/(a1L*a2L);
    uL=bL/(a1L*a2L);
    return (psi(sR+h,uR,vL,w2R)-psi(sR,uR,vL,w2R)+psi(sL,uL,vL,w2L)-psi(sL+h,uL,vL,w2L))*a2L;
};
/* rectangular, p=0 */
inline static dbl pot_patch2(dbl a,dbl h,dbl x,dbl y,dbl z2) {
    dbl sR,sL,w2L,w2R,uL,uR;
    dbl bL,bR;
    bL=-x;bR=-x+a;
    sR=-y;
    sL=-y;
    w2R=z2+bR*bR;
    w2L=z2+bL*bL;
    uR=bR;
    uL=bL;
    return psi0(sR+h,uR,w2R)-psi0(sR,uR,w2R)+psi0(sL,uL,w2L)-psi0(sL+h,uL,w2L);
};
inline static dbl pot_patch0(dbl a,dbl p,dbl h,dbl x,dbl y,dbl z2) {
    dbl sR,sL,w2L,w2R,uL,uR,vL,vR;
    dbl a1L,a2L,a1R,a2R;
    dbl aL,bL,aR,bR;
    aL=p/h;aR=(p-a)/h;
    bL=aL*y-x;bR=aR*y-x+a;
    a2R=1.0+aR*aR;
    a1R=sqrt(a2R);
    a2L=1.0+aL*aL;
    a1L=sqrt(a2L);
    sR=aR*bR/a2R-y;
    sL=aL*bL/a2L-y;
    w2R=z2/a2R+bR*bR/(a2R*a2R);
    w2L=z2/a2L+bL*bL/(a2L*a2L);
    vR=aR/a1R;
    vL=aL/a1L;
    uR=bR/(a1R*a2R);
    uL=bL/(a1L*a2L);
    return h*log(a1R/a1L)+(psi(sR+h,uR,vR,w2R)-psi(sR,uR,vR,w2R))*a2R+(psi(sL,uL,vL,w2L)-psi(sL+h,uL,vL,w2L))*a2L;
};





dbl pot_patch_t(v3_flt pnt,patch_flt patch) {
    v3_flt k,n,m;
    dbl a,h,p,r2;
    dbl x,y,z,zz;
    a=sqrt(v3_dot_flt(patch[1],patch[1]));
    r2=v3_dot_flt(patch[2],patch[2]);
    k[0]=patch[1][0]/a;k[1]=patch[1][1]/a;k[2]=patch[1][2]/a;
    p=patch[2][0]*k[0]+patch[2][1]*k[1]+patch[2][2]*k[2];
    h=sqrt(r2-p*p);
    m[0]=(patch[2][0]-p*k[0])/h;
    m[1]=(patch[2][1]-p*k[1])/h;
    m[2]=(patch[2][2]-p*k[2])/h; 
    n[0]=k[1]*m[2]-k[2]*m[1];n[1]=k[2]*m[0]-k[0]*m[2];n[2]=k[0]*m[1]-k[1]*m[0];
    x=v3_dot_flt(pnt,k)-v3_dot_flt(patch[0],k);
    y=v3_dot_flt(pnt,m)-v3_dot_flt(patch[0],m);
    z=v3_dot_flt(pnt,n)-v3_dot_flt(patch[0],n);
    zz=z*z;
    if(zz<ZZERO) zz=ZZERO;
    return pot_patch0(a,p,h,x,y,zz);
}
/* parallelogram */
inline static dbl pot_patch1(dbl a,dbl p,dbl h,dbl x,dbl y,dbl z2) {
    dbl sR,sL,w2L,w2R,uL,uR,vL;
    dbl a1L,a2L;
    dbl aL,bL,bR;
    aL=p/h;
    bL=aL*y-x;bR=bL+a;
    a2L=1.0+aL*aL;
    a1L=sqrt(a2L);
    sR=aL*bR/a2L-y;
    sL=aL*bL/a2L-y;
    w2R=z2/a2L+bR*bR/(a2L*a2L);
    w2L=z2/a2L+bL*bL/(a2L*a2L);
    vL=aL/a1L;
    uR=bR/(a1L*a2L);
    uL=bL/(a1L*a2L);
    return (psi(sR+h,uR,vL,w2R)-psi(sR,uR,vL,w2R)+psi(sL,uL,vL,w2L)-psi(sL+h,uL,vL,w2L))*a2L;
};
/* rectangular, p=0 */
inline static dbl pot_patch2(dbl a,dbl h,dbl x,dbl y,dbl z2) {
    dbl sR,sL,w2L,w2R,uL,uR;
    dbl bL,bR;
    bL=-x;bR=-x+a;
    sR=-y;
    sL=-y;
    w2R=z2+bR*bR;
    w2L=z2+bL*bL;
    uR=bR;
    uL=bL;
    return psi0(sR+h,uR,w2R)-psi0(sR,uR,w2R)+psi0(sL,uL,w2L)-psi0(sL+h,uL,w2L);
};
/* parallelogramm */
dbl pot_patch_p(v3_flt pnt,patch_flt patch) {
    v3_flt k,n,m;
    dbl a,h,p,r2;
    dbl x,y,z,zz;
    a=sqrt(v3_dot_flt(patch[1],patch[1]));
    r2=v3_dot_flt(patch[2],patch[2]);
    k[0]=patch[1][0]/a;k[1]=patch[1][1]/a;k[2]=patch[1][2]/a;
    p=patch[2][0]*k[0]+patch[2][1]*k[1]+patch[2][2]*k[2];
    h=sqrt(r2-p*p);
    m[0]=(patch[2][0]-p*k[0])/h;
    m[1]=(patch[2][1]-p*k[1])/h;
    m[2]=(patch[2][2]-p*k[2])/h; 
    n[0]=k[1]*m[2]-k[2]*m[1];n[1]=k[2]*m[0]-k[0]*m[2];n[2]=k[0]*m[1]-k[1]*m[0];
    x=v3_dot_flt(pnt,k)-v3_dot_flt(patch[0],k);
    y=v3_dot_flt(pnt,m)-v3_dot_flt(patch[0],m);
    z=v3_dot_flt(pnt,n)-v3_dot_flt(patch[0],n);
    zz=z*z;
    if(zz<ZZERO) zz=ZZERO;
    return pot_patch1(a,p,h,x,y,zz);
}
/* rectangular p=0*/
dbl pot_patch_r(v3_flt pnt,patch_flt patch) {
    v3_flt k,n,m;
    dbl a,h,p,r2;
    dbl x,y,z,zz;
    a=sqrt(v3_dot_flt(patch[1],patch[1]));
    r2=v3_dot_flt(patch[2],patch[2]);
    k[0]=patch[1][0]/a;k[1]=patch[1][1]/a;k[2]=patch[1][2]/a;
    p=patch[2][0]*k[0]+patch[2][1]*k[1]+patch[2][2]*k[2];
    h=sqrt(r2-p*p);
    m[0]=(patch[2][0]-p*k[0])/h;
    m[1]=(patch[2][1]-p*k[1])/h;
    m[2]=(patch[2][2]-p*k[2])/h; 
    n[0]=k[1]*m[2]-k[2]*m[1];n[1]=k[2]*m[0]-k[0]*m[2];n[2]=k[0]*m[1]-k[1]*m[0];
    x=v3_dot_flt(pnt,k)-v3_dot_flt(patch[0],k);
    y=v3_dot_flt(pnt,m)-v3_dot_flt(patch[0],m);
    z=v3_dot_flt(pnt,n)-v3_dot_flt(patch[0],n);
    zz=z*z;
    if(zz<ZZERO) zz=ZZERO;
    return pot_patch2(a,h,x,y,zz);
}



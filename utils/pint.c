/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "pint.h"
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
/* slower version of pot_patch_t, larger diff with trgint0 ???? */
dbl patch_pot1(v3_flt pnt,patch_flt patch) {
    dbl x0,x1,x2,y0,y1,y2,z0,z1,z2;
    dbl k0,k1,k2,m0,m1,m2;
    //v3_flt k,n,m;
    dbl a,h,p,r2;
    dbl x,y,zz;
    x0=pnt[0]-patch[0][0];y0=pnt[1]-patch[0][1];z0=pnt[2]-patch[0][2];
    x1=patch[1][0];y1=patch[1][1];z1=patch[1][2];
    x2=patch[2][0];y2=patch[2][1];z2=patch[2][2];
    a=sqrt(x1*x1+y1*y1+z1*z1);
    r2=x2*x2+y2*y2+z2*z2;
    k0=x1/a;k1=y1/a;k2=z1/a;
    p=x2*k0+y2*k1+z2*k2;
    h=sqrt(r2-p*p);
    m0=(x2-p*k0)/h;m1=(y2-p*k1)/h;m2=(z2-p*k2)/h; 
    x=x0*k0+y0*k1+z0*k0;
    y=x0*m0+y0*m1+z0*m0;
    zz=x0*x0+y0*y0+z0*z0-x*x-y*y;
    if(zz<ZZERO) zz=ZZERO;
    return pot_patch0(a,p,h,x,y,zz);
}
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
inline static dbl dpsidw2(dbl s,dbl u,dbl v,dbl w2) {
    dbl vv=(1.0-v*v);
    dbl dd=w2*vv-u*u;
    dbl d=sqrt(dd);
    dbl r=sqrt(s*s+w2);
    dbl b=r+v*s+u;
    dbl a=d*(r*v+s);
    dbl c=dd+u*b;
    return 1/d*atan2(a,c)*vv/2.0
	-1/d*(-d*s*vv*r*c*c-d*s*vv*r*a*a+d*s*s*vv*c*c
		+d*s*s*vv*a*a+d*v*u*r*c*c+d*v*u*r*a*a
		-d*v*u*s*c*c-d*v*u*s*a*a+u*b*d*c*c
		+u*b*d*a*a-b*a*vv*c*r*r+b*a*vv*c*r*s
		-b*d*d*d*v*c*r+b*d*d*d*v*c*s-2.0*b*a*d*d*vv*r*r+
		2.0*b*a*d*d*vv*r*s+b*a*d*d*u*r
		-b*a*d*d*u*s)/r/b/(r-s)/(c*c+a*a)/2.0;
};
inline static dbl dpsids(dbl s,dbl u,dbl v,dbl w2) {
    dbl vv=(1.0-v*v);
    dbl r=sqrt(s*s+w2);
    dbl b=r+v*s+u;
    return vv*log(b);
};
inline static dbl dpsidu(dbl s,dbl u,dbl v,dbl w2) {
    dbl vv=(1.0-v*v);
    dbl dd=w2*vv-u*u;
    dbl d=sqrt(dd);
    dbl r=sqrt(s*s+w2);
    dbl b=r+v*s+u;
    dbl a=d*(r*v+s);
    dbl c=dd+u*b;
    return -v*log(b)-log(r-s)-1/d*atan2(a,c)*u+1/d*(d*s*vv*c*c
	    +d*s*vv*a*a-d*v*u*c*c-d*v*u*a*a-b*a*u*c-b*a*d*d*r
	    -b*a*d*d*v*s)/b/(c*c+a*a);
};



     
      



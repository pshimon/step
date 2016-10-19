/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static inline double int_ijk(double x,double w,double v,double u,double aa,double bb,double cc) {
    double v1p=1.0+v;
    double v1m=1.0-v;
    double vv=v1p*v1m;
    double dd=w*w*vv-u*u;
    double d=sqrt(dd);
    double s=sqrt(x*x+w*w);
    double sp=s+x;
    double sm=s-x;
    double b=s+v*x+u;
    double a=d*(s*v+x);
    double c=dd+u*b;
    double ln1=log(sp);
    double ln2=log(sm);
    double ln3=log(b);
    double ii= 0.5*(x*s+w*w*ln1);
    double at=atan2(a,c);
    double jj= -x+(x-v*u/vv)*ln3-u/vv*ln2+d/vv*at;
    double dduu=d*d-u*u;
    double v1m2=1.0/(v1m*v1m);
    double v1p2=1.0/(v1p*v1p);
    double v4=1.0/(vv*vv);
    double kk= -0.25*x*x
	+0.25*u*(sm/v1m+sp/v1p)
	+0.25*(dduu*(1.0+v*v)*v4+2.0*x*x+w*w)*ln3
	+0.125*dduu*(ln2*v1m2+ln1*v1p2)
	+d*u*v*v4*at;
    return aa*ii+bb*jj+cc*kk;
}

double int_trg(double x,double y,double z,
             double x0,double y0,double z0,double q0,
             double xx1,double yy1,double zz1,double q1,
             double xx2,double yy2,double zz2,double q2) {

    double rr11,rr12,rr13,rr22,rr23,rr33,m,mm;
    double nrm,w1,w2,r1,r2,v1,v2,u1,u2,s1u,s2u,s1d,s2d;
    double a1,a2,b1,b2,b,c1,c2,c;
    double r3m,r1m,r2m;
    double x1=xx1-x0;
    double y1=yy1-y0;
    double z1=zz1-z0;
    double x2=xx2-x0;
    double y2=yy2-y0;
    double z2=zz2-z0;
    double x3=x-x0;
    double y3=y-y0;
    double z3=z-z0;
    rr11=x1*x1+y1*y1+z1*z1;r1=sqrt(rr11);
    rr12=x1*x2+y1*y2+z1*z2;
    rr13=x1*x3+y1*y3+z1*z3;
    rr22=x2*x2+y2*y2+z2*z2;r2=sqrt(rr22);
    rr23=x2*x3+y2*y3+z2*z3;
    rr33=x3*x3+y3*y3+z3*z3;
    mm=0.25*(rr11+rr22-2.0*rr12);
    m=sqrt(mm);
    nrm=0.25*sqrt(rr11*rr22-rr12*rr12)/mm;//S/2/m^2
    s1d=-rr13/r1;
    s2d=-rr23/r2;
    s1u=s1d+r1;
    s2u=s2d+r2;
    w1=sqrt(rr33-rr13*rr13/rr11);
    w2=sqrt(rr33-rr23*rr23/rr22);
    a1=(q1-q2)/r1;
    a2=(q1-q2)/r2;
    c=(q1+q2-2.0*q0)*m-0.25*(q1-q2)*(rr11-rr22)/m;
    b=2.0*q0*m+(q1-q2)*0.5*(rr13-rr23)/m;
    b1=(b+rr13/rr11*c)/r1;
    b2=(b+rr23/rr22*c)/r2;
    c1=c/rr11;
    c2=c/rr22;
    r1m=0.5*(rr11-rr12);
    r2m=0.5*(rr12-rr22);
    r3m=0.5*(rr13-rr23);
    v1=r1m/(m*r1);
    v2=r2m/(m*r2);
    u1=(r1m*rr13/rr11-r3m)/m;
    u2=(r2m*rr23/rr22-r3m)/m;
     return nrm*(int_ijk(s1u,w1,v1,u1,a1,b1,c1)-int_ijk(s1d,w1,v1,u1,a1,b1,c1)
	    -int_ijk(s2u,w2,v2,u2,a2,b2,c2)+int_ijk(s2d,w2,v2,u2,a2,b2,c2));
}


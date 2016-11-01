/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include <math.h>
#define ZZERO 1.0e-12
#define EPS 1.0e-6
static inline double newg(double u,double v,double q,double zz,
		double a,double b,double c) {
	double vv,z,vz,svq,s,t,ww,qq,vq,ss,tt,bqv;
	double g1,g2,g3,g4,g5,L1,L2,L3,g6;
	double vv2,vv1,bvvqq,ts,s2t2,cvts;
	vv=v*v;
	vv2=(vv+1.0);
	vv1=vv2-2.0;
	b/=vv2;
	z=sqrt(zz);
	qq=q*q;
	ww=zz+qq;
	s=sqrt(u*u+ww)-u;
	t=s+2.0*u;
	vz=v*z;
	vq=v*q;
	svq=s+vq;
	g1=2.0*c*z*vv2*atan2(svq,vz);
	ss=s*s;
	tt=t*t;
	ts=t-s;
	cvts=c*v*ts;
	s2t2=0.5*(ss+tt)*vv;
	bqv=b*q*v;
	g2=0.25*a*v*(ss-tt)+0.5*b*s2t2-bqv*vv*t-bqv*s+cvts;
	L1=log(svq*svq+vz*vz);
	L2=log(s);
	L3=log(vv2);
	bvvqq=2.0*b*vv*qq;
	g3=L2*(a*v*ww+b*ww+2.0*q*c);
	g4=(L3+L2-L1)*(b*s2t2-bqv*vv1*ts+cvts);
	g5=L1*(-0.5*(vv*vv+1.0)*ww*b+c*vv1*q);
	g6=(L1-L2)*bvvqq;
	return g1+g2+g3+g4+g5+g6;
}

static inline double ff(double wtt,double wbt,double wtb,double wbb,
		double qt,double qb,double Vt,double Vb,
		double zz,double a,double b,double c) {
	double res;
	res=(-newg(wtt,Vt,qt,zz,a,b,c)+newg(wbt,Vt,qt,zz,a,b,c))/(Vt*Vt+1.0)
		+(newg(wtb,Vb,qb,zz,a,b,c)-newg(wbb,Vb,qb,zz,a,b,c))/(Vb*Vb+1.0);
	return res;
}
double trgint1(double x,double y,double z,
             double x0,double y0,double z0,
             double x1,double y1,double z1,
             double x2,double y2,double z2,
	     double q0,double q1,double q2) {
	double x10,y10,z10,x20,y20,z20;
	double nx,ny,nz,kx,ky,kz,n2x,n2y,n2z;
	double d1,d2,c2,s2,u2,u1,v2,zz,u3,v3,v4,a,b,c;
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
	a=(q1-q0)/d1;
	b=(d1*(q2-q0)-d2*c2*(q1-q0))/(d1*v2);
	c=q0-a*u3-b*v3;
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
	res=ff(wtt,wbt,wtb,wbb,qt,qb,Vt,Vb,zz,a,b,c);
    return res;

}

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


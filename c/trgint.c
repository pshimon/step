/**************************************************************************
 * Author: Shimon Panfil                                                  *
 * http://panfil.org                                                      *
 * Copyright (c) 2007 Shimon Panfil: Industrial Physics and Simulations   *
 * http://industrialphys.com                                              * 
 *                                                                        *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       *
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    * 
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*
 *  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  *
 *  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  *
 *  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     * 
 *  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                *
 *************************************************************************/ 

#include <math.h>
#define ZZERO 1.0e-12
//#define ONEOVER4PI 0.07957747154594767

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

double trgint(double x,double y,double z,
             double x0,double y0,double z0,double q0,
             double x1,double y1,double z1,double q1,
             double x2,double y2,double z2,double q2) {
	double x00,y00,z00;
	double x10,y10,z10,x20,y20,z20;
	double nx,ny,nz,kx,ky,kz,n2x,n2y,n2z;
	double d1,d2,c2,s2,u2,u1,v2,zz,u3,v3,v4,a,b,c;
	double pt1,pt,pb,qt,qb,pb1,ppt,ppb,wtt,wbb,wtb,wbt,Vb,Vt,res;
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
/*
	x0-=x;
	y0-=y;
	z0-=z;
	
	u3=nx*x0+ny*y0+nz*z0;
	v3=kx*x0+ky*y0+kz*z0;	
	zz=x0*x0+y0*y0+z0*z0-u3*u3-v3*v3;

*/
	x00=x-x0;
	y00=y-x0;
	z00=z-x0;


	u3=-(nx*x00+ny*y00+nz*z00);
	v3=-(kx*x00+ky*y00+kz*z00);	
	zz=x00*x00+y00*y00+z00*z00-u3*u3-v3*v3;

	if(zz<ZZERO) zz=ZZERO;
	v4=v2+v3;
	a=(q1-q0)/d1;
	b=(d1*(q2-q0)-d2*c2*(q1-q0))/(d1*v2);
	c=q0-a*u3-b*v3;
	pt1=(u2-u1)/v2;
	pb1=u2/v2;
	ppt=sqrt(1.0+pt1*pt1);
	ppb=sqrt(1.0+pb1*pb1);
	qt=(u3+u1-pt1*v3)/ppt;
	qb=(u3-pb1*v3)/ppb;
	wtt=ppt*v4+pt1*qt;
	wtb=ppb*v4+pb1*qb;
	wbt=ppt*v3+pt1*qt;
	wbb=ppb*v3+pb1*qb;
	pt=pt1/ppt;
	pb=pb1/ppb;
	Vt=sqrt((1.0+pt)/(1.0-pt));
	Vb=sqrt((1.0+pb)/(1.0-pb));  
	res=ff(wtt,wbt,wtb,wbb,qt,qb,Vt,Vb,zz,a,b,c);
    return res;
}



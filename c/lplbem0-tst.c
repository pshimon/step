/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "lplbem.h"
#include "timers.h"
#include "lalg.h"
/* constant charge distribution */
int main(int argc,char * argv[]) {
    TSurf s;   
    double time_start,time_stop;
    int ret,n,t,i;
    Flt * cnt=0;
    Flt * qt0=0;
    Dbl *lm0=0;
    int * ipiv=0;
    Dbl * e=0;
    Dbl * xp=0;
    Dbl a,b,phi,z,dlt,mdlt,fct;
    
    if(argc!=2) {
	fprintf(stderr,"usage: %s surf \n",argv[0]);
	exit(1);
    }
    initTSurf(&s);
    time_start=cpuClock();
    ret=readTSurf(&s,argv[1]);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"readTSurf returns %d\n",ret);
	exit(1);
    }
    printf("reading surf takes %e s\n",time_stop-time_start);
    n=s.nv;
    t=s.nt;
    printf("nt=%d, nv=%d\n",t,n);
    cnt=ALLOC_MEM(Flt,3*t);
    time_start=cpuClock();
    ret=mkCenters(cnt,&s);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkCenters returns %d\n",ret);
	exit(1);
    }
    printf("mkCenters takes %e s\n",time_stop-time_start);
    qt0=ALLOC_MEM(Flt,t);
    time_start=cpuClock();
    ret=mkQtot0(qt0,&s);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkQtot0 returns %d\n",ret);
	exit(1);
    }
    printf("mkQtot0 takes %e s\n",time_stop-time_start);
    lm0=ALLOC_MEM(Dbl,t*t);
    time_start=cpuClock();
    ret=mkSAMat0(lm0,cnt,&s,lplGfC2);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkSAMat0 with lplGfC1 returns %d\n",ret);
	exit(1);
    }
    printf("mkSAMat0 with lplGfC1 takes %e s\n",time_stop-time_start);
    ipiv=ALLOC_MEM(int,t);
    time_start=cpuClock();
    ret=getrfDbl(t,lm0,ipiv);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"getrfDbl  returns %d\n",ret);
	exit(1);
    }
    printf("getrfDbl takes %e s\n",time_stop-time_start);
    e=ALLOC_MEM(Dbl,t);
    xp=ALLOC_MEM(Dbl,t);
    for(i=0;i<t;i++) {
	xp[i]=-cnt[3*i+2];
	e[i]=1.0;
    }
    ret=getrsDbl(t,lm0,1,xp,ipiv);
    if(ret) {
	fprintf(stderr,"getrsDbl returns %d\n",ret);
	exit(1);
    }
    ret=getrsDbl(t,lm0,1,e,ipiv);
    if(ret) {
	fprintf(stderr,"getrsDbl returns %d\n",ret);
	exit(1);
    }
    a=0.0;
    b=0.0;
    for(i=0;i<t;i++) {
	a+=qt0[i]*xp[i];
	b+=qt0[i]*e[i];
    }
    phi=a/b;
    printf("phi=%e\n",phi);
    a=0.0;
    for(i=0;i<t;i++) {
	e[i]=phi*e[i]-xp[i];
	a+=qt0[i]*e[i];
	
    }
    printf("qtot=%e\n",a);
    fct=0.75f/M_PI;
    dlt=0.0f;
    mdlt=0.0f;
    for(i=0;i<t;i++) {
	z=s.vvec[3*i+2];
	b=(e[i]-z*fct)*(e[i]-z*fct);
	dlt+=b;
	b=sqrtf(b);
	if(b>mdlt) mdlt=b;
    }
    dlt=sqrtf(dlt/t);
    printf("maxerr=%e,std=%e\n",mdlt,dlt);
    FREE_MEM(xp);
    FREE_MEM(e);
    FREE_MEM(ipiv);
    FREE_MEM(lm0);
    FREE_MEM(qt0);
    FREE_MEM(cnt);
    cleanTSurf(&s);
    return 0;
}

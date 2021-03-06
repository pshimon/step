/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "lplbem.h"
#include "timers.h"
#include "lalg.h"

Dbl xpot(Flt z, Dbl phi) {
    return -z+phi;
}
Dbl lplgfl1 (Vec3Flt dst ,Vec3Flt vrt0,Vec3Flt vrt1,Vec3Flt vrt2,Dbl q0,Dbl q1,Dbl q2); 





/* linear charge distribution */
int main(int argc,char * argv[]) {
    TSurf s;   
    double time_start,time_stop;
    int ret,n,t,i,maxi,j;
    Flt * cnt=0;
    Flt * qt1=0;
    Int * ntc=0;
    CList *tcvec=0;
    Dbl *lm1=0; 
    //Dbl *lm2=0;
    int * ipiv=0;
    Dbl * e=0;
    Dbl * ee=0;
    Dbl * xp=0;
    Dbl a,b,phi,z,dlt,mdlt,fct;
    Dbl phi0=1.0;
    char str[100];


    if(argc!=3) {
	fprintf(stderr,"usage: %s surf lbl\n",argv[0]);
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
    ntc=ALLOC_MEM(int,n);
    tcvec=ALLOC_MEM(CList,n);
    time_start=cpuClock();
    ret=getTrgCon(ntc,tcvec,&s);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"getTrgCon returns %d\n",ret);
	exit(1);
    }
    printf("getTrgCon takes %e s\n",time_stop-time_start);
    cnt=ALLOC_MEM(Flt,3*t);
    time_start=cpuClock();
    ret=mkCenters(cnt,&s);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkCenters returns %d\n",ret);
	exit(1);
    }
    printf("mkCenters takes %e s\n",time_stop-time_start);
    qt1=ALLOC_MEM(Flt,n);
    time_start=cpuClock();
    ret=mkQtot1(qt1,ntc,tcvec,&s);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkQtot1 returns %d\n",ret);
	exit(1);
    }
    printf("mkQtot1 takes %e s\n",time_stop-time_start);
 /*   lm2=ALLOC_MEM(Dbl,n*n);
    time_start=cpuClock();
    ret=mkSAMat1(lm2,ntc,tcvec,&s,lplGfL2);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkSAMat1 with lplGfL2 returns %d\n",ret);
	exit(1);
    }
    printf("mkSAMat1 with lplGfL2 takes %e s\n",time_stop-time_start);
*/
    lm1=ALLOC_MEM(Dbl,n*n);
    time_start=cpuClock();
    ret=mkSAMat1(lm1,ntc,tcvec,&s,lplGfL1);
  /*  for(i=0;i<n;i++) {
	for(j=0;j<n;j++) {
	    printf("dst=%d src=%d %le\n",j,i,lm1[j+n*i]);
	    if(lm1[j+n*i]<0.0) printf("***************************\n");
	}
    }*/
    //ret=mkLM(lm1,ntc,tcvec,&s);
    //ret=mkSAMat1(lm1,ntc,tcvec,&s,lplgfl1);
     time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkSAMat1 with lplGfL1 returns %d\n",ret);
	exit(1);
    }
    printf("mkSAMat1 with lplGfL1 takes %e s\n",time_stop-time_start);
    sprintf(str,"%s-lm1.bin",argv[2]);
    ret=write2DataBufDbl(lm1,n,n,str);

/*    dlt=0.0;
    for(i=0;i<n*n;i++)  {
	a=fabs(lm1[i]-lm2[i]);
	if(a>dlt) dlt=a;
    }
    printf("maxdiff between lm1 and lm2 is %e\n",dlt);
   FREE_MEM(lm2);
   */
    ipiv=ALLOC_MEM(int,n);
    time_start=cpuClock();
    ret=getrfDbl(n,lm1,ipiv);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"getrfDbl  returns %d\n",ret);
	exit(1);
    }
    printf("getrfDbl takes %e s\n",time_stop-time_start);
    e=ALLOC_MEM(Dbl,n);
    xp=ALLOC_MEM(Dbl,n);
    for(i=0;i<n;i++) {
	xp[i]=xpot(s.vvec[3*i+2],phi0);
	e[i]=1.0;
    }
    ret=getrsDbl(n,lm1,1,xp,ipiv);
    if(ret) {
	fprintf(stderr,"getrsDbl returns %d\n",ret);
	exit(1);
    }
    ret=getrsDbl(n,lm1,1,e,ipiv);
    if(ret) {
	fprintf(stderr,"getrsDbl returns %d\n",ret);
	exit(1);
    }
    a=0.0;
    b=0.0;
    for(i=0;i<n;i++) {
	a+=qt1[i]*xp[i];
	b+=qt1[i]*e[i];
    }
    phi=a/b;
    printf("phi=%e\n",phi);
    a=0.0;
    for(i=0;i<n;i++) {
	e[i]=phi*e[i]-xp[i];
	a+=qt1[i]*e[i];
	
    }
    printf("qtot=%e\n",a);
    ee=ALLOC_MEM(Dbl,n);
    fct=0.75/M_PI;
    for(i=0;i<n;i++) {
	z=s.vvec[3*i+2];
	ee[i]=e[i]-z*fct;
    }

    
    dlt=0.0;
    mdlt=0.0;
    maxi=0;
    for(i=0;i<n;i++) {
	b=ee[i]*ee[i];
	dlt+=b;
	b=sqrtf(b);
	if(b>mdlt) {mdlt=b;maxi=i;}
    }
    dlt=sqrt(dlt/n);
    printf("maxerr=%e at %d (%f %f %f) std=%e\n",mdlt,maxi,s.vvec[3*maxi+0],s.vvec[3*maxi+1],s.vvec[3*maxi+2],dlt);
    sprintf(str,"%s-ee1.bin",argv[2]);
    ret=write1DataBufDbl(ee,n,str);
    
    FREE_MEM(xp);
    FREE_MEM(ee);
    FREE_MEM(e);
    FREE_MEM(ipiv);
    FREE_MEM(lm1);
    FREE_MEM(qt1);
    FREE_MEM(cnt);
    FREE_MEM(tcvec);
    FREE_MEM(ntc);
    cleanTSurf(&s);
    return 0;
}

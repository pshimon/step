/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "lplbem.h"
#include "timers.h"
 

/* linear charge distribution */
int main(int argc,char * argv[]) {
    TSurf s;   
    double time_start,time_stop;
    int ret,n,t;
    Dbl * cnt=0;
    Dbl * qt1=0;
    Int * ntc=0;
    CList *tcvec=0;
    Dbl *lm1=0; 
    //Dbl *lm2=0;
    char str[100];
  

    if(argc!=3) {
	fprintf(stderr,"usage: %s surf lbl\n",argv[0]);
	exit(1);
    }
    initTSurf(&s);
  //  time_start=cpuClock();
    ret=readTSurf(&s,argv[1]);
  //  time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"readTSurf returns %d\n",ret);
	exit(1);
    }
    //printf("reading surf takes %e s\n",time_stop-time_start);
    n=s.nv;
    t=s.nt;
    printf("nt=%d, nv=%d\n",t,n);
    ntc=ALLOC_MEM(int,n);
    tcvec=ALLOC_MEM(CList,n);
 //   time_start=cpuClock();
    ret=getTrgCon(ntc,tcvec,&s);
  //   time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"getTrgCon returns %d\n",ret);
	exit(1);
    }
    //printf("getTrgCon takes %e s\n",time_stop-time_start);
    /*
    cnt=ALLOC_MEM(Dbl,3*t);
    time_start=cpuClock();
    ret=mkCenters(cnt,&s);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkCenters returns %d\n",ret);
	exit(1);
    }
    printf("mkCenters takes %e s\n",time_stop-time_start);
    */
    /*
    qt1=ALLOC_MEM(Dbl,n);
    time_start=cpuClock();
    ret=mkQtot1(qt1,ntc,tcvec,&s);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkQtot1 returns %d\n",ret);
	exit(1);
    }
    printf("mkQtot1 takes %e s\n",time_stop-time_start);
    */
    lm1=ALLOC_MEM(Dbl,n*n);
    time_start=cpuClock();
    ret=mkSAMat1(lm1,ntc,tcvec,&s,lplGfL1);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkSAMat1 with lplGfL1 returns %d\n",ret);
	exit(1);
    }
    printf("mkSAMat1 with lplGfL1 takes %e s\n",time_stop-time_start);
    sprintf(str,"%s-lm1.bin",argv[2]);
    ret=write2DataBufDbl(lm1,n,n,str);
    time_start=cpuClock();
    ret=mkSAMat1Tst(lm1,ntc,tcvec,&s); //"wrong" order of loops  as in ref
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkSAMat1 with trgint returns %d\n",ret);
	exit(1);
    }
    printf("mkSAMat1 with trgint takes %e s\n",time_stop-time_start);
    sprintf(str,"%s-lm1-ref.bin",argv[2]);
    ret=write2DataBufDbl(lm1,n,n,str);

    
    FREE_MEM(lm1);
    FREE_MEM(qt1);
    FREE_MEM(cnt);
    FREE_MEM(tcvec);
    FREE_MEM(ntc);
    cleanTSurf(&s);
    return 0;
}

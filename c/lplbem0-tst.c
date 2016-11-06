/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "lplbem.h"
#include "timers.h"
/* constant charge distribution */
int main(int argc,char * argv[]) {
    TSurf s;   
    double time_start,time_stop;
    int ret,n,t;
    Flt * cnt=0;
    Flt * qt0=0;
    Dbl *lm0=0;
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
    ret=mkSAMat0(lm0,cnt,&s,lplGfC1);
    time_stop=cpuClock();
    if(ret) {
	fprintf(stderr,"mkSAMat0 with lplGfC1 returns %d\n",ret);
	exit(1);
    }
    printf("mkSAMat0 with lplGfC1 takes %e s\n",time_stop-time_start);
    FREE_MEM(lm0);
    FREE_MEM(qt0);
    FREE_MEM(cnt);
    cleanTSurf(&s);
    return 0;
}

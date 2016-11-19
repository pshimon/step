/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "tsurf.h"

int main(int argc,char * argv[]) {
    int ret,j;
    TSurf s1,s2; 
    CList * vcvec;
    int * nvc;
    TSurf *sold,*snew,*tmp;
    char fname[32];
    initTSurf(&s1);
    initTSurf(&s2);
    sold=&s1;
    snew=&s2;
    ret=mkIcosahedron(sold);
    for(j=0;j<=5;j++) {
	sprintf(fname,"usph_%d.tsb",j);
	printf("%s %d %d\n",fname,sold->nt,sold->nv);
	printStat(sold);
	ret=writeTSurf(sold,fname);
	nvc=ALLOC_MEM(int,sold->nv);
	vcvec=ALLOC_MEM(CList,sold->nv);
	ret= getVrtCon(nvc,vcvec,sold);
	if(ret) {
	    fprintf(stderr,"get_vrt_con returns %d\n",ret);
	    exit(1);
	}
	ret=refineTSurf2(snew,sold,nvc,vcvec); 
	if(ret) {
	    fprintf(stderr,"refine2 returns %d\n",ret);
	    exit(1);
	}
	mkUnitSphere(snew);
	FREE_MEM(vcvec);
	FREE_MEM(nvc);
	tmp=sold;
	sold=snew;
	snew=tmp;
    }

    cleanTSurf(&s1);
    cleanTSurf(&s2);
    return 0;
}

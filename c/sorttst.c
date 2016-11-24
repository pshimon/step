/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "sortlib.h"
#include "timers.h"
int check(Dbl *a,int n) {
    int ret,i;
    ret=0;
    for(i=1;i<n;i++) if(a[i]<a[i-1]) {
	ret=1;
	fprintf(stderr,"a[%d]=%e <a[%d]=%e\n",i,a[i],i-1,a[i-1]);
	break;}
    return ret;
}

int main() {
    Dbl *a;
    int n,n1,i;
    int left,right,m=4,j;
    //time_t t;
    Dbl fct,t1,t0,t2,t3;
    n=100;
    n1=10;
 
    srand((unsigned) time(0));
    fct=1.0/RAND_MAX;
    for(i=0;i<m;i++) {
	a=ALLOC_MEM(Dbl,n);
	for(j=0;j<n;j++) a[j]=fct*rand();
	left=0;
	right=n-1;
	t0=cpuClock();
	sortHeapDbl(a,left,right);
//	sortShellDbl(a,left,right);
//	sortInsDbl(a,left,right);
//	sortSelDbl(a,left,right);
//	sortQuickNrDbl(a,left,right);
//	sortQuickDbl(a,left,right);
	t1=cpuClock();
	if(check(a,n)) {
	    fprintf(stderr,"array not sorted!\n");
	    exit(1);
	}
// now sort the sorted list
	t2=cpuClock();
	sortHeapDbl(a,left,right);
//	sortInsDbl(a,left,right);	
//	sortInsDbl(a,left,right);	
//	sortSelDbl(a,left,right);
//	sortQuickNrDbl(a,left,right);
//	sortQuickDbl(a,left,right);
    	t3=cpuClock();
	if(check(a,n)) {
	    fprintf(stderr,"array not sorted!\n");
	    exit(1);
	}
	printf("%d %e %e\n",n,t1-t0,t3-t2);
	n*=n1;
	FREE_MEM(a);
    }
    return 0;
}

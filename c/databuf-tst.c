/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "databuf.h"


int main(int argc,char * argv[]) {
    int s[ARRAY_MAX_RANK],r;
    int ret,i;
    float maxdiff,diff;
    char str[30];
    DataBufDbl a,b;
    initDataBufDbl(&a);
    initDataBufDbl(&b);
    for(i=0;i<ARRAY_MAX_RANK;i++) s[i]=i+2;
    for(r=1;r<=ARRAY_MAX_RANK;r++) {
	printf("rank=%d\n",r);
	makeDataBufDbl(&a,r,s);
	printf("initial arr\n");
	for(i=0;i<ARRAY_SHAPE_LENGTH;i++) printf("%d ",a.shape[i]);
	printf("\n");
	for(i=0;i<ARRAY_SHAPE_LENGTH;i++) printf("%d ",a.stride[i]);
	printf("\n");
	for(i=0;i<a.stride[r];i++) a.data[i]=0.001*i;
	sprintf(str,"DataBufDbl%d.bin",r);
	ret=writeDataBufDbl(&a,str);
	if(ret) {
	    fprintf(stderr,"writeDataBufDbl returns %d\n",ret);
	    exit(1);
	}
	ret=readDataBufDbl(&b,str);
	if(ret) {
	    fprintf(stderr,"readDataBufDbl returns %d\n",ret);
	    exit(1);
	}
	printf("final arr\n");
	for(i=0;i<ARRAY_SHAPE_LENGTH;i++) printf("%d ",b.shape[i]);
	printf("\n");
	for(i=0;i<ARRAY_SHAPE_LENGTH;i++) printf("%d ",b.stride[i]);
	printf("\n");
	maxdiff=0.0f;
	for(i=0;i<a.stride[r];i++) {
	    diff=fabs(a.data[i]-b.data[i]);
	    if (diff>maxdiff) maxdiff=diff;
	}
	printf("maxdiff %f\n",maxdiff);
    }
    printf("%f\n",AEL7(&b,1,1,1,1,1,1,1));
    return 0;
}

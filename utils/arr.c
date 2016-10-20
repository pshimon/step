/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#include "arr.h"
/* float arrays */

void initAFlt(AFlt *a) {
    int i;
    a->data=0;
    for(i=0;i<=ARRAY_MAX_RANK;i++)a->shape[i]=0;
    for(i=0;i<ARRAY_MAX_RANK-1;i++)a->stride[i]=0;
}
void cleanAFlt(AFlt *a) {
    FREE_MEM(a->data);
    initAFlt(a);
}

int makeAFlt(AFlt *a,int rank,int ind[]) {//ind[] should be of length rank at least
    int i,oldsize,newsize;
    oldsize=1;for(i=1;i<=rank;i++) oldsize*=a->shape[i];
    if(oldsize<0) return -1;//bad (uninitialized?) Array
    if((rank<1)||(rank>ARRAY_MAX_RANK)) return -2; //bad rank
    newsize=1;for(i=0;i<rank;i++) {
	if(ind[i]<1) return -3; //bad index
	newsize*=ind[i];
    }
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(Flt,newsize);
	if(a->data==0) return -4;//memory problem
    }
    for(i=0;i<=ARRAY_MAX_RANK;i++)a->shape[i]=0;
    for(i=0;i<ARRAY_MAX_RANK-1;i++)a->stride[i]=0;
    a->shape[0]=rank;
    for(i=1;i<=rank;i++) a->shape[i]=ind[i-1];
    if(1==rank) return 0;//no need to compute strides
    a->stride[0]=a->shape[1];
    for(i=1;i<rank-1;i++) a->stride[i]=a->stride[i-1]*a->shape[i+1];
    return 0;
}

int writeAFlt(AFlt *a,char *file) {
    FILE* fp;
    size_t n,r;
    int ret,i;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    r=a->shape[0];
    n=1;
    for(i=1;i<=r;i++) n*=a->shape[i];
    if(fwrite(a->shape,sizeof(int),r+1,fp)!=r+1) {ret=-2;goto abend;} //shape
    if(fwrite(a->data,sizeof(Flt),n,fp)!=n) {ret=-3;goto abend;}//data
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int readAFlt(ArrFlt *a,char *file) {
    int ret,n,i,s[ARRAY_MAX_RANK],r;
    FILE* fp;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&r,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}   //rank
    if(fread(s,sizeof(int),r,fp)!=r) {ret=-2;goto abend;}   //shape
    if(makeAFlt1(a,r,s)){ret=-3;goto abend;};
    for(i=1;i<=r;i++) n*=a->shape[i];
    if(fread(a->data,sizeof(float),n,fp)!=n) {ret=-3; goto abend;}//data
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 


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
    for(i=0;i<ARRAY_SHAPE_LENGTH;i++)a->shape[i]=1;
    for(i=0;i<ARRAY_SHAPE_LENGTH;i++)a->stride[i]=0;
}
void cleanAFlt(AFlt *a) {
    FREE_MEM(a->data);
    initAFlt(a);
}

int makeAFlt(AFlt *a,int rank,int ind[]) {//ind[] should be of length rank at least
    int i,oldsize,newsize,oldrank;
    oldrank=a->shape[0]-Flt_LBL;
    oldsize=a->stride[oldrank];
    if(oldsize<0) return -1;//bad (uninitialized?) Array
    if((rank<1)||(rank>ARRAY_MAX_RANK)) return -2; //bad rank
    a->shape[0]=rank+Flt_LBL;
    a->stride[0]=1;
    for(i=1;i<=rank;i++) {
	if(ind[i-1]<1) return -3; //bad index
	a->shape[i]=ind[i-1];
	a->stride[i]=a->stride[i-1]*a->shape[i];
    }
    newsize=a->stride[rank];
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(Flt,newsize);
	if(a->data==0) return -4;//memory problem
    }
    return 0;
}

int writeAFlt(AFlt *a,char *file) {
    FILE* fp;
    size_t n;
    int ret,rank;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    n=ARRAY_SHAPE_LENGTH;   
    if(fwrite(a->shape,sizeof(int),n,fp)!=n) {ret=-2;goto abend;} //shape
    rank=a->shape[0]-Flt_LBL;
    n=a->stride[rank];//total size
    if(fwrite(a->data,sizeof(Flt),n,fp)!=n) {ret=-3;goto abend;}//data
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int readAFlt(AFlt *a,char *file) {
    int ret,s[ARRAY_SHAPE_LENGTH],rank;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    n=ARRAY_SHAPE_LENGTH;
    if(fread(s,sizeof(int),n,fp)!=n) {ret=-2;goto abend;}   //shape
    rank=s[0]-Flt_LBL;
    if(makeAFlt(a,rank,s+1)){ret=-3;goto abend;};
    n=a->stride[rank];//total size
    if(fread(a->data,sizeof(Flt),n,fp)!=n) {ret=-4;goto abend;}//data
    ret=0;
abend:    
    fclose(fp);
    return  ret;
} 


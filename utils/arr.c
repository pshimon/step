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
    int i,oldsize,newsize;
    oldsize=a->stride[ARRAY_MAX_RANK];
    if(oldsize<0) return -1;//bad (uninitialized?) Array
    if((rank<1)||(rank>ARRAY_MAX_RANK)) return -2; //bad rank
    
    //newsize=1;
    a->shape[0]=rank;
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
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    n=ARRAY_SHAPE_LENGTH;   
    if(fwrite(a->shape,sizeof(int),n,fp)!=n) {ret=-2;goto abend;} //shape
    n=a->stride[a->shape[0]];//total size
    if(fwrite(a->data,sizeof(Flt),n,fp)!=n) {ret=-3;goto abend;}//data
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int readAFlt(AFlt *a,char *file) {
    int ret,s[ARRAY_SHAPE_LENGTH];
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    n=ARRAY_SHAPE_LENGTH;
    if(fread(s,sizeof(int),n,fp)!=n) {ret=-2;goto abend;}   //shape
    if(makeAFlt(a,s[0],s+1)){ret=-3;goto abend;};
    n=a->stride[a->shape[0]];//total size
    if(fread(a->data,sizeof(Flt),n,fp)!=n) {ret=-4;goto abend;}//data
    ret=0;
abend:    
    fclose(fp);
    return  ret;
} 


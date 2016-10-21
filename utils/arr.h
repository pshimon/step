/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef ARR_H
#define ARR_H
#include "my_cdefs.h"
#define ARRAY_MAX_RANK 7		// no need to change that	
#define ARRAY_SHAPE_LENGTH (ARRAY_MAX_RANK+1) //8 no need to change that

/* C-style index: i1=0,n1-1,i2=0,n2-1,... col-major */

#define AEL1(a,i1) (a)->data[(a)->stride[0]*(i1)]

#define AEL2(a,i1,i2) (a)->data[(a)->stride[0]*(i1)+(a)->stride[1]*(i2)]

#define AEL3(a,i1,i2,i3) (a)->data[(a)->stride[0]*(i1)+(a)->stride[1]*(i2)\
			    +(a)->stride[2]*(i3)]

#define AEL4(a,i1,i2,i3,i4) (a)->data[(a)->stride[0]*(i1)+(a)->stride[1]*(i2)\
			    +(a)->stride[2]*(i3)+(a)->stride[3]*(i4)]

#define AEL5(a,i1,i2,i3,i4,i5)	(a)->data[(a)->stride[0]*(i1)+(a)->stride[1]*(i2)\
				+(a)->stride[2]*(i3)+(a)->stride[3]*(i4)\
				+(a)->stride[4]*(i5)]

#define AEL6(a,i1,i2,i3,i4,i5,i6) (a)->data[(a)->stride[0]*(i1)+(a)->stride[1]*(i2)\
				    +(a)->stride[2]*(i3)+(a)->stride[3]*(i4)\
				    +(a)->stride[4]*(i5)+(a)->stride[5]*(i6)]

#define AEL7(a,i1,i2,i3,i4,i5,i6,i7) (a)->data[(a)->stride[0]*(i1)+(a)->stride[1]*(i2)\
					+(a)->stride[2]*(i3)+(a)->stride[3]*(i4)\
					+(a)->stride[4]*(i5)+(a)->stride[5]*(i6)\
					+(a)->stride[6]*(i7)]

typedef struct _AFlt {
    Flt *data;
    int shape[ARRAY_SHAPE_LENGTH];//shape[0]=rank
    int stride[ARRAY_SHAPE_LENGTH];//stride[0] is normally 1,stride[rank] is total size
} AFlt;

void initAFlt(AFlt *a);
void cleanAFlt(AFlt *a);
int makeAFlt(AFlt *a,int rank,int ind[]);//ind[] should be of length rank at least
int writeAFlt(AFlt *a,char *file);
int readAFlt(AFlt *a,char *file);

#endif

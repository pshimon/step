/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef ARR_H
#define ARR_H
#include "my_cdefs.h"
#define ARRAY_MAX_RANK 7

#define AEL1(a,i1) (a)->data[(i1)]
/* C-style index: i1=0,n1-1,i2=0,n2-1 col-major */
#define AEL2(a,i1,i2)	    (a)->data[(i1)\
			    +(a)->stride[0]*(i2)]

#define AEL3(a,i1,i2,i3)    (a)->data[(i1)\
			    +(a)->stride[0]*(i2)\
			    +(a)->stride[1]*(i3)]

#define AEL4(a,i1,i2,i3,i4) (a)->data[(i1)\
			    +(a)->stride[0]*(i2)\
			    +(a)->stride[1]*(i3)\
			    +(a)->stride[2]*(i4)]

#define AEL5(a,i1,i2,i3,i4,i5)	(a)->data[(i1)\
				+(a)->stride[0]*(i2)\
				+(a)->stride[1]*(i3)\
				+(a)->stride[2]*(i4)\
				+(a)->stride[3]*(i5)]

#define AEL6(a,i1,i2,i3,i4,i5,i6)   (a)->data[(i1)\
				    +(a)->stride[0]*(i2)\
				    +(a)->stride[1]*(i3)\
				    +(a)->stride[2]*(i4)\
				    +(a)->stride[3]*(i5)\
				    +(a)->stride[4]*(i6)]

#define AEL6(a,i1,i2,i3,i4,i5,i6,i7)	(a)->data[(i1)\
					+(a)->stride[0]*(i2)\
					+(a)->stride[1]*(i3)\
					+(a)->stride[2]*(i4)\
					+(a)->stride[3]*(i5)\
					+(a)->stride[4]*(i6)\
					+(a)->stride[5]*(i7)]

typedef struct _AFlt {
    Flt *data;
    int shape[ARRAY_MAX_RANK+1];//shape[0]=rank
    int stride[ARRAY_MAX_RANK-1];//first stride is always 1, not stored
} AFlt;

void initAFlt(AFlt *a);
void cleanAFlt(AFlt *a);
int makeAFlt(AFlt *a,int rank,int ind[]);//ind[] should be of length rank at least
int writeAFlt(AFlt *a,char *file);
int readAFlt(ArrFlt *a,char *file);

#endif

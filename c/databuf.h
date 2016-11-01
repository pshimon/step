/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef DATABUF_H
#define DATABUF_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

/* 
 * Mathematical constant
 */ 
# define MATH_PI           3.14159265358979323846  /* pi */
# define MATH_PI_2         1.57079632679489661923  /* pi/2 */
# define MATH_PI_4         0.78539816339744830962  /* pi/4 */
# define MATH_1_PI         0.31830988618379067154  /* 1/pi */
# define MATH_2_PI         0.63661977236758134308  /* 2/pi */
# define MATH_2PI	    (2.0*MATH_PI)          /* 2*pi */
# define MATH_4PI	    (4.0*MATH_PI)          /* 4*pi */
# define MATH_SQRT2PI	    2.506628274631 /* (2*pi)^(1/2) */
# define MATH_1_SQRT2PI	    0.39894228040143 /* (2*pi)^(-1/2) */
# define MATH_E            2.7182818284590452354   /* e */
# define MATH_LOG2E        1.4426950408889634074   /* log_2 e */
# define MATH_LOG10E       0.43429448190325182765  /* log_10 e */
# define MATH_LN2          0.69314718055994530942  /* log_e 2 */
# define MATH_LN10         2.30258509299404568402  /* log_e 10 */
# define MATH_RAD2DEG	    (180.0*MATH_1_PI)     
# define MATH_DEG2RAD	    (0.00555555555556 *MATH_PI)   /* pi/180 */  
# define MATH_HUGE_NUMBER    1.0e300  /* practically infinity */
# define MATH_TINY_NUMBER    1.0e-300 /* practically 0 */
# define MATH_LARGE_NUMBER    1.0e30  /* practically infinity */
# define MATH_SMALL_NUMBER    1.0e-30 /* practically 0 */


#define ALLOC_MEM(type,n) ( type *)malloc((size_t) (n)*sizeof( type))
#define ALLOC_MEM0(type,n) ( type *)calloc((size_t) (n),sizeof( type))
#define ALLOC_MEM1(type) ( type *) calloc((size_t) 1, sizeof( type))
#define REALLOC_MEM(pointer,type,N) ( type *) realloc((void *)(pointer),(size_t) ((N)*sizeof( type)))
#define FREE_MEM(pnt) if(pnt) {free(pnt);pnt=0;}
#define ERROR_MESSAGE(msg) fprintf(stderr,"ERROR: %s src: %s line:%d\n",msg,__FILE__,__LINE__)

#define ARRAY_MAX_RANK 7		// no need to change that	
#define ARRAY_SHAPE_LENGTH (ARRAY_MAX_RANK+1) //8 no need to change that
#define Chr_LBL	    0
#define Flt_LBL	    16
#define Dbl_LBL	    32
#define Cflt_LBL    48
#define Cdbl_LBL    64
#define Int_LBL	    80

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

/* basic types */
typedef double	Dbl;
typedef float	Flt;
typedef int	Int;
typedef char	Chr;

/* small derived types */
typedef Flt Cfl[2];
typedef Dbl Cdb[2];

typedef Flt Vec2Flt[2];
typedef Flt Vec3Flt[3];
typedef Flt Vec4Flt[4];
typedef Flt Mat3Flt[9];  /* column major matrix:
			* m[0] m[3] m[6]
			* m[1] m[4] m[7]
			* m[2] m[5] m[8]
			*/ 
typedef Flt Mat4Flt[16];  /* column major matrix:
			* m[0] m[4] m[8]  m[12]
			* m[1] m[5] m[9]  m[13]
			* m[2] m[6] m[10] m[14]
			* m[3] m[7] m[11] m[15]
			*/

/* data buffer */
typedef struct  {
    Flt *data;
    int shape[ARRAY_SHAPE_LENGTH];//shape[0]=rank+data_lbl
    int stride[ARRAY_SHAPE_LENGTH];//stride[0] is normally 1,stride[rank] is total size
} DataBufFlt;

void initDataBufFlt(DataBufFlt *a);
void cleanDataBufFlt(DataBufFlt *a);
int makeDataBufFlt(DataBufFlt *a,int rank,int ind[]);//ind[] should be of length rank at least
int writeDataBufFlt(DataBufFlt *a,char *file);
int readDataBufFlt(DataBufFlt *a,char *file);




#endif

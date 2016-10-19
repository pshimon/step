/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef MY_CDEFS_H
#define MY_CDEFS_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

#define ALLOC_MEM(type,n) ( type *)malloc((size_t) (n)*sizeof( type))
#define ALLOC_MEM0(type,n) ( type *)calloc((size_t) (n),sizeof( type))
#define ALLOC_MEM1(type) ( type *) calloc((size_t) 1, sizeof( type))
#define REALLOC_MEM(pointer,type,N) ( type *) realloc((void *)(pointer),(size_t) ((N)*sizeof( type)))
#define FREE_MEM(pnt) if(pnt) {free(pnt);pnt=0;}
#define ERROR_MESSAGE(msg) fprintf(stderr,"ERROR: %s src: %s line:%d\n",msg,__FILE__,__LINE__)

#define SQR(a) ((a)*(a))
#define SIGN(x) ((x)>=0?1:-1)
#define SWAP(a,b) { swap_tmp=a;a=b;b=swap_tmp;}
#define HPT(a,b) (sqrt(SQR(a)+SQR(b)))
#define ARG(a,b) atan2(a,b)

/* used in ogl */
typedef float t_v3[3];
typedef float t_v4[4];
typedef float t_m3[9];	/* column major matrix:
			* m[0] m[3] m[6]
			* m[1] m[4] m[7]
			* m[2] m[5] m[8]
			*/ 
typedef float t_m4[16];/* column major matrix:
			* m[0] m[4] m[8]  m[12]
			* m[1] m[5] m[9]  m[13]
			* m[2] m[6] m[10] m[14]
			* m[3] m[7] m[11] m[15]
			*/


#if 0
/* data types */
typedef float t_c32[2];
typedef double t_c64[2];

typedef float Cfl[2];
typedef double Cdb[2];
typedef unsigned int t_u32;
typedef double	t_f64;
typedef float	t_f32;
typedef char	t_i08;
typedef unsigned char t_u08;
typedef int  t_i32;
typedef unsigned long long t_u64;

/* high resolution timer */
double cpu_clock();
#define STRLENGTH 80
#define TINY_F64 1.0e-20
#endif
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

#endif



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
# define MATH_HUGE_NUMBER_DBL    1.0e300  /* practically infinity */
# define MATH_TINY_NUMBER_DBL    1.0e-300 /* practically 0 */
# define MATH_HUGE_NUMBER_FLT    1.0e30  /* practically infinity */
# define MATH_TINY_NUMBER_FLT    1.0e-30 /* practically 0 */
# define MATH_LARGE_NUMBER_DBL    1.0e200  /* practically infinity */
# define MATH_SMALL_NUMBER_DBL    1.0e-200 /* practically 0 */
# define MATH_LARGE_NUMBER_FLT    1.0e20  /* practically infinity */
# define MATH_SMALL_NUMBER_FLT    1.0e-20 /* practically 0 */


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

#define indf3(i,j) ((i)+3*(j))
#define indf4(i,j) ((i)+4*(j))


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



inline static Flt dot3Flt(Vec3Flt u,Vec3Flt v) {/* u.v */
    return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]);
}
inline static Flt dot4Flt(Vec4Flt u, Vec4Flt v) {/* u.v */
    return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]+u[3]*v[3]);
}

inline static Flt norm3Flt(Vec3Flt u) {/* |u| */
    return sqrtf(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
}

inline static Flt norm4Flt(Vec4Flt u) {/* |u| */
    return sqrtf(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
}
inline static Flt dist3Flt(Vec3Flt u,Vec3Flt v) {/* |u-v|*/
    Vec3Flt w;
    w[0]=u[0]-v[0]; w[1]=u[1]-v[1]; w[2]=u[2]-v[2];
    return norm3Flt(w);
}

inline static void cross3Flt(Vec3Flt w,Vec3Flt u,Vec3Flt v) {/* w=uxv */
    w[0]=u[1]*v[2]-u[2]*v[1];
    w[1]=u[2]*v[0]-u[0]*v[2];
    w[2]=u[0]*v[1]-u[1]*v[0];
}

inline static void vxm3Flt(Vec3Flt v,Mat3Flt m,Vec3Flt u) {/* v=u*m */
    int j,i;
    for(i=0;i<3;i++) {
	v[i]=0.0f;
	for(j=0;j<3;j++) v[i]+=u[j]*m[indf3(j,i)];
    }
}
inline static void mxv3Flt(Vec3Flt v,Mat3Flt m,Vec3Flt u) {/* v=m*u */
    int j,i;
    for(i=0;i<3;i++) {
	v[i]=0.0f;
	for(j=0;j<3;j++) v[i]+=u[j]*m[indf3(i,j)];
    }
}
inline static void mxm3Flt(Mat3Flt c,Mat3Flt a,Mat3Flt b) {/* c=a*b */
    int j,i,k;
    for(i=0;i<3*3;i++) c[i]=0.0f; 
    for(j=0;j<3;j++) 
	for(k=0;k<3;k++)
	    for(i=0;i<3;i++)
		c[indf3(i,j)]+=a[indf3(i,k)]*b[indf3(k,j)];
}
/* a in rads */
inline static void rotm3Flt(Mat3Flt m,Vec3Flt v,Flt a) {
    int i,j;
    Vec3Flt u;
    Flt c,s;
    c=cosf(a);
    s=sinf(a);
    for(j=0;j<3;j++) u[j]=v[j]/norm3Flt(v);
    for(j=0;j<3;j++)
	for(i=0;i<3;i++)
	    m[indf3(i,j)]=u[i]*u[j]*(1.0f-c);
    for(j=0;j<3;j++) m[indf3(i,i)]+=c;
    m[indf3(1,0)]+= s*u[2];
    m[indf3(2,0)]+=-s*u[1];
    m[indf3(0,1)]+=-s*u[2];
    m[indf3(2,1)]+= s*u[0];
    m[indf3(0,2)]+= s*u[1];
    m[indf3(1,2)]+=-s*u[0];
}

inline static void vxm4Flt(Vec4Flt v,Mat4Flt m,Vec4Flt u) {/* v=u*m */
    int j,i;
    for(i=0;i<4;i++) {
	v[i]=0.0f;
	for(j=0;j<4;j++) v[i]+=u[j]*m[indf4(j,i)];
    }
}
inline static void mxv4Flt(Vec4Flt v,Mat4Flt m,Vec4Flt u) {/* v=m*u */
    int j,i;
    for(i=0;i<4;i++) {
	v[i]=0.0f;
	for(j=0;j<4;j++) v[i]+=u[j]*m[indf4(i,j)];
    }
}
inline static void mxm4Flt(Mat4Flt c,Mat4Flt a,Mat4Flt b) {/* c=a*b */
    int j,i,k;
    for(i=0;i<4*4;i++) c[i]=0.0f; 
    for(j=0;j<4;j++) 
	for(k=0;k<4;k++)
	    for(i=0;i<4;i++)
		c[indf4(i,j)]+=a[indf4(i,k)]*b[indf4(k,j)];
}

inline static void setv4Flt(Vec4Flt v4,Vec3Flt v3,Flt a) {/* v4=(v3,a) */
    int i;
    for(i=0;i<3;i++) v4[i]=v3[i];
    v4[3]=a;
}
/* a in rads */
inline static void rotm4Flt(Mat4Flt m,Vec3Flt v,Flt a) {
    int i,j;
    Vec3Flt u;
    Flt c,s;
    c=cosf(a);
    s=sinf(a);  
    for(j=0;j<3;j++) u[j]=v[j]/
	norm3Flt(v);
    for(j=0;j<3;j++)
	for(i=0;i<3;i++)
	    m[indf4(i,j)]=u[i]*u[j]*(1.0f-c);
    for(j=0;j<3;j++) m[indf4(j,j)]+=c;
    m[indf4(1,0)]+= s*u[2];
    m[indf4(2,0)]+=-s*u[1];
    m[indf4(0,1)]+=-s*u[2];
    m[indf4(2,1)]+= s*u[0];
    m[indf4(0,2)]+= s*u[1];
    m[indf4(1,2)]+=-s*u[0];
    for(j=0;j<3;j++) {
	m[indf4(j,3)]=0.0f;
	m[indf4(3,j)]=0.0f;
    }
    m[indf4(3,3)]=1.0f;
}
inline static void translate4Flt(Mat4Flt m,Vec3Flt v) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    for(j=0;j<4;j++) m[indf4(j,j)]=1.0f;
    for(j=0;j<3;j++) m[indf4(j,3)]=v[j];
}
inline static void scale4Flt(Mat4Flt m,Vec3Flt v) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    for(j=0;j<3;j++) m[indf4(j,j)]=v[j];
    m[indf4(3,3)]=1.0f;
}
inline static void frustumFlt(Mat4Flt m,Flt l,Flt r,Flt b,Flt t,Flt n,Flt f) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[indf4(0,0)]=2.0f*n/(r-l);
    m[indf4(1,1)]=2.0f*n/(t-b);
    m[indf4(0,2)]=(r+l)/(r-l);
    m[indf4(1,2)]=(t+b)/(t-b);
    m[indf4(2,2)]=-(f+n)/(f-n);
    m[indf4(3,2)]=-1.0f;
    m[indf4(2,3)]=-2.0f*f*n/(f-n);
}  
inline static void frustum_invFlt(Mat4Flt m,Flt l,Flt r,Flt b,Flt t,Flt n,Flt f) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[indf4(0,0)]=(r-l)/(2.0f*n);
    m[indf4(1,1)]=(t-b)/(2.0f*n);
    m[indf4(3,2)]=-(f-n)/(2.0f*f*n);
    m[indf4(0,3)]=(r+l)/(2.0f*n);
    m[indf4(1,3)]=(t+b)/(2.0f*n);
    m[indf4(2,3)]=-1.0f;
    m[indf4(3,3)]=(f+n)/(2.0f*f*n);
} 

inline static void orthoFlt(Mat4Flt m,Flt l,Flt r,Flt b,Flt t,Flt n,Flt f) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[indf4(0,0)]=2.0f/(r-l);
    m[indf4(1,1)]=2.0f/(t-b);
    m[indf4(2,2)]=-2.0/(f-n);
    m[indf4(0,3)]=-(r+l)/(r-l);
    m[indf4(1,3)]=-(t+b)/(t-b);
    m[indf4(2,3)]=-(f+n)/(f-n);
    m[indf4(3,3)]=1.0;
}  
inline static void ortho_invFlt(Mat4Flt m,Flt l,Flt r,Flt b,Flt t,Flt n,Flt f) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[indf4(0,0)]=(r-l)/2.0f;
    m[indf4(1,1)]=(t-b)/2.0f;
    m[indf4(3,2)]=-(f-n)/2.0f;
    m[indf4(0,3)]=(r+l)/2.0f;
    m[indf4(1,3)]=(t+b)/2.0f;
    m[indf4(2,3)]=(f+n)/2.0f;
    m[indf4(3,3)]=1.0f;
} 
inline static void proj_matFlt(Mat4Flt m,Flt fovy,Flt aspect_ratio, Flt near_plane,Flt far_plane) {
    Flt y_scale=1.0f/tan(M_PI/180*fovy/2);
    Flt x_scale=y_scale / aspect_ratio;
    Flt frustum_length = far_plane - near_plane;
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[0] = x_scale;
    m[5] = y_scale;	
    m[10] = -((far_plane + near_plane) / frustum_length);
    m[11] = -1;
    m[14] = -((2 * near_plane * far_plane) / frustum_length);
}
inline static void lookatFlt(Mat4Flt m,Vec3Flt eye,Vec3Flt cnt,Vec3Flt up) { 
    Vec3Flt f,s,u;
    int i;
    Flt a;
    for(i=0;i<3;i++) f[i]=cnt[i]-eye[i];
    a=norm3Flt(f);
    for(i=0;i<3;i++) f[i]=f[i]/a;
    cross3Flt(s,f,up);
    a=norm3Flt(s);
    for(i=0;i<3;i++) s[i]=s[i]/a;
    cross3Flt(u,s,f);
    for(i=0;i<16;i++) m[i]=0.0f;
    for(i=0;i<3;i++) {
	m[indf4(i,0)]=s[i];
	m[indf4(i,1)]=u[i];
	m[indf4(i,2)]=-f[i];
    }
    m[indf4(3,0)] =-dot3Flt(s, eye);
    m[indf4(3,1)] =-dot3Flt(u, eye);
    m[indf4(3,2)] =-dot3Flt(f, eye);
}

inline static void m3_from_m4Flt(Mat3Flt m1,Mat4Flt m2) {
    int i,j;
    for(i=0;i<3;i++)
	for(j=0;j<3;j++)
	    m1[indf3(i,j)]=m2[indf4(i,j)];
}



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

/* doubles */

typedef Dbl Vec2Dbl[2];
typedef Dbl Vec3Dbl[3];
typedef Dbl Vec4Dbl[4];
typedef Dbl Mat3Dbl[9];  /* column major matrix:
			* m[0] m[3] m[6]
			* m[1] m[4] m[7]
			* m[2] m[5] m[8]
			*/ 
typedef Dbl Mat4Dbl[16];  /* column major matrix:
			* m[0] m[4] m[8]  m[12]
			* m[1] m[5] m[9]  m[13]
			* m[2] m[6] m[10] m[14]
			* m[3] m[7] m[11] m[15]
			*/

inline static Dbl dot3Dbl(Vec3Dbl u,Vec3Dbl v) {/* u.v */
    return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]);
}
inline static Dbl dot4Dbl(Vec4Dbl u, Vec4Dbl v) {/* u.v */
    return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]+u[3]*v[3]);
}

inline static Dbl norm3Dbl(Vec3Dbl u) {/* |u| */
    return sqrtf(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
}

inline static Dbl norm4Dbl(Vec4Dbl u) {/* |u| */
    return sqrtf(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
}
inline static Dbl dist3Dbl(Vec3Dbl u,Vec3Dbl v) {/* |u-v|*/
    Vec3Dbl w;
    w[0]=u[0]-v[0]; w[1]=u[1]-v[1]; w[2]=u[2]-v[2];
    return norm3Dbl(w);
}

inline static void cross3Dbl(Vec3Dbl w,Vec3Dbl u,Vec3Dbl v) {/* w=uxv */
    w[0]=u[1]*v[2]-u[2]*v[1];
    w[1]=u[2]*v[0]-u[0]*v[2];
    w[2]=u[0]*v[1]-u[1]*v[0];
}

inline static void vxm3Dbl(Vec3Dbl v,Mat3Dbl m,Vec3Dbl u) {/* v=u*m */
    int j,i;
    for(i=0;i<3;i++) {
	v[i]=0.0f;
	for(j=0;j<3;j++) v[i]+=u[j]*m[indf3(j,i)];
    }
}
inline static void mxv3Dbl(Vec3Dbl v,Mat3Dbl m,Vec3Dbl u) {/* v=m*u */
    int j,i;
    for(i=0;i<3;i++) {
	v[i]=0.0f;
	for(j=0;j<3;j++) v[i]+=u[j]*m[indf3(i,j)];
    }
}
inline static void mxm3Dbl(Mat3Dbl c,Mat3Dbl a,Mat3Dbl b) {/* c=a*b */
    int j,i,k;
    for(i=0;i<3*3;i++) c[i]=0.0f; 
    for(j=0;j<3;j++) 
	for(k=0;k<3;k++)
	    for(i=0;i<3;i++)
		c[indf3(i,j)]+=a[indf3(i,k)]*b[indf3(k,j)];
}
/* a in rads */
inline static void rotm3Dbl(Mat3Dbl m,Vec3Dbl v,Dbl a) {
    int i,j;
    Vec3Dbl u;
    Dbl c,s;
    c=cosf(a);
    s=sinf(a);
    for(j=0;j<3;j++) u[j]=v[j]/norm3Dbl(v);
    for(j=0;j<3;j++)
	for(i=0;i<3;i++)
	    m[indf3(i,j)]=u[i]*u[j]*(1.0f-c);
    for(j=0;j<3;j++) m[indf3(i,i)]+=c;
    m[indf3(1,0)]+= s*u[2];
    m[indf3(2,0)]+=-s*u[1];
    m[indf3(0,1)]+=-s*u[2];
    m[indf3(2,1)]+= s*u[0];
    m[indf3(0,2)]+= s*u[1];
    m[indf3(1,2)]+=-s*u[0];
}

inline static void vxm4Dbl(Vec4Dbl v,Mat4Dbl m,Vec4Dbl u) {/* v=u*m */
    int j,i;
    for(i=0;i<4;i++) {
	v[i]=0.0f;
	for(j=0;j<4;j++) v[i]+=u[j]*m[indf4(j,i)];
    }
}
inline static void mxv4Dbl(Vec4Dbl v,Mat4Dbl m,Vec4Dbl u) {/* v=m*u */
    int j,i;
    for(i=0;i<4;i++) {
	v[i]=0.0f;
	for(j=0;j<4;j++) v[i]+=u[j]*m[indf4(i,j)];
    }
}
inline static void mxm4Dbl(Mat4Dbl c,Mat4Dbl a,Mat4Dbl b) {/* c=a*b */
    int j,i,k;
    for(i=0;i<4*4;i++) c[i]=0.0f; 
    for(j=0;j<4;j++) 
	for(k=0;k<4;k++)
	    for(i=0;i<4;i++)
		c[indf4(i,j)]+=a[indf4(i,k)]*b[indf4(k,j)];
}

inline static void setv4Dbl(Vec4Dbl v4,Vec3Dbl v3,Dbl a) {/* v4=(v3,a) */
    int i;
    for(i=0;i<3;i++) v4[i]=v3[i];
    v4[3]=a;
}
/* a in rads */
inline static void rotm4Dbl(Mat4Dbl m,Vec3Dbl v,Dbl a) {
    int i,j;
    Vec3Dbl u;
    Dbl c,s;
    c=cosf(a);
    s=sinf(a);  
    for(j=0;j<3;j++) u[j]=v[j]/
	norm3Dbl(v);
    for(j=0;j<3;j++)
	for(i=0;i<3;i++)
	    m[indf4(i,j)]=u[i]*u[j]*(1.0f-c);
    for(j=0;j<3;j++) m[indf4(j,j)]+=c;
    m[indf4(1,0)]+= s*u[2];
    m[indf4(2,0)]+=-s*u[1];
    m[indf4(0,1)]+=-s*u[2];
    m[indf4(2,1)]+= s*u[0];
    m[indf4(0,2)]+= s*u[1];
    m[indf4(1,2)]+=-s*u[0];
    for(j=0;j<3;j++) {
	m[indf4(j,3)]=0.0f;
	m[indf4(3,j)]=0.0f;
    }
    m[indf4(3,3)]=1.0f;
}
inline static void translate4Dbl(Mat4Dbl m,Vec3Dbl v) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    for(j=0;j<4;j++) m[indf4(j,j)]=1.0f;
    for(j=0;j<3;j++) m[indf4(j,3)]=v[j];
}
inline static void scale4Dbl(Mat4Dbl m,Vec3Dbl v) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    for(j=0;j<3;j++) m[indf4(j,j)]=v[j];
    m[indf4(3,3)]=1.0f;
}
inline static void frustumDbl(Mat4Dbl m,Dbl l,Dbl r,Dbl b,Dbl t,Dbl n,Dbl f) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[indf4(0,0)]=2.0f*n/(r-l);
    m[indf4(1,1)]=2.0f*n/(t-b);
    m[indf4(0,2)]=(r+l)/(r-l);
    m[indf4(1,2)]=(t+b)/(t-b);
    m[indf4(2,2)]=-(f+n)/(f-n);
    m[indf4(3,2)]=-1.0f;
    m[indf4(2,3)]=-2.0f*f*n/(f-n);
}  
inline static void frustum_invDbl(Mat4Dbl m,Dbl l,Dbl r,Dbl b,Dbl t,Dbl n,Dbl f) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[indf4(0,0)]=(r-l)/(2.0f*n);
    m[indf4(1,1)]=(t-b)/(2.0f*n);
    m[indf4(3,2)]=-(f-n)/(2.0f*f*n);
    m[indf4(0,3)]=(r+l)/(2.0f*n);
    m[indf4(1,3)]=(t+b)/(2.0f*n);
    m[indf4(2,3)]=-1.0f;
    m[indf4(3,3)]=(f+n)/(2.0f*f*n);
} 

inline static void orthoDbl(Mat4Dbl m,Dbl l,Dbl r,Dbl b,Dbl t,Dbl n,Dbl f) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[indf4(0,0)]=2.0f/(r-l);
    m[indf4(1,1)]=2.0f/(t-b);
    m[indf4(2,2)]=-2.0/(f-n);
    m[indf4(0,3)]=-(r+l)/(r-l);
    m[indf4(1,3)]=-(t+b)/(t-b);
    m[indf4(2,3)]=-(f+n)/(f-n);
    m[indf4(3,3)]=1.0;
}  
inline static void ortho_invDbl(Mat4Dbl m,Dbl l,Dbl r,Dbl b,Dbl t,Dbl n,Dbl f) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[indf4(0,0)]=(r-l)/2.0f;
    m[indf4(1,1)]=(t-b)/2.0f;
    m[indf4(3,2)]=-(f-n)/2.0f;
    m[indf4(0,3)]=(r+l)/2.0f;
    m[indf4(1,3)]=(t+b)/2.0f;
    m[indf4(2,3)]=(f+n)/2.0f;
    m[indf4(3,3)]=1.0f;
} 
inline static void proj_matDbl(Mat4Dbl m,Dbl fovy,Dbl aspect_ratio, Dbl near_plane,Dbl far_plane) {
    Dbl y_scale=1.0f/tan(M_PI/180*fovy/2);
    Dbl x_scale=y_scale / aspect_ratio;
    Dbl frustum_length = far_plane - near_plane;
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[0] = x_scale;
    m[5] = y_scale;	
    m[10] = -((far_plane + near_plane) / frustum_length);
    m[11] = -1;
    m[14] = -((2 * near_plane * far_plane) / frustum_length);
}
inline static void lookatDbl(Mat4Dbl m,Vec3Dbl eye,Vec3Dbl cnt,Vec3Dbl up) { 
    Vec3Dbl f,s,u;
    int i;
    Dbl a;
    for(i=0;i<3;i++) f[i]=cnt[i]-eye[i];
    a=norm3Dbl(f);
    for(i=0;i<3;i++) f[i]=f[i]/a;
    cross3Dbl(s,f,up);
    a=norm3Dbl(s);
    for(i=0;i<3;i++) s[i]=s[i]/a;
    cross3Dbl(u,s,f);
    for(i=0;i<16;i++) m[i]=0.0f;
    for(i=0;i<3;i++) {
	m[indf4(i,0)]=s[i];
	m[indf4(i,1)]=u[i];
	m[indf4(i,2)]=-f[i];
    }
    m[indf4(3,0)] =-dot3Dbl(s, eye);
    m[indf4(3,1)] =-dot3Dbl(u, eye);
    m[indf4(3,2)] =-dot3Dbl(f, eye);
}

inline static void m3_from_m4Dbl(Mat3Dbl m1,Mat4Dbl m2) {
    int i,j;
    for(i=0;i<3;i++)
	for(j=0;j<3;j++)
	    m1[indf3(i,j)]=m2[indf4(i,j)];
}




/* data buffer */
typedef struct  {
    Dbl *data;
    int shape[ARRAY_SHAPE_LENGTH];//shape[0]=rank+data_lbl
    int stride[ARRAY_SHAPE_LENGTH];//stride[0] is normally 1,stride[rank] is total size
} DataBufDbl;

void initDataBufDbl(DataBufDbl *a);
void cleanDataBufDbl(DataBufDbl *a);
int makeDataBufDbl(DataBufDbl *a,int rank,int ind[]);//ind[] should be of length rank at least
int writeDataBufDbl(DataBufDbl *a,char *file);
int readDataBufDbl(DataBufDbl *a,char *file);
int write1DataBufDbl(Dbl *a,int n1,char * fname);
int write2DataBufDbl(Dbl *a,int n1,int n2,char * fname);

#endif

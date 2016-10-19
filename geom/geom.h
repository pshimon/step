/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef GEOM_H
#define GEOM_H
#include "my_cdefs.h"
#define indf3(i,j) ((i)+3*(j))
#define indf4(i,j) ((i)+4*(j))


inline static float dot3(float * u,float * v) {/* u.v */
    return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]);
}
inline static float dot4(float * u, float *v) {/* u.v */
    return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]+u[3]*v[3]);
}

inline static float norm3(float * u) {/* |u| */
    return sqrtf(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
}

inline static float norm4(float * u) {/* |u| */
    return sqrtf(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
}
inline static float dist3(float *u,float  *v) {/* |u-v|*/
    t_v3 w;
    w[0]=u[0]-v[0]; w[1]=u[1]-v[1]; w[2]=u[2]-v[2];
    return norm3(w);
}

inline static void cross3(float * w,float * u,float * v) {/* w=uxv */
    w[0]=u[1]*v[2]-u[2]*v[1];
    w[1]=u[2]*v[0]-u[0]*v[2];
    w[2]=u[0]*v[1]-u[1]*v[0];
}

inline static void vxm3(float * v,float * m,float * u) {/* v=u*m */
    int j,i;
    for(i=0;i<3;i++) {
	v[i]=0.0f;
	for(j=0;j<3;j++) v[i]+=u[j]*m[indf3(j,i)];
    }
}
inline static void mxv3(float * v,float * m,float * u) {/* v=m*u */
    int j,i;
    for(i=0;i<3;i++) {
	v[i]=0.0f;
	for(j=0;j<3;j++) v[i]+=u[j]*m[indf3(i,j)];
    }
}
inline static void mxm3(float * c,float * a,float * b) {/* c=a*b */
    int j,i,k;
    for(i=0;i<3*3;i++) c[i]=0.0f; 
    for(j=0;j<3;j++) 
	for(k=0;k<3;k++)
	    for(i=0;i<3;i++)
		c[indf3(i,j)]+=a[indf3(i,k)]*b[indf3(k,j)];
}
/* a in rads */
inline static void rotm3(float * m,float * v,float a) {
    int i,j;
    t_v3 u;
    float c,s;
    c=cosf(a);
    s=sinf(a);
    for(j=0;j<3;j++) u[j]=v[j]/norm3(v);
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

inline static void vxm4(float * v,float * m,float * u) {/* v=u*m */
    int j,i;
    for(i=0;i<4;i++) {
	v[i]=0.0f;
	for(j=0;j<4;j++) v[i]+=u[j]*m[indf4(j,i)];
    }
}
inline static void mxv4(float *v,float * m,float * u) {/* v=m*u */
    int j,i;
    for(i=0;i<4;i++) {
	v[i]=0.0f;
	for(j=0;j<4;j++) v[i]+=u[j]*m[indf4(i,j)];
    }
}
inline static void mxm4(float * c,float * a,float * b) {/* c=a*b */
    int j,i,k;
    for(i=0;i<4*4;i++) c[i]=0.0f; 
    for(j=0;j<4;j++) 
	for(k=0;k<4;k++)
	    for(i=0;i<4;i++)
		c[indf4(i,j)]+=a[indf4(i,k)]*b[indf4(k,j)];
}

inline static void setv4(float * v4,float * v3,float a) {/* v4=(v3,a) */
    int i;
    for(i=0;i<3;i++) v4[i]=v3[i];
    v4[3]=a;
}
/* a in rads */
inline static void rotm4(float * m,float * v,float a) {
    int i,j;
    t_v3 u;
    float c,s;
    c=cosf(a);
    s=sinf(a);  
    for(j=0;j<3;j++) u[j]=v[j]/norm3(v);
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
inline static void translate4(float * m,float * v) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    for(j=0;j<4;j++) m[indf4(j,j)]=1.0f;
    for(j=0;j<3;j++) m[indf4(j,3)]=v[j];
}
inline static void scale4(float * m,float * v) {
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    for(j=0;j<3;j++) m[indf4(j,j)]=v[j];
    m[indf4(3,3)]=1.0f;
}
inline static void frustum(float * m,float l,float r,float b,float t,float n,float f) {
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
inline static void frustum_inv(float * m,float l,float r,float b,float t,float n,float f) {
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

inline static void ortho(float * m,float l,float r,float b,float t,float n,float f) {
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
inline static void ortho_inv(float * m,float l,float r,float b,float t,float n,float f) {
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
inline static void proj_mat(float * m,float fovy,float aspect_ratio, float near_plane,float far_plane) {
    float y_scale=1.0f/tan(M_PI/180*fovy/2);
    float x_scale=y_scale / aspect_ratio;
    float frustum_length = far_plane - near_plane;
    int j;
    for(j=0;j<16;j++) m[j]=0.0f;
    m[0] = x_scale;
    m[5] = y_scale;	
    m[10] = -((far_plane + near_plane) / frustum_length);
    m[11] = -1;
    m[14] = -((2 * near_plane * far_plane) / frustum_length);
}
inline static void lookat(float * m,t_v3 eye,t_v3 cnt,t_v3 up) { 
    t_v3 f,s,u;
    int i;
    float a;
    for(i=0;i<3;i++) f[i]=cnt[i]-eye[i];
    a=norm3(f);
    for(i=0;i<3;i++) f[i]=f[i]/a;
    cross3(s,f,up);
    a=norm3(s);
    for(i=0;i<3;i++) s[i]=s[i]/a;
    cross3(u,s,f);
    for(i=0;i<16;i++) m[i]=0.0f;
    for(i=0;i<3;i++) {
	m[indf4(i,0)]=s[i];
	m[indf4(i,1)]=u[i];
	m[indf4(i,2)]=-f[i];
    }
    m[indf4(3,0)] =-dot3(s, eye);
    m[indf4(3,1)] =-dot3(u, eye);
    m[indf4(3,2)] =-dot3(f, eye);
}

inline static void m3_from_m4(t_m3 m1,t_m4 m2) {
    int i,j;
    for(i=0;i<3;i++)
	for(j=0;j<3;j++)
	    m1[indf3(i,j)]=m2[indf4(i,j)];
}
#endif 

/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef SURF_H
#define SURF_H
#include "geom.h"
/* if distance between points is less than MIN_DIST 
 * they are considered coinciding
 */
#define MIN_DIST 1.0e-8
/* if norm of vector is less than MIN_NORM it is zero */
#define MIN_NORM 1.0e-12 
/* maximum number of connected vertices */
#define MAX_CONNECT 8
typedef int CList[MAX_CONNECT];
typedef struct {
    int * tvec;/*trg vector(3*nt);*/
    float * vvec;/* vertices vec(3*nv)*/
    float * nvec;/* normal vec(3*nv)*/
    int nt; /* number of triangles */
    int nv; /* number of vertices */
} TSurf;

/* unintialized variable may contain garbage */
int initTSurf(TSurf* s);
/* to avoid memory leaks */
void cleanTSurf(TSurf* s); 
/* creates empty surface */
int makeTSurf(TSurf* s,int T,int N);
int cpTSurf(TSurf* dst,TSurf*src);
/* binary */
int writeTSurf(TSurf *s,char * fname);
int readTSurf(TSurf *s,char * fname);
int printTSurf(TSurf *s,char * fname);
int mkTetrahedron(TSurf *s);
int mkHexahedron(TSurf *s);
int mkOctahedron(TSurf *s);
int mkDodecahedron(TSurf *s);
int mkIcosahedron(TSurf *s);

/* rerturns area of trg */
float trgNorm(Vec3Flt w,TSurf *s,int t);
/* array lst must be of  MAX_CONNECT length at least*/
//int get_trgs(int * lst,TSurf *s,int v);
//int get_trg_pair(int * first,int * second,int * lst,int nc,TSurf *s,int v);
//int remove_repeated_vertices(TSurf *s, int vstart);
int getTrgCon(int *ntc,CList *tcvec,TSurf *s);
int getVrtCon(int *nvc,CList *vcvec,TSurf *s);
int getVrtBetw(int *nvc,CList *vcvec,TSurf *s,int k1,int k2);
int refine2(TSurf *s,TSurf *sold,int *nvc,CList *vcvec);
int checkNormals(float * d,TSurf *s);
void mkUnitSphere(TSurf *s);
#endif

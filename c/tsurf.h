/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef SURF_H
#define SURF_H
#include "databuf.h"

#define MIN_DIST 1.0e-8 /* if distance between points is less 
			   than MIN_DIST they are considered coinciding */

#define MIN_NORM 1.0e-12/* if norm of vector is less than 
			   MIN_NORM it is zero */
 
#define MAX_CONNECT 8 /* maximum number of connected vertices */

typedef int CList[MAX_CONNECT];

typedef struct {
    int * tvec;/*trg vector(3*nt);*/
    Dbl * vvec;/* vertices vec(3*nv)*/
    Dbl * nvec;/* normal vec(3*nv)*/
    int nt; /* number of triangles */
    int nv; /* number of vertices */
} TSurf; /* triangulated surface */

int initTSurf(TSurf* s);
void cleanTSurf(TSurf* s); 
int makeTSurf(TSurf* s,int T,int N);
int cpTSurf(TSurf* dst,TSurf*src);
int writeTSurf(TSurf *s,char * fname);
int readTSurf(TSurf *s,char * fname);
int printTSurf(TSurf *s,char * fname);
int mkTetrahedron(TSurf *s);
int mkHexahedron(TSurf *s);
int mkOctahedron(TSurf *s);
int mkDodecahedron(TSurf *s);
int mkIcosahedron(TSurf *s);


Dbl trgNorm(Vec3Dbl w,TSurf *s,int t);
int getTrgCon(int *ntc,CList *tcvec,TSurf *s);
int getVrtCon(int *nvc,CList *vcvec,TSurf *s);
int getVrtBetw(int *nvc,CList *vcvec,TSurf *s,int k1,int k2);
int refineTSurf2(TSurf *s,TSurf *sold,int *nvc,CList *vcvec);
int checkNormals(Dbl * d,TSurf *s);
void mkUnitSphere(TSurf *s);
void printStat(TSurf *s);
#endif

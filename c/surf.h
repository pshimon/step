/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef SURF_H
#define SURF_H
#include "my_cdefs.h" 
#include "geom.h"
/* if distance between points is less than MIN_DIST 
 * they are considered coinciding
 */
#define MIN_DIST 1.0e-8
/* if norm of vector is less than MIN_NORM it is zero */
#define MIN_NORM 1.0e-12 
/* maximum number of connected vertices */
#define MAX_CONNECT 8
typedef int t_clst[MAX_CONNECT];
typedef struct {
    int * tvec;/*trg vector(3*nt);*/
    float * vvec;/* vertices vec(3*nv)*/
    float * nvec;/* normal vec(3*nv)*/
    int nt; /* number of triangles */
    int nv; /* number of vertices */
} TSurf;

/* unintialized variable may contain garbage */
int iniTSurf(TSurf* s);
/* to avoid memory leaks */
void cleanTSurf(TSurf* s); 
/* creates empty surface */
int makeTSurf(TSurf* s,int T,int N);
int cpTSurf(TSurf* dst,TSurf*src);
/* binary */
int writeTSurf(TSurf *s,char * fname);
int readTSurf(TSurf *s,char * fname);
int prinTSurf(TSurf *s,char * fname);
int mk_tetrahedron(TSurf *s);
int mk_hexahedron(TSurf *s);
int mk_octahedron(TSurf *s);
int mk_dodecahedron(TSurf *s);
int mk_icosahedron(TSurf *s);

/* rerturns are of trg */
float trg_norm(Vec3F w,TSurf *s,int t);
/* array lst must be of  MAX_CONNECT length at least*/
int get_trgs(int * lst,TSurf *s,int v);
int get_trg_pair(int * first,int * second,int * lst,int nc,TSurf *s,int v);
int remove_repeated_vertices(TSurf *s, int vstart);
int get_trg_con(int *ntc,t_clst *tcvec,TSurf *s);
int get_vrt_con(int *nvc,t_clst *vcvec,TSurf *s);
int get_vrt_betw(int *nvc,t_clst *vcvec,TSurf *s,int k1,int k2);
int refine2(TSurf *s,TSurf *sold,int *nvc,t_clst *vcvec);
int check_normals(float * d,TSurf *s);
void mk_unit_sphere(TSurf *s);
#endif

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
    t_f32 * vvec;/* vertices vec(3*nv)*/
    t_f32 * nvec;/* normal vec(3*nv)*/
    int nt; /* number of triangles */
    int nv; /* number of vertices */
} t_surf;

/* unintialized variable may contain garbage */
int init_surf(t_surf* s);
/* to avoid memory leaks */
void clean_surf(t_surf* s); 
/* creates empty surface */
int make_surf(t_surf* s,int T,int N);
int cp_surf(t_surf* dst,t_surf*src);
/* binary */
int write_surf(t_surf *s,char * fname);
int read_surf(t_surf *s,char * fname);
int print_surf(t_surf *s,char * fname);
int mk_tetrahedron(t_surf *s);
int mk_hexahedron(t_surf *s);
int mk_octahedron(t_surf *s);
int mk_dodecahedron(t_surf *s);
int mk_icosahedron(t_surf *s);

/* rerturns are of trg */
t_f32 trg_norm(t_v3 w,t_surf *s,int t);
/* array lst must be of  MAX_CONNECT length at least*/
int get_trgs(int * lst,t_surf *s,int v);
int get_trg_pair(int * first,int * second,int * lst,int nc,t_surf *s,int v);
int remove_repeated_vertices(t_surf *s, int vstart);
int get_trg_con(int *ntc,t_clst *tcvec,t_surf *s);
int get_vrt_con(int *nvc,t_clst *vcvec,t_surf *s);
int get_vrt_betw(int *nvc,t_clst *vcvec,t_surf *s,int k1,int k2);
int refine2(t_surf *s,t_surf *sold,int *nvc,t_clst *vcvec);
int check_normals(t_f32 * d,t_surf *s);
void mk_unit_sphere(t_surf *s);
#endif

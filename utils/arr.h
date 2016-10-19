/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef ARR_H
#define ARR_H
#include "my_cdefs.h"
/* fortran-style index: i1=1,n1 */
#define FEL1(a,i1) (a)->data[(i1)-1]
/* fortran-style index: i1=1,n1,i2=1,n2 col-major */
#define FEL2(a,i1,i2) (a)->data[(i1)-1+(a)->n1*((i2)-1)]
/* fortran-style index: i1=1,n1,i2=1,n2,i3=1,n3 col-major */
#define FEL3(a,i1,i2,i3) (a)->data[(i1)-1+(a)->n1*((i2)-1+(a)->n2*((i3)-1))]
/* fortran-style index: i1=1,n1,i2=1,n2,i3=1,n3,i4=1,n4 col-major */
#define FEL4(a,i1,i2,i3,i4) (a)->data[(i1)-1+(a)->n1*((i2)-1+(a)->n2*((i3)-1+(a)->n3*((i4)-1)))]
/* C-style index: i1=0,n1-1 */
#define CEL1(a,i1) (a)->data[(i1)]
/* C-style index: i1=0,n1-1,i2=0,n2-1 col-major */
#define CEL2(a,i1,i2) (a)->data[(i1)+(a)->n1*((i2))]
/* C-style index: i1=0,n1-1,i2=0,n2-1,i0=1,n3-1 col-major */
#define CEL3(a,i1,i2,i3) (a)->data[(i1)+(a)->n1*((i2)+(a)->n2*((i3)))]
/* C-style index: i1=0,n1-1,i2=0,n2-1,i3=0,n3-1,i4=0,n4-1 col-major */
#define CEL4(a,i1,i2,i3,i4) (a)->data[(i1)+(a)->n1*((i2)+(a)->n2*((i3)+(a)->n3*((i4))))]

typedef float Cfl[2];
typedef double Cdb[2];

/* float arrays*/
/* one dimensional  array */

typedef struct  {
    float *data;
    int n1;
} ArrFlt1;
void init_ArrFlt1(ArrFlt1 *a);
void clean_ArrFlt1(ArrFlt1 *a);
int make_ArrFlt1(ArrFlt1 *a,int n1);
int write_ArrFlt1(ArrFlt1 *a,char *file);
int read_ArrFlt1(ArrFlt1 *a,char *file);

/* two dimensional  array, matrix */

typedef struct  {
    float *data;
    int n1,n2;
} ArrFlt2;
void init_ArrFlt2(ArrFlt2 *a);
void clean_ArrFlt2(ArrFlt2 *a);
int make_ArrFlt2(ArrFlt2 *a,int n1,int n2);
int write_ArrFlt2(ArrFlt2 *a,char *file);
int read_ArrFlt2(ArrFlt2 *a,char *file);
int read_txArrFlt2(ArrFlt2 *a,char *file);
int write_txArrFlt2(ArrFlt2 *a,char *file);
/* 3d  array */

typedef struct  {
    float *data;
    int n1,n2,n3;
} ArrFlt3;
void init_ArrFlt3(ArrFlt3 *a);
void clean_ArrFlt3(ArrFlt3 *a);
int make_ArrFlt3(ArrFlt3 *a,int n1,int n2,int n3);
int write_ArrFlt3(ArrFlt3 *a,char *file);
int read_ArrFlt3(ArrFlt3 *a,char *file);

/* 4d  array */

typedef struct  {
    float *data;
    int n1,n2,n3,n4;
} ArrFlt4;
void init_ArrFlt4(ArrFlt4 *a);
void clean_ArrFlt4(ArrFlt4 *a);
int make_ArrFlt4(ArrFlt4 *a,int n1,int n2,int n3,int n4);
int write_ArrFlt4(ArrFlt4 *a,char *file);
int read_ArrFlt4(ArrFlt4 *a,char *file);

/* double arrays*/
/* one dimensional  array */

typedef struct  {
    double *data;
    int n1;
} ArrDbl1;
void init_ArrDbl1(ArrDbl1 *a);
void clean_ArrDbl1(ArrDbl1 *a);
int make_ArrDbl1(ArrDbl1 *a,int n1);
int write_ArrDbl1(ArrDbl1 *a,char *file);
int read_ArrDbl1(ArrDbl1 *a,char *file);

/* two dimensional  array, matrix */

typedef struct  {
    double *data;
    int n1,n2;
} ArrDbl2;
void init_ArrDbl2(ArrDbl2 *a);
void clean_ArrDbl2(ArrDbl2 *a);
int make_ArrDbl2(ArrDbl2 *a,int n1,int n2);
int write_ArrDbl2(ArrDbl2 *a,char *file);
int read_ArrDbl2(ArrDbl2 *a,char *file);

/* 3d  array */

typedef struct  {
    double *data;
    int n1,n2,n3;
} ArrDbl3;
void init_ArrDbl3(ArrDbl3 *a);
void clean_ArrDbl3(ArrDbl3 *a);
int make_ArrDbl3(ArrDbl3 *a,int n1,int n2,int n3);
int write_ArrDbl3(ArrDbl3 *a,char *file);
int read_ArrDbl3(ArrDbl3 *a,char *file);

/* 4d  array */

typedef struct  {
    double *data;
    int n1,n2,n3,n4;
} ArrDbl4;
void init_ArrDbl4(ArrDbl4 *a);
void clean_ArrDbl4(ArrDbl4 *a);
int make_ArrDbl4(ArrDbl4 *a,int n1,int n2,int n3,int n4);
int write_ArrDbl4(ArrDbl4 *a,char *file);
int read_ArrDbl4(ArrDbl4 *a,char *file);

/* int arrays*/
/* one dimensional  array */

typedef struct  {
    int *data;
    int n1;
} ArrInt1;
void init_ArrInt1(ArrInt1 *a);
void clean_ArrInt1(ArrInt1 *a);
int make_ArrInt1(ArrInt1 *a,int n1);
int write_ArrInt1(ArrInt1 *a,char *file);
int read_ArrInt1(ArrInt1 *a,char *file);

/* two dimensional  array, matrix */

typedef struct  {
    int *data;
    int n1,n2;
} ArrInt2;
void init_ArrInt2(ArrInt2 *a);
void clean_ArrInt2(ArrInt2 *a);
int make_ArrInt2(ArrInt2 *a,int n1,int n2);
int write_ArrInt2(ArrInt2 *a,char *file);
int read_ArrInt2(ArrInt2 *a,char *file);

/* 3d  array */

typedef struct  {
    int *data;
    int n1,n2,n3;
} ArrInt3;
void init_ArrInt3(ArrInt3 *a);
void clean_ArrInt3(ArrInt3 *a);
int make_ArrInt3(ArrInt3 *a,int n1,int n2,int n3);
int write_ArrInt3(ArrInt3 *a,char *file);
int read_ArrInt3(ArrInt3 *a,char *file);

/* 4d  array */

typedef struct  {
    int *data;
    int n1,n2,n3,n4;
} ArrInt4;
void init_ArrInt4(ArrInt4 *a);
void clean_ArrInt4(ArrInt4 *a);
int make_ArrInt4(ArrInt4 *a,int n1,int n2,int n3,int n4);
int write_ArrInt4(ArrInt4 *a,char *file);
int read_ArrInt4(ArrInt4 *a,char *file);

/* complex float arrays*/
/* one dimensional  array */
typedef struct  {
    Cfl *data;
    int n1;
} ArrCfl1;
void init_ArrCfl1(ArrCfl1 *a);
void clean_ArrCfl1(ArrCfl1 *a);
int make_ArrCfl1(ArrCfl1 *a,int n1);
int write_ArrCfl1(ArrCfl1 *a,char *file);
int read_ArrCfl1(ArrCfl1 *a,char *file);
int re_ArrCfl1(ArrFlt1 *rea,ArrCfl1 *a);
int im_ArrCfl1(ArrFlt1 *ima,ArrCfl1 *a);
int mod_ArrCfl1(ArrFlt1 *mda,ArrCfl1 *a);
int arg_ArrCfl1(ArrFlt1 *ara,ArrCfl1 *a);
int zip_ArrCfl1(ArrCfl1 *a,ArrFlt1 *rea,ArrFlt1 *ima);
/* two dimensional  array, matrix */

typedef struct  {
    Cfl *data;
    int n1,n2;
} ArrCfl2;
void init_ArrCfl2(ArrCfl2 *a);
void clean_ArrCfl2(ArrCfl2 *a);
int make_ArrCfl2(ArrCfl2 *a,int n1,int n2);
int write_ArrCfl2(ArrCfl2 *a,char *file);
int read_ArrCfl2(ArrCfl2 *a,char *file);

/* 3d  array */

typedef struct  {
    Cfl *data;
    int n1,n2,n3;
} ArrCfl3;
void init_ArrCfl3(ArrCfl3 *a);
void clean_ArrCfl3(ArrCfl3 *a);
int make_ArrCfl3(ArrCfl3 *a,int n1,int n2,int n3);
int write_ArrCfl3(ArrCfl3 *a,char *file);
int read_ArrCfl3(ArrCfl3 *a,char *file);


/* 4d  array */

typedef struct  {
    Cfl *data;
    int n1,n2,n3,n4;
} ArrCfl4;
void init_ArrCfl4(ArrCfl4 *a);
void clean_ArrCfl4(ArrCfl4 *a);
int make_ArrCfl4(ArrCfl4 *a,int n1,int n2,int n3,int n4);
int write_ArrCfl4(ArrCfl4 *a,char *file);
int read_ArrCfl4(ArrCfl4 *a,char *file);

/* complex double arrays*/
/* one dimensional  array */

typedef struct  {
    Cdb *data;
    int n1;
} ArrCdb1;
void init_ArrCdb1(ArrCdb1 *a);
void clean_ArrCdb1(ArrCdb1 *a);
int make_ArrCdb1(ArrCdb1 *a,int n1);
int write_ArrCdb1(ArrCdb1 *a,char *file);
int read_ArrCdb1(ArrCdb1 *a,char *file);

/* two dimensional  array, matrix */

typedef struct  {
    Cdb *data;
    int n1,n2;
} ArrCdb2;
void init_ArrCdb2(ArrCdb2 *a);
void clean_ArrCdb2(ArrCdb2 *a);
int make_ArrCdb2(ArrCdb2 *a,int n1,int n2);
int write_ArrCdb2(ArrCdb2 *a,char *file);
int read_ArrCdb2(ArrCdb2 *a,char *file);

/* 3d  array */

typedef struct  {
    Cdb *data;
    int n1,n2,n3;
} ArrCdb3;
void init_ArrCdb3(ArrCdb3 *a);
void clean_ArrCdb3(ArrCdb3 *a);
int make_ArrCdb3(ArrCdb3 *a,int n1,int n2,int n3);
int write_ArrCdb3(ArrCdb3 *a,char *file);
int read_ArrCdb3(ArrCdb3 *a,char *file);

/* 4d  array */

typedef struct  {
    Cdb *data;
    int n1,n2,n3,n4;
} ArrCdb4;
void init_ArrCdb4(ArrCdb4 *a);
void clean_ArrCdb4(ArrCdb4 *a);
int make_ArrCdb4(ArrCdb4 *a,int n1,int n2,int n3,int n4);
int write_ArrCdb4(ArrCdb4 *a,char *file);
int read_ArrCdb4(ArrCdb4 *a,char *file);


#endif

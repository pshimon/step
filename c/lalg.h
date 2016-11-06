/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef LALG_H
#define LALG_H
#include <mkl.h>
#include "databuf.h"
/*
inline static int getrfDbl(int n,Dbl *a,int *ipiv) {
    return LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,a,n,ipiv);
}

inline static int getrsDbl(int n,Dbl *a,int k,Dbl *b,int *ipiv) {
    return LAPACKE_dgetrs(LAPACK_COL_MAJOR,'N',n,k,a,n,ipiv,b,n);
}
*/
inline static int getrfDbl(int n,Dbl *a,int *ipiv) {
    int info;
    dgetrf_(&n,&n,a,&n,ipiv,&info);
    return info;
}
inline static int getrsDbl(int n,Dbl *a,int k,Dbl *b,int *ipiv) {
    int info;
    dgetrs_("N",&n,&k,a,&n,ipiv,b,&n,&info);
    return info;
}


#endif


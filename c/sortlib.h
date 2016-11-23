/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#ifndef SORTLIB_H
#define SORTLIB_H
#include "databuf.h"
#define ORDER(A, B) if (BEFORE(B, A)) SWAP(A, B)
#define SWITCHSORT 10
void sortSelDbl(Dbl a[], int l, int r);
void sortInsDbl(Dbl a[], int l, int r);
void sortShellDbl(Dbl a[], int l, int r);
void sortQuickDbl(Dbl a[], int l, int r);
void sortQuickNrDbl(Dbl a[], int l, int r);
#endif

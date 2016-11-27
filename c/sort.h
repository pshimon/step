/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef SORT_H
#define SORT_H
#include "my_cdefs.h"
/* Based upon but not copied from Numerical Recipes, 3 ed */
/* Sort an array  a[n]  in ascending order */
void heapsort(double * a,int n);
void heapsortf(float * a,int n);
/* Sort an array  a[n]  in ascending order, remove duplicates */
int hsort1(double * a,int *n,int * mult);
int hsort1f(float * a,int *n,int * mult);
/* Sort an array  a[n]  in ascending order with indexation 
 */
void heapsort_ind(double * a,int n,int * ind);
void heapsort_indf(float * a,int n,int * ind); 

#endif

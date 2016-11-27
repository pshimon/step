/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#include "sort.h"
void sift_down(Dbl * ra, int l,int r){
    int j,jold;
    Dbl a;
    a=ra[l];
    jold=l;
    j=2*l+1;
    while (j <= r) {
	if (j < r && ra[j] < ra[j+1]) j++; 
	if (a >= ra[j]) break; 
	ra[jold]=ra[j];
	jold=j;
	j=2*j+1;
    }
    ra[jold]=a; 
}
void heapsort(Dbl * ra,int n) {
    int i;
    Dbl tmp;
    for (i=n/2-1; i>=0; i--) sift_down(ra,i,n-1);
    for (i=n-1; i>0; i--) {
	tmp=ra[0];
	ra[0]=ra[i];
	ra[i]=tmp;
	sift_down(ra,0,i-1); 
    }
}
void sift_downf(float * ra, int l,int r){
    int j,jold;
    float a;
    a=ra[l];
    jold=l;
    j=2*l+1;
    while (j <= r) {
	if (j < r && ra[j] < ra[j+1]) j++; 
	if (a >= ra[j]) break; 
	ra[jold]=ra[j];
	jold=j;
	j=2*j+1;
    }
    ra[jold]=a; 
}
void heapsortf(float * ra,int n) {
    int i;
    float tmp;
    for (i=n/2-1; i>=0; i--) sift_downf(ra,i,n-1);
    for (i=n-1; i>0; i--) {
	tmp=ra[0];
	ra[0]=ra[i];
	ra[i]=tmp;
	sift_downf(ra,0,i-1); 
    }
}


int hsort1(Dbl * a,int *n,int * mult) {
    Dbl * tmp;
    int i,k;
    tmp=ALLOC_MEM(Dbl,*n);
    if(!tmp) return -1;
    for(i=0;i<*n;i++) {
	tmp[i]=a[i];
	mult[i]=0;
    }
    heapsort(tmp,*n);
    a[0]=tmp[0];
    mult[0]++;
    k=1;
    for(i=1;i<*n;i++) 
	if(tmp[i]>tmp[i-1]) 
	    a[k++]=tmp[i];
	else
	    mult[k]++;
    *n=k;
    FREE_MEM(tmp);
    return 0;	    

}
int hsort1f(float * a,int *n,int * mult) {
    float * tmp;
    int i,k;
    tmp=ALLOC_MEM(float,*n);
    if(!tmp) return -1;
    for(i=0;i<*n;i++) {
	tmp[i]=a[i];
	mult[i]=0;
    }
    heapsortf(tmp,*n);
    a[0]=tmp[0];
    mult[0]++;
    k=1;
    for(i=1;i<*n;i++) 
	if(tmp[i]>tmp[i-1]) 
	    a[k++]=tmp[i];
	else
	    mult[k]++;
    *n=k;
    FREE_MEM(tmp);
    return 0;	    

}

void sift_down_ind(Dbl * ra, int l,int r,int * ind){
    int j,jold,i;
    Dbl a;
    a=ra[l];
    i=ind[l];
    jold=l;
    j=2*l+1;
    while (j <= r) {
	if (j < r && ra[j] < ra[j+1]) j++; 
	if (a >= ra[j]) break; 
	ra[jold]=ra[j];
	ind[jold]=ind[j];
	jold=j;
	j=2*j+1;
    }
    ra[jold]=a;
    ind[jold]=i;
}

void heapsort_ind(Dbl * ra,int n,int * ind) {
    int i,ii;
    Dbl tmp;
    for (i=n/2-1; i>=0; i--) sift_down_ind(ra,i,n-1,ind);
    for (i=n-1; i>0; i--) {
	tmp=ra[0];
	ra[0]=ra[i];
	ra[i]=tmp;
	ii=ind[0];
	ind[0]=ind[i];
	ind[i]=ii;
	sift_down_ind(ra,0,i-1,ind); 
    }
}
void sift_down_indf(float * ra, int l,int r,int * ind){
    int j,jold,i;
    float a;
    a=ra[l];
    i=ind[l];
    jold=l;
    j=2*l+1;
    while (j <= r) {
	if (j < r && ra[j] < ra[j+1]) j++; 
	if (a >= ra[j]) break; 
	ra[jold]=ra[j];
	ind[jold]=ind[j];
	jold=j;
	j=2*j+1;
    }
    ra[jold]=a;
    ind[jold]=i;
}

void heapsort_indf(float * ra,int n,int * ind) {
    int i,ii;
    float tmp;
    for (i=n/2-1; i>=0; i--) sift_down_indf(ra,i,n-1,ind);
    for (i=n-1; i>0; i--) {
	tmp=ra[0];
	ra[0]=ra[i];
	ra[i]=tmp;
	ii=ind[0];
	ind[0]=ind[i];
	ind[i]=ii;
	sift_down_indf(ra,0,i-1,ind); 
    }
}




/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "sortlib.h"

/* adopted from NR and Sedgewick */
/* function names chosen to avoid any connotation with actual ordering
 * say before instead of BEFORE (it may be more actually)
 * first instead of minimal etc.
 */
#define BEFORE(A, B) ((A) < (B))
#define SWAP(A, B) { Dbl t = A; A = B; B = t; }

void sortSelDbl(Dbl a[], int l, int r) { 
    int i, j,min;
    for (i = l; i < r; i++) {  
	min = i;
	for (j = i+1; j <= r; j++) if (BEFORE(a[j], a[min])) min = j;
	SWAP(a[i], a[min]);
    } 
}

void sortInsDbl(Dbl a[], int l, int r) { 
    int i,j,min=l;
    Dbl v;
    for (i = l; i <= r; i++) if (BEFORE(a[i], a[min])) min = i;/* find location of element to be the first */
    if(min>l) SWAP(a[l], a[min]); 
    /* above step is not necessary, but useful: only one condition to check in while loop below */
    for (i = l+2; i <= r; i++) {  
	j = i; v = a[i]; 
	while (BEFORE(v, a[j-1])) { a[j] = a[j-1]; j--; }
	a[j] = v; 
    } 
}

void sortShellDbl(Dbl a[], int l, int r) { 
    int i, j, h; 
    Dbl v;
    for (h = 1; h <= (r-l)/9; h = 3*h+1);/* finding initial increment */
    for ( ; h > 0; h /= 3) {
	for (i = l+h; i <= r; i++) { 
	    j = i;  
	    v = a[i]; 
	    while (j >= l+h && BEFORE(v, a[j-h])) { a[j] = a[j-h]; j -= h; }
	    a[j] = v; 
	} 
    }
}

/* quick sort */

static int partition(Dbl a[], int l, int r) { 
    int i = l-1, j = r; 
    Dbl v = a[r];

    SWAP(a[(l+r)/2], a[r-1]);
    ORDER(a[l], a[r-1]); 
    ORDER(a[l], a[r]); 
    ORDER(a[r-1], a[r]);

    for (;;) { 
	while (BEFORE(a[++i], v)) ;
	while (BEFORE(v, a[--j])) if (j == l) break;
	if (i >= j) break;
	SWAP(a[i], a[j]);
    }
    SWAP(a[i], a[r]);
    return i;
}
/* partitions less than SWITCHSORT remain usnsorted */
static void quickSortDbl(Dbl a[], int l, int r) { 
    int i; 
    if (r-l <= SWITCHSORT) return;
    i = partition(a, l+1, r-1);
    quickSortDbl(a, l, i-1);
    quickSortDbl(a, i+1, r);
} 
void sortQuickDbl(Dbl a[], int l, int r) {
    quickSortDbl(a,l,r);
    sortInsDbl(a,l,r);	  
}
/* nor recursive implementation */

#define STACK_SIZE 32 /* enough for 32 bit   index */
int stackArr[STACK_SIZE];
int stackTop;
#define STACK_INIT stackTop=-1 
#define STACK_NOT_EMPTY stackTop>-1
#define PUSH(i) stackArr[++stackTop]=(i)
#define PUSH2(A, B)  PUSH(B); PUSH(A)
#define POP stackArr[stackTop--]

static void quickSortNrDbl(Dbl a[], int l, int r) { 
    int i;
    STACK_INIT; 
    PUSH2(l, r);
    while (STACK_NOT_EMPTY) {
	l = POP; r = POP; 
	if (r-l <= SWITCHSORT) continue;
	i = partition(a, l, r);
	if (i-l > r-i) 	{ /* push pointers to larger part on stack  first */
	    PUSH2(l, i-1); 
	    PUSH2(i+1, r);
	} else 	{ 
	    PUSH2(i+1, r); 
	    PUSH2(l, i-1); 
	}
    }
}

void sortQuickNrDbl(Dbl a[], int l, int r) {
    quickSortNrDbl(a,l,r);
    sortInsDbl(a,l,r);
}

/* heap sort */
static void siftDownDbl(Dbl  ra[], int l,int r){
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
void sortHeapDbl(Dbl  ra[],int l,int r) {
    int i;
    Dbl tmp;
    for (i=(r+1)/2-1; i>=l; i--) siftDownDbl(ra,i,r);
    for (i=r; i>l; i--) {
	tmp=ra[l];
	ra[l]=ra[i];
	ra[i]=tmp;
	siftDownDbl(ra,l,i-1); 
    }
}


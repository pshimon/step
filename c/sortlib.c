/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "databuf.h"

/* adopted from NR and Sedgewick */
/* function names chosen to avoid any connotation with actual ordering
 * say before instead of BEFORE (it may be more actually)
 * first instead of minimal etc.
 */
typedef Dbl Item;
#define KEY(A) (A)
#define BEFORE(A, B) (KEY(A) < KEY(B))
#define SWAP(A, B) { Item t = A; A = B; B = t; } 
#define ORDER(A, B) if (BEFORE(B, A)) SWAP(A, B)

void sortSel(Item a[], int l, int r) { 
    int i, j,min;
    for (i = l; i < r; i++) {  
	min = i;
	for (j = i+1; j <= r; j++) if (BEFORE(a[j], a[min])) min = j;/* find location of element to be the first */
	SWAP(a[i], a[min]);
    } 
}

void sortIns(Item a[], int l, int r) { 
    int i,j,min=l;
    Item v;
    for (i = l; i <= r; i++) if (BEFORE(a[i], a[min])) min = i;/* find location of element to be the first */
    if(min>l) SWAP(a[l], a[min]); 
    /* above step is not necessary, but useful: only one condition to check in while loop below */
    for (i = l+2; i <= r; i++) {  
	j = i; v = a[i]; 
	while (BEFORE(v, a[j-1])) { a[j] = a[j-1]; j--; }
	a[j] = v; 
    } 
}

void sortShell(Item a[], int l, int r) { 
    int i, j, h; 
    Item v;
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

static int partition(Item a[], int l, int r) { 
    int i = l-1, j = r; 
    Item v = a[r];

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

#define M 10
void sortQuick(Item a[], int l, int r) { 
    int i; 
    if (r-l <= M) return;
    i = partition(a, l+1, r-1);
    sortQuick(a, l, i-1);
    sortQuick(a, i+1, r);
    sortIns(a,l,r);
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

void sortQuickNr(Item a[], int l, int r) { 
    int i;
    STACK_INIT; 
    PUSH2(l, r);
    while (STACK_NOT_EMPTY) {
	l = POP; r = POP; 
	if (r-l <= M) continue;
	i = partition(a, l, r);
	if (i-l > r-i) 	{ /* push pointers to larger part on stack  first */
	    PUSH2(l, i-1); 
	    PUSH2(i+1, r);
	} else 	{ 
	    PUSH2(i+1, r); 
	    PUSH2(l, i-1); 
	}
    }
    sortIns(a,l,r);
}
/* usage :
 
void sort(Item a[], int l, int r)
  { 
    quicksort(a, l, r);
    insertion(a, l, r);
  }
*/

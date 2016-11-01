#include <stdio.h>
#include <conio.h>

#define MAXELT          100
#define INFINITY        32760         // numbers in list should not exceed
                                      // this. change the value to suit your
                                      // needs
#define SMALLSIZE       10            // not less than 3
#define STACKSIZE       100           // should be ceiling(lg(MAXSIZE)+1)

int list[MAXELT+1];                   // one extra, to hold INFINITY

struct {                              // stack element.
        int a,b;
} stack[STACKSIZE];

int top=-1;                           // initialise stack

int main()                           // overhead!
{
    int i=-1,j,n;
    char t[10];
    void quicksort(int);
        
    do {
        if (i!=-1)
            list[i++]=n;
        else
            i++;
        printf("Enter the numbers <End by #>: ");
        fflush(stdin);
        scanf("%[^\n]",t);
        if (sscanf(t,"%d",&n)<1)
        break;
    } while (1);
        
    quicksort(i-1);
        
    printf("\nThe list obtained is ");
    for (j=0;j<i;j++)
        printf("\n %d",list[j]);

    printf("\n\nProgram over.");
    getch();
    return 0;       // successful termination.
}

void interchange(int *x,int *y)        // swap
{
    int temp;
        
    temp=*x;
    *x=*y;
    *y=temp;
}

void split(int first,int last,int *splitpoint)
{
    int x,i,j,s,g;
    
    // here, atleast three elements are needed
    if (list[first]<list[(first+last)/2]) {  // find median
        s=first;
        g=(first+last)/2;
    }
    else {
        g=first;
        s=(first+last)/2;
    }
    if (list[last]<=list[s]) 
        x=s;
    else if (list[last]<=list[g])
        x=last;
    else
        x=g;
    interchange(&list[x],&list[first]);      // swap the split-point element
                                             // with the first
    x=list[first];
    i=first+1;                               // initialise
    j=last+1;
    while (i<j) {
        do {                                 // find j 
            j--;
        } while (list[j]>x);
        do {
            i++;                             // find i
        } while (list[i]<x);
        interchange(&list[i],&list[j]);      // swap
    }
    interchange(&list[i],&list[j]);          // undo the extra swap
    interchange(&list[first],&list[j]);      // bring the split-point 
                                             // element to the first
    *splitpoint=j;
}

void push(int a,int b)                        // push
{
    top++;
    stack[top].a=a;
    stack[top].b=b;
}

void pop(int *a,int *b)                       // pop
{
    *a=stack[top].a;
    *b=stack[top].b;
    top--;
}

void insertion_sort(int first,int last)
{
    int i,j,c;
        
    for (i=first;i<=last;i++) {
        j=list[i];
        c=i;
        while ((list[c-1]>j)&&(c>first)) {
            list[c]=list[c-1];
            c--;
        }
        list[c]=j;
    }
}

void quicksort(int n)
{
    int first,last,splitpoint;
        
    push(0,n);
    while (top!=-1) {
        pop(&first,&last);
        for (;;) {
            if (last-first>SMALLSIZE) {
                // find the larger sub-list
                split(first,last,&splitpoint);
                // push the smaller list
                if (last-splitpoint<splitpoint-first) {
                    push(first,splitpoint-1);
                    first=splitpoint+1;
                }
                else {
                    push(splitpoint+1,last);
                    last=splitpoint-1;
                }
            }
            else {  // sort the smaller sub-lists
                    // through insertion sort
                insertion_sort(first,last);
                break;
            }
        }
    }                        // iterate for larger list
}



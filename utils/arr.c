/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#include "arr.h"
/* float arrays */

void init_ArrFlt1(ArrFlt1 *a) {
    a->data=0;
    a->n1=0;
}

void clean_ArrFlt1(ArrFlt1 *a) {
    FREE_MEM(a->data);
    init_ArrFlt1(a);
}

int make_ArrFlt1(ArrFlt1 *a,int n1) {
    int oldsize=a->n1;
    if(n1<1) return 1;
    if(n1!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(float,n1);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    return 0;
}
	
int write_ArrFlt1(ArrFlt1 *a,char *file) {
    FILE* fp;
    size_t n=a->n1;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(float),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrFlt1(ArrFlt1 *a,char *file) {
    int ret,n1;
    FILE* fp;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}    
    if(make_ArrFlt1(a,n1)){ret=-3;goto abend;};
    if(fread(a->data,sizeof(float),n1,fp)!=n1) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 


void init_ArrFlt2(ArrFlt2 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
}

void clean_ArrFlt2(ArrFlt2 *a) {
    FREE_MEM(a->data);
    init_ArrFlt2(a);
}

int make_ArrFlt2(ArrFlt2 *a,int n1,int n2) {
    int oldsize=a->n1*a->n2;
    int newsize=n1*n2;
    if((n1<1)||(n2<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(float,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    return 0;	

}
	
int write_ArrFlt2(ArrFlt2 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(float),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrFlt2(ArrFlt2 *a,char *file) {
    int ret,n1,n2;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrFlt2(a,n1,n2)){ret=-3;goto abend;};
    n=n1*n2;
    if(fread(a->data,sizeof(float),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 
int read_txArrFlt2(ArrFlt2 *a,char *file) {
    FILE* fp;
    int j,k,cnt;
    float v;
    fp=fopen(file,"r");
    if(0==fp) return -1;
    if(fscanf(fp,"%d %d",&j,&k)!=2) return -2;
    if(make_ArrFlt2(a,j,k)) return -4;
    cnt=0;
    for(j=0;j<a->n1;j++) {
	for(k=0;k<a->n2;k++) {
	    if(fscanf(fp,"%f",&v)!=1) break;
	    CEL2(a,j,k)=v;
	    cnt++;
	}
    }
    if(cnt!=a->n1*a->n2)return -3; 
    fclose(fp);
    return 0;
}
int write_txArrFlt2(ArrFlt2 *a,char *file) {
    FILE* fp;
    int j,k;
    //float v;
    fp=fopen(file,"w");
    if(0==fp) return -1;
    fprintf(fp,"%d %d\n",a->n1,a->n2);
    for(j=0;j<a->n1;j++) {
	for(k=0;k<a->n2;k++) fprintf(fp,"%e\t",CEL2(a,j,k));
	fprintf(fp,"\n");
    }
    fclose(fp);
    return 0;
}

void init_ArrFlt3(ArrFlt3 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
    a->n3=0;
}

void clean_ArrFlt3(ArrFlt3 *a) {
    FREE_MEM(a->data);
    init_ArrFlt3(a);
}

int make_ArrFlt3(ArrFlt3 *a,int n1,int n2,int n3) {
    int oldsize=a->n1*a->n2*a->n3;
    int newsize=n1*n2*n3;
    if((n1<1)||(n2<1)||(n3<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(float,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    a->n3=n3;
    return 0;	

}
	
int write_ArrFlt3(ArrFlt3 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2*a->n3;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(float),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrFlt3(ArrFlt3 *a,char *file) {
    int ret,n1,n2,n3;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrFlt3(a,n1,n2,n3)){ret=-3;goto abend;};
    n=n1*n2*n3;
    if(fread(a->data,sizeof(float),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 

void init_ArrFlt4(ArrFlt4 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
    a->n3=0;
    a->n4=0;
}

void clean_ArrFlt4(ArrFlt4 *a) {
    FREE_MEM(a->data);
    init_ArrFlt4(a);
}

int make_ArrFlt4(ArrFlt4 *a,int n1,int n2,int n3,int n4) {
    int oldsize=a->n1*a->n2*a->n3*a->n4;
    int newsize=n1*n2*n3*n4;
    if((n1<1)||(n2<1)||(n3<1)||(n4<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(float,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    a->n3=n3;
    a->n4=n4;
    return 0;	

}
	
int write_ArrFlt4(ArrFlt4 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2*a->n3*a->n4;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n4,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(float),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrFlt4(ArrFlt4 *a,char *file) {
    int ret,n1,n2,n3,n4;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n4,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrFlt4(a,n1,n2,n3,n4)){ret=-3;goto abend;};
    n=n1*n2*n3*n4;
    if(fread(a->data,sizeof(float),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 

/* double arrays*/

void init_ArrDbl1(ArrDbl1 *a) {
    a->data=0;
    a->n1=0;
}

void clean_ArrDbl1(ArrDbl1 *a) {
    FREE_MEM(a->data);
    init_ArrDbl1(a);
}

int make_ArrDbl1(ArrDbl1 *a,int n1) {
    int oldsize=a->n1;
    if(n1<1) return 1;
    if(n1!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(double,n1);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    return 0;
}
	
int write_ArrDbl1(ArrDbl1 *a,char *file) {
    FILE* fp;
    size_t n=a->n1;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(double),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrDbl1(ArrDbl1 *a,char *file) {
    int ret,n1;
    FILE* fp;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}    
    if(make_ArrDbl1(a,n1)){ret=-3;goto abend;};
    if(fread(a->data,sizeof(double),n1,fp)!=n1) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 


void init_ArrDbl2(ArrDbl2 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
}

void clean_ArrDbl2(ArrDbl2 *a) {
    FREE_MEM(a->data);
    init_ArrDbl2(a);
}

int make_ArrDbl2(ArrDbl2 *a,int n1,int n2) {
    int oldsize=a->n1*a->n2;
    int newsize=n1*n2;
    if((n1<1)||(n2<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(double,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    return 0;	

}
	
int write_ArrDbl2(ArrDbl2 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(double),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrDbl2(ArrDbl2 *a,char *file) {
    int ret,n1,n2;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrDbl2(a,n1,n2)){ret=-3;goto abend;};
    n=n1*n2;
    if(fread(a->data,sizeof(double),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 


void init_ArrDbl3(ArrDbl3 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
    a->n3=0;
}

void clean_ArrDbl3(ArrDbl3 *a) {
    FREE_MEM(a->data);
    init_ArrDbl3(a);
}

int make_ArrDbl3(ArrDbl3 *a,int n1,int n2,int n3) {
    int oldsize=a->n1*a->n2*a->n3;
    int newsize=n1*n2*n3;
    if((n1<1)||(n2<1)||(n3<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(double,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    a->n3=n3;
    return 0;	

}
	
int write_ArrDbl3(ArrDbl3 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2*a->n3;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(double),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrDbl3(ArrDbl3 *a,char *file) {
    int ret,n1,n2,n3;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrDbl3(a,n1,n2,n3)){ret=-3;goto abend;};
    n=n1*n2*n3;
    if(fread(a->data,sizeof(double),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 

void init_ArrDbl4(ArrDbl4 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
    a->n3=0;
    a->n4=0;
}

void clean_ArrDbl4(ArrDbl4 *a) {
    FREE_MEM(a->data);
    init_ArrDbl4(a);
}

int make_ArrDbl4(ArrDbl4 *a,int n1,int n2,int n3,int n4) {
    int oldsize=a->n1*a->n2*a->n3*a->n4;
    int newsize=n1*n2*n3*n4;
    if((n1<1)||(n2<1)||(n3<1)||(n4<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(double,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    a->n3=n3;
    a->n4=n4;
    return 0;	

}
	
int write_ArrDbl4(ArrDbl4 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2*a->n3*a->n4;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n4,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(double),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrDbl4(ArrDbl4 *a,char *file) {
    int ret,n1,n2,n3,n4;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n4,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrDbl4(a,n1,n2,n3,n4)){ret=-3;goto abend;};
    n=n1*n2*n3*n4;
    if(fread(a->data,sizeof(double),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 

/* complex float arrays */

void init_ArrCfl1(ArrCfl1 *a) {
    a->data=0;
    a->n1=0;
}

void clean_ArrCfl1(ArrCfl1 *a) {
    FREE_MEM(a->data);
    init_ArrCfl1(a);
}

int make_ArrCfl1(ArrCfl1 *a,int n1) {
    int oldsize=a->n1;
    if(n1<1) return 1;
    if(n1!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(Cfl,n1);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    return 0;
}
	
int write_ArrCfl1(ArrCfl1 *a,char *file) {
    FILE* fp;
    size_t n=a->n1;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(Cfl),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrCfl1(ArrCfl1 *a,char *file) {
    int ret,n1;
    FILE* fp;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}    
    if(make_ArrCfl1(a,n1)){ret=-3;goto abend;};
    if(fread(a->data,sizeof(Cfl),n1,fp)!=n1) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 
int re_ArrCfl1(ArrFlt1 *rea,ArrCfl1 *a) {
    int i,n=a->n1;
    if(make_ArrFlt1(rea,n)) return -1;
    for(i=0;i<n;i++) rea->data[i]=a->data[i][0];
    return 0;
}
int im_ArrCfl1(ArrFlt1 *ima,ArrCfl1 *a) {
    int i,n=a->n1;
    if(make_ArrFlt1(ima,n)) return -1;
    for(i=0;i<n;i++) ima->data[i]=a->data[i][1];
    return 0;
}
int mod_ArrCfl1(ArrFlt1 *mda,ArrCfl1 *a) {
    int i,n=a->n1;
    if(make_ArrFlt1(mda,n)) return -1;
    for(i=0;i<n;i++) mda->data[i]=HPT(a->data[i][0],a->data[i][1]);
    return 0;
}
int arg_ArrCfl1(ArrFlt1 *ara,ArrCfl1 *a) {
    int i,n=a->n1;
    if(make_ArrFlt1(ara,n)) return -1;
    for(i=0;i<n;i++) ara->data[i]=ARG(a->data[i][0],a->data[i][1]);
    return 0;
}

int zip_ArrCfl1(ArrCfl1 *a,ArrFlt1 *rea,ArrFlt1 *ima) {
    int i,n=rea->n1;
    if(ima->n1!=n) return -2;
    if(make_ArrCfl1(a,n)) return -1;
    for(i=0;i<n;i++) a->data[i][0]=a->data[i][1]=0.0;
    if(rea) for(i=0;i<n;i++) a->data[i][0]=rea->data[i];
    if(ima) for(i=0;i<n;i++) a->data[i][1]=ima->data[i];
    return 0;
}
void init_ArrCfl2(ArrCfl2 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
}

void clean_ArrCfl2(ArrCfl2 *a) {
    FREE_MEM(a->data);
    init_ArrCfl2(a);
}

int make_ArrCfl2(ArrCfl2 *a,int n1,int n2) {
    int oldsize=a->n1*a->n2;
    int newsize=n1*n2;
    if((n1<1)||(n2<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(Cfl,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    return 0;	

}
	
int write_ArrCfl2(ArrCfl2 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(Cfl),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrCfl2(ArrCfl2 *a,char *file) {
    int ret,n1,n2;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrCfl2(a,n1,n2)){ret=-3;goto abend;};
    n=n1*n2;
    if(fread(a->data,sizeof(Cfl),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 


void init_ArrCfl3(ArrCfl3 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
    a->n3=0;
}

void clean_ArrCfl3(ArrCfl3 *a) {
    FREE_MEM(a->data);
    init_ArrCfl3(a);
}

int make_ArrCfl3(ArrCfl3 *a,int n1,int n2,int n3) {
    int oldsize=a->n1*a->n2*a->n3;
    int newsize=n1*n2*n3;
    if((n1<1)||(n2<1)||(n3<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(Cfl,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    a->n3=n3;
    return 0;	

}
	
int write_ArrCfl3(ArrCfl3 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2*a->n3;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(Cfl),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrCfl3(ArrCfl3 *a,char *file) {
    int ret,n1,n2,n3;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrCfl3(a,n1,n2,n3)){ret=-3;goto abend;};
    n=n1*n2*n3;
    if(fread(a->data,sizeof(Cfl),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 

void init_ArrCfl4(ArrCfl4 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
    a->n3=0;
    a->n4=0;
}

void clean_ArrCfl4(ArrCfl4 *a) {
    FREE_MEM(a->data);
    init_ArrCfl4(a);
}

int make_ArrCfl4(ArrCfl4 *a,int n1,int n2,int n3,int n4) {
    int oldsize=a->n1*a->n2*a->n3*a->n4;
    int newsize=n1*n2*n3*n4;
    if((n1<1)||(n2<1)||(n3<1)||(n4<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(Cfl,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    a->n3=n3;
    a->n4=n4;
    return 0;	

}
	
int write_ArrCfl4(ArrCfl4 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2*a->n3*a->n4;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n4,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(Cfl),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrCfl4(ArrCfl4 *a,char *file) {
    int ret,n1,n2,n3,n4;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n4,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrCfl4(a,n1,n2,n3,n4)){ret=-3;goto abend;};
    n=n1*n2*n3*n4;
    if(fread(a->data,sizeof(Cfl),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 

/* comlex double arrays */

void init_ArrCdb1(ArrCdb1 *a) {
    a->data=0;
    a->n1=0;
}

void clean_ArrCdb1(ArrCdb1 *a) {
    FREE_MEM(a->data);
    init_ArrCdb1(a);
}

int make_ArrCdb1(ArrCdb1 *a,int n1) {
    int oldsize=a->n1;
    if(n1<1) return 1;
    if(n1!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(Cdb,n1);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    return 0;
}
	
int write_ArrCdb1(ArrCdb1 *a,char *file) {
    FILE* fp;
    size_t n=a->n1;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(Cdb),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrCdb1(ArrCdb1 *a,char *file) {
    int ret,n1;
    FILE* fp;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}    
    if(make_ArrCdb1(a,n1)){ret=-3;goto abend;};
    if(fread(a->data,sizeof(Cdb),n1,fp)!=n1) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 


void init_ArrCdb2(ArrCdb2 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
}

void clean_ArrCdb2(ArrCdb2 *a) {
    FREE_MEM(a->data);
    init_ArrCdb2(a);
}

int make_ArrCdb2(ArrCdb2 *a,int n1,int n2) {
    int oldsize=a->n1*a->n2;
    int newsize=n1*n2;
    if((n1<1)||(n2<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(Cdb,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    return 0;	

}
	
int write_ArrCdb2(ArrCdb2 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(Cdb),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrCdb2(ArrCdb2 *a,char *file) {
    int ret,n1,n2;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrCdb2(a,n1,n2)){ret=-3;goto abend;};
    n=n1*n2;
    if(fread(a->data,sizeof(Cdb),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 


void init_ArrCdb3(ArrCdb3 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
    a->n3=0;
}

void clean_ArrCdb3(ArrCdb3 *a) {
    FREE_MEM(a->data);
    init_ArrCdb3(a);
}

int make_ArrCdb3(ArrCdb3 *a,int n1,int n2,int n3) {
    int oldsize=a->n1*a->n2*a->n3;
    int newsize=n1*n2*n3;
    if((n1<1)||(n2<1)||(n3<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(Cdb,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    a->n3=n3;
    return 0;	

}
	
int write_ArrCdb3(ArrCdb3 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2*a->n3;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(Cdb),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrCdb3(ArrCdb3 *a,char *file) {
    int ret,n1,n2,n3;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrCdb3(a,n1,n2,n3)){ret=-3;goto abend;};
    n=n1*n2*n3;
    if(fread(a->data,sizeof(Cdb),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 

void init_ArrCdb4(ArrCdb4 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
    a->n3=0;
    a->n4=0;
}

void clean_ArrCdb4(ArrCdb4 *a) {
    FREE_MEM(a->data);
    init_ArrCdb4(a);
}

int make_ArrCdb4(ArrCdb4 *a,int n1,int n2,int n3,int n4) {
    int oldsize=a->n1*a->n2*a->n3*a->n4;
    int newsize=n1*n2*n3*n4;
    if((n1<1)||(n2<1)||(n3<1)||(n4<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(Cdb,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    a->n3=n3;
    a->n4=n4;
    return 0;	

}
	
int write_ArrCdb4(ArrCdb4 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2*a->n3*a->n4;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n4,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(Cdb),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrCdb4(ArrCdb4 *a,char *file) {
    int ret,n1,n2,n3,n4;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n4,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrCdb4(a,n1,n2,n3,n4)){ret=-3;goto abend;};
    n=n1*n2*n3*n4;
    if(fread(a->data,sizeof(Cdb),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 

/* int arrays */

void init_ArrInt1(ArrInt1 *a) {
    a->data=0;
    a->n1=0;
}

void clean_ArrInt1(ArrInt1 *a) {
    FREE_MEM(a->data);
    init_ArrInt1(a);
}

int make_ArrInt1(ArrInt1 *a,int n1) {
    int oldsize=a->n1;
    if(n1<1) return 1;
    if(n1!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(int,n1);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    return 0;
}
	
int write_ArrInt1(ArrInt1 *a,char *file) {
    FILE* fp;
    size_t n=a->n1;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(int),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrInt1(ArrInt1 *a,char *file) {
    int ret,n1;
    FILE* fp;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}    
    if(make_ArrInt1(a,n1)){ret=-3;goto abend;};
    if(fread(a->data,sizeof(int),n1,fp)!=n1) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 


void init_ArrInt2(ArrInt2 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
}

void clean_ArrInt2(ArrInt2 *a) {
    FREE_MEM(a->data);
    init_ArrInt2(a);
}

int make_ArrInt2(ArrInt2 *a,int n1,int n2) {
    int oldsize=a->n1*a->n2;
    int newsize=n1*n2;
    if((n1<1)||(n2<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(int,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    return 0;	

}
	
int write_ArrInt2(ArrInt2 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(int),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrInt2(ArrInt2 *a,char *file) {
    int ret,n1,n2;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrInt2(a,n1,n2)){ret=-3;goto abend;};
    n=n1*n2;
    if(fread(a->data,sizeof(int),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 


void init_ArrInt3(ArrInt3 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
    a->n3=0;
}

void clean_ArrInt3(ArrInt3 *a) {
    FREE_MEM(a->data);
    init_ArrInt3(a);
}

int make_ArrInt3(ArrInt3 *a,int n1,int n2,int n3) {
    int oldsize=a->n1*a->n2*a->n3;
    int newsize=n1*n2*n3;
    if((n1<1)||(n2<1)||(n3<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(int,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    a->n3=n3;
    return 0;	

}
	
int write_ArrInt3(ArrInt3 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2*a->n3;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(int),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrInt3(ArrInt3 *a,char *file) {
    int ret,n1,n2,n3;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrInt3(a,n1,n2,n3)){ret=-3;goto abend;};
    n=n1*n2*n3;
    if(fread(a->data,sizeof(int),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 

void init_ArrInt4(ArrInt4 *a) {
    a->data=0;
    a->n1=0;
    a->n2=0;
    a->n3=0;
    a->n4=0;
}

void clean_ArrInt4(ArrInt4 *a) {
    FREE_MEM(a->data);
    init_ArrInt4(a);
}

int make_ArrInt4(ArrInt4 *a,int n1,int n2,int n3,int n4) {
    int oldsize=a->n1*a->n2*a->n3*a->n4;
    int newsize=n1*n2*n3*n4;
    if((n1<1)||(n2<1)||(n3<1)||(n4<1)) return 1;
    if(newsize!=oldsize) {
	FREE_MEM(a->data);
	a->data=ALLOC_MEM(int,newsize);
	if(a->data==0) return -1;
    }
    a->n1=n1;
    a->n2=n2;
    a->n3=n3;
    a->n4=n4;
    return 0;	

}
	
int write_ArrInt4(ArrInt4 *a,char *file) {
    FILE* fp;
    size_t n=a->n1*a->n2*a->n3*a->n4;
    int ret;
    fp=fopen(file,"wb");
    if(!fp) return -1;
    if(fwrite(&a->n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(&a->n4,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fwrite(a->data,sizeof(int),n,fp)!=n) {ret=-3;goto abend;}
    ret=0;
abend:    
    fclose(fp);
    return ret;
}

int read_ArrInt4(ArrInt4 *a,char *file) {
    int ret,n1,n2,n3,n4;
    FILE* fp;
    size_t n;
    fp=fopen(file,"rb");
    if(!fp) return -1;
    if(fread(&n1,sizeof(int),1,fp)!=1) {ret=-2;goto abend;} 
    if(fread(&n2,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n3,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(fread(&n4,sizeof(int),1,fp)!=1) {ret=-2;goto abend;}
    if(make_ArrInt4(a,n1,n2,n3,n4)){ret=-3;goto abend;};
    n=n1*n2*n3*n4;
    if(fread(a->data,sizeof(int),n,fp)!=n) {ret=-3; goto abend;}
    ret=0;
abend:	
    fclose(fp);
    return  ret;
} 


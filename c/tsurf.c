/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "tsurf.h"

int initTSurf(TSurf* s) {
/* unintialized variable may contain garbage */    
    if(s==0) return -1;
    s->vvec=0;s->nvec=0;s->nv=0;
    s->tvec=0;s->nt=0;
    return 0;
}

void cleanTSurf(TSurf* s) {
/* to avoid memory leaks */    
    if(s==0) return;
    FREE_MEM(s->vvec);
    FREE_MEM(s->tvec);
    FREE_MEM(s->nvec);
    s->nt=0;s->nv=0;
}

int makeTSurf(TSurf* s,int T,int N) {
/* creates empty surface */    
    int ret;
    if(0==s) return -1;
    cleanTSurf(s);
    if(N<3) {ret=-2;goto abend;}
    if(T<1) {ret=-3;goto abend;}
    s->nt=T;
    s->tvec=ALLOC_MEM(int,3*T);
    if(0==s->tvec) {ret=-4;goto abend;}
    s->nv=N;
    s->vvec=ALLOC_MEM(Dbl,3*N);
    if(0==s->vvec) {ret=-5;goto abend;}
    s->nvec=ALLOC_MEM(Dbl,3*N);
    if(0==s->nvec) {ret=-6;goto abend;}
    return 0;
abend:
     cleanTSurf(s);
    return ret;   
}

int cpTSurf(TSurf* dst,TSurf* src) {
    if(makeTSurf(dst,src->nt,src->nv)) return -1;
    memcpy(dst->tvec,src->tvec,3*src->nt*sizeof(int));
    memcpy(dst->vvec,src->vvec,3*src->nv*sizeof(Dbl));
    memcpy(dst->nvec,src->nvec,3*src->nv*sizeof(Dbl));
    return 0;
}
int writeTSurf(TSurf *s,char * fname) {
    FILE *fp=fopen(fname,"wb");
    size_t l;
    if(fp==0) return -1;
    fwrite(&s->nt,sizeof(int),1,fp);
    fwrite(&s->nv,sizeof(int),1,fp);
    l=3*s->nt;
    if(fwrite(s->tvec,sizeof(int),l,fp)!= l) return -2;
    l=3*s->nv;
    if(fwrite(s->vvec,sizeof(Dbl),l,fp)!= l) return -3;
    if(fwrite(s->nvec,sizeof(Dbl),l,fp)!= l) return -4;
    fclose(fp);
    return 0;
}
int readTSurf(TSurf *s,char * fname) {
    FILE *fp=fopen(fname,"rb");
    int n,t,ret;
    size_t l;
    if(fp==0) return -1;
    if(fread(&t,sizeof(int),1,fp)!=1) {ret=-5;goto abend;}
    if(fread(&n,sizeof(int),1,fp)!=1) {ret=-5;goto abend;}
    if(makeTSurf(s,t,n)) {ret=-2;goto abend;}
    l=3*t;
    if(fread(s->tvec,sizeof(int), l,fp)!= l) {ret=-3;goto abend;}
    l=3*n;
    if(fread(s->vvec,sizeof(Dbl), l,fp)!= l) {ret=-4;goto abend;}
    if(fread(s->nvec,sizeof(Dbl), l,fp)!= l) {ret=-4;goto abend;}
    ret=0;
abend:
    fclose(fp);    
    return ret;
}

int printTSurf(TSurf *s,char * fname) {
    FILE *fp=fopen(fname,"w");
    int j;
    if(fp==0) return -1;
    fprintf(fp,"\t%d\t%d\n",s->nt,s->nv);
    for(j=0;j<s->nt;j++) fprintf(fp,"\t%d\t%d\t%d\n",s->tvec[3*j+0],
	    s->tvec[3*j+1],s->tvec[3*j+2]);
    for(j=0;j<s->nv;j++) fprintf(fp,"%f %f %f\n",s->vvec[3*j+0],
	    s->vvec[3*j+1],s->vvec[3*j+2]);
    for(j=0;j<s->nv;j++) fprintf(fp,"%f %f %f\n",s->nvec[3*j+0],
	    s->nvec[3*j+1],s->nvec[3*j+2]);
    fclose(fp);
    return 0;
}

int  mkTetrahedron(TSurf *s) {
    int nt=4,nv=4,j;
    int vt[]={0,1,2,0,3,1,0,2,3,1,3,2};
    Dbl vv[]={ 0.577350259,0.577350259,0.577350259,
	0.577350259,     -0.577350259,     -0.577350259,
	-0.577350259 ,     0.577350259,     -0.577350259,
	-0.577350259,     -0.577350259,      0.577350259};
    Dbl vn[]=  {0.577350318,      0.577350318,      0.577350318,    
	0.577350318,     -0.577350318,     -0.577350318,    
	-0.577350318,      0.577350318,     -0.577350318,    
	-0.577350318,     -0.577350318,      0.577350318};    
    if(makeTSurf(s,nt,nv)) return -1;
    for(j=0;j<3*nt;j++)    s->tvec[j]=vt[j];
    for(j=0;j<3*nv;j++)  {
	s->vvec[j]=vv[j];
	s->nvec[j]=vn[j];
    }
    return 0;
}
int mkHexahedron(TSurf *s) {
    int nt=24,nv=14,j;
    int vt[]={8,           0,           2,
           8,           2,           3,
           8,           3,           1,
           8,           1,           0,
           9,           4,           5,
           9,           5,           7,
           9,           7,           6,
           9,           6,           4,
          10,           0,           1,
          10,           1,           5,
          10,           5,           4,
          10,           4,           0,
          11,           2,           6,
          11,           6,           7,
          11,           7,           3,
          11,           3,           2,
          12,           0,           4,
          12,           4,           6,
          12,           6,           2,
          12,           2,           0,
          13,           3,           7,
          13,           7,           5,
          13,           5,           1,
          13,           1,           3};
    Dbl vv[]={   0.577350259     ,  0.577350259     ,  0.577350259     ,
  0.577350259     ,  0.577350259     , -0.577350259     ,
  0.577350259     , -0.577350259     ,  0.577350259     ,
  0.577350259     , -0.577350259     , -0.577350259     ,
 -0.577350259     ,  0.577350259     ,  0.577350259     ,
 -0.577350259     ,  0.577350259     , -0.577350259     ,
 -0.577350259     , -0.577350259     ,  0.577350259     ,
 -0.577350259     , -0.577350259     , -0.577350259     ,
  0.577350259     ,   0.00000000     ,   0.00000000     ,
 -0.577350259     ,   0.00000000     ,   0.00000000     ,
   0.00000000     ,  0.577350259     ,   0.00000000     ,
   0.00000000     , -0.577350259     ,   0.00000000     ,
   0.00000000     ,   0.00000000     ,  0.577350259     ,
   0.00000000     ,   0.00000000     , -0.577350259     
};
    Dbl vn[]=  { 0.577350318     ,  0.577350318     ,  0.577350318     ,
  0.577350318     ,  0.577350318     , -0.577350318     ,
  0.577350318     , -0.577350318     ,  0.577350318     ,
  0.577350318     , -0.577350318     , -0.577350318     ,
 -0.577350318     ,  0.577350318     ,  0.577350318     ,
 -0.577350318     ,  0.577350318     , -0.577350318     ,
 -0.577350318     , -0.577350318     ,  0.577350318     ,
 -0.577350318     , -0.577350318     , -0.577350318     ,
   1.00000000     ,   0.00000000     ,   0.00000000     ,
  -1.00000000     ,   0.00000000     ,   0.00000000     ,
   0.00000000     ,   1.00000000     ,   0.00000000     ,
   0.00000000     ,  -1.00000000     ,   0.00000000     ,
   0.00000000     ,   0.00000000     ,   1.00000000     ,
   0.00000000     ,   0.00000000     ,  -1.00000000     };    
    if(makeTSurf(s,nt,nv)) return -1;
    for(j=0;j<3*nt;j++)    s->tvec[j]=vt[j];
    for(j=0;j<3*nv;j++)  {
	s->vvec[j]=vv[j];
	s->nvec[j]=vn[j];
    }
    return 0;
}
 
int mkOctahedron(TSurf *s) {
    int nt=8,nv=6,j;
    int vt[]={          4 ,           0 ,           2 ,
           5 ,           2 ,           0 ,
           4 ,           3 ,           0 ,
           5 ,           0 ,           3 ,
           4 ,           2 ,           1 ,
           5 ,           1 ,           2 ,
           4 ,           1 ,           3 ,
           5 ,           3 ,           1 };
    Dbl vv[]={  1.00000000     ,   0.00000000     ,   0.00000000     ,
  -1.00000000     ,   0.00000000     ,   0.00000000     ,
   0.00000000     ,   1.00000000     ,   0.00000000     ,
   0.00000000     ,  -1.00000000     ,   0.00000000     ,
   0.00000000     ,   0.00000000     ,   1.00000000     ,
   0.00000000     ,   0.00000000     ,  -1.00000000 };
    Dbl vn[]={ 1.00000000     ,   0.00000000     ,   0.00000000     ,
  -1.00000000     ,   0.00000000     ,   0.00000000     ,
   0.00000000     ,   1.00000000     ,   0.00000000     ,
   0.00000000     ,  -1.00000000     ,   0.00000000     ,
   0.00000000     ,   0.00000000     ,   1.00000000     ,
   0.00000000     ,   0.00000000     ,  -1.00000000     };
    if(makeTSurf(s,nt,nv)) return -1;
    for(j=0;j<3*nt;j++)    s->tvec[j]=vt[j];
    for(j=0;j<3*nv;j++)  {
	s->vvec[j]=vv[j];
	s->nvec[j]=vn[j];
    }
    return 0;
}

int mkDodecahedron(TSurf *s) {
    int nt=60,nv=32,j;
    int vt[]={ 
          20 ,          7 ,           14 ,
          20 ,          14 ,          15 ,
          20 ,          15 ,          16 ,
          20 ,           16 ,          8 ,
	  20 ,           8 ,           7 ,
          21 ,           5 ,           6 ,
          21 ,           6 ,           7 ,
          21 ,           7 ,           8 ,
          21 ,           8 ,           9 ,
          21 ,           9 ,           5 ,
          22 ,           0 ,           4 ,
          22 ,           4 ,           3 ,
          22 ,           3 ,           2 ,
          22 ,           2 ,           1 ,
          22 ,           1 ,           0 ,
          23 ,           0 ,           1 ,
          23 ,           1 ,          11 ,
          23 ,          11 ,          10 ,
          23 ,          10 ,          19 ,
          23 ,          19 ,           0 ,
          24 ,           3 ,           4 ,
          24 ,           4 ,          17 ,
          24 ,          17 ,          16 ,
          24 ,          16 ,          15 ,
          24 ,          15 ,           3 ,
          25 ,           5 ,           9 ,
          25 ,           9 ,          18 ,
          25 ,          18 ,          19 ,
          25 ,          19 ,          10 ,
          25 ,          10 ,           5 ,
          26 ,           2 ,           3 ,
          26 ,           3 ,          15 ,
          26 ,          15 ,          14 ,
          26 ,          14 ,          13 ,
          26 ,          13 ,           2 ,
          27 ,           5 ,          10 ,
          27 ,          10 ,          11 ,
          27 ,          11 ,          12 ,
          27 ,          12 ,           6 ,
          27 ,           6 ,           5 ,
          28 ,           8 ,          16 ,
          28 ,          16 ,          17 ,
          28 ,          17 ,          18 ,
          28 ,          18 ,           9 ,
          28 ,           9 ,           8 ,
          29 ,           0 ,          19 ,
          29 ,          19 ,          18 ,
          29 ,          18 ,          17 ,
          29 ,          17 ,           4 ,
          29 ,           4 ,           0 ,
          30 ,           6 ,          12 ,
          30 ,          12 ,          13 ,
          30 ,          13 ,          14 ,
          30 ,          14 ,           7 ,
          30 ,           7 ,           6 ,
          31 ,           1 ,           2 ,
          31 ,           2 ,          13 ,
          31 ,          13 ,          12 ,
          31 ,          12 ,          11 ,
          31 ,          11 ,           1 };
    Dbl vv[]={ 0.356822103     , -0.934172392     ,   0.00000000     ,
 -0.356822103     , -0.934172392     ,   0.00000000     ,
 -0.577350259     , -0.577350259     ,  0.577350259     ,
   0.00000000     , -0.356822103     ,  0.934172392     ,
  0.577350259     , -0.577350259     ,  0.577350259     ,
   0.00000000     ,  0.356822103     , -0.934172392     ,
 -0.577350259     ,  0.577350259     , -0.577350259     ,
 -0.356822103     ,  0.934172392     ,   0.00000000     ,
  0.356822103     ,  0.934172392     ,   0.00000000     ,
  0.577350259     ,  0.577350259     , -0.577350259     ,
   0.00000000     , -0.356822103     , -0.934172392     ,
 -0.577350259     , -0.577350259     , -0.577350259     ,
 -0.934172392     ,   0.00000000     , -0.356822103     ,
 -0.934172392     ,   0.00000000     ,  0.356822103     ,
 -0.577350259     ,  0.577350259     ,  0.577350259     ,
   0.00000000     ,  0.356822103     ,  0.934172392     ,
  0.577350259     ,  0.577350259     ,  0.577350259     ,
  0.934172392     ,   0.00000000     ,  0.356822103     ,
  0.934172392     ,   0.00000000     , -0.356822103     ,
  0.577350259     , -0.577350259     , -0.577350259     ,
   0.00000000     ,  0.675973475     ,  0.417774558     ,
   0.00000000     ,  0.675973535     , -0.417774558     ,
   0.00000000     , -0.675973475     ,  0.417774558     ,
   0.00000000     , -0.675973475     , -0.417774558     ,
  0.417774558     ,   0.00000000     ,  0.675973475     ,
  0.417774558     ,   0.00000000     , -0.675973475     ,
 -0.417774558     ,   0.00000000     ,  0.675973475     ,
 -0.417774558     ,   0.00000000     , -0.675973475     ,
  0.675973535     ,  0.417774558     ,   0.00000000     ,
  0.675973535     , -0.417774558     ,   0.00000000     ,
 -0.675973475     ,  0.417774558     ,   0.00000000     ,
 -0.675973535     , -0.417774558     ,   0.00000000  };
    Dbl vn[]={  0.356822103     , -0.934172392     ,   0.00000000     ,
 -0.356822103     , -0.934172392     ,   0.00000000     ,
 -0.577350318     , -0.577350318     ,  0.577350318     ,
   0.00000000     , -0.356822103     ,  0.934172392     ,
  0.577350318     , -0.577350318     ,  0.577350318     ,
   0.00000000     ,  0.356822103     , -0.934172392     ,
 -0.577350318     ,  0.577350318     , -0.577350318     ,
 -0.356822103     ,  0.934172392     ,   0.00000000     ,
  0.356822103     ,  0.934172392     ,   0.00000000     ,
  0.577350318     ,  0.577350318     , -0.577350318     ,
   0.00000000     , -0.356822103     , -0.934172392     ,
 -0.577350318     , -0.577350318     , -0.577350318     ,
 -0.934172392     ,   0.00000000     , -0.356822103     ,
 -0.934172392     ,   0.00000000     ,  0.356822103     ,
 -0.577350318     ,  0.577350318     ,  0.577350318     ,
   0.00000000     ,  0.356822103     ,  0.934172392     ,
  0.577350318     ,  0.577350318     ,  0.577350318     ,
  0.934172392     ,   0.00000000     ,  0.356822103     ,
  0.934172392     ,   0.00000000     , -0.356822103     ,
  0.577350318     , -0.577350318     , -0.577350318     ,
   0.00000000     ,  0.850650787     ,  0.525731087     ,
   0.00000000     ,  0.850650787     , -0.525731027     ,
   0.00000000     , -0.850650787     ,  0.525731087     ,
   0.00000000     , -0.850650787     , -0.525731087     ,
  0.525731087     ,   0.00000000     ,  0.850650787     ,
  0.525731087     ,   0.00000000     , -0.850650787     ,
 -0.525731087     ,   0.00000000     ,  0.850650787     ,
 -0.525731087     ,   0.00000000     , -0.850650787     ,
  0.850650787     ,  0.525731027     ,   0.00000000     ,
  0.850650787     , -0.525731027     ,   0.00000000     ,
 -0.850650787     ,  0.525731087     ,   0.00000000     ,
 -0.850650787     , -0.525731027     ,   0.00000000    };
    if(makeTSurf(s,nt,nv)) return -1;
    for(j=0;j<3*nt;j++)    s->tvec[j]=vt[j];
    for(j=0;j<3*nv;j++)  {
	s->vvec[j]=vv[j];
	s->nvec[j]=vn[j];
    }
    return 0;
}
int mkIcosahedron(TSurf *s) {
    int nt=20,nv=12,j;
    int vt[]={         0 ,           4 ,           8 ,
           0 ,           8 ,           1 ,
           0 ,           1 ,          10 ,
           0 ,          10 ,           6 ,
           0 ,           6 ,           4 ,
           1 ,           8 ,           5 ,
           1 ,           5 ,           7 ,
           1 ,           7 ,          10 ,
           2 ,           3 ,           9 ,
           2 ,           9 ,           4 ,
           2 ,           4 ,           6 ,
           2 ,           6 ,          11 ,
           2 ,          11 ,           3 ,
           3 ,          11 ,           7 ,
           3 ,           7 ,           5 ,
           3 ,           5 ,           9 ,
           4 ,           9 ,           8 ,
           5 ,           8 ,           9 ,
           6 ,          10 ,          11 ,
           7 ,          11 ,          10 };
    Dbl vv[]={   0.00000000     ,  0.850650847     ,  0.525731146     ,
   0.00000000     ,  0.850650847     , -0.525731146     ,
   0.00000000     , -0.850650847     ,  0.525731146     ,
   0.00000000     , -0.850650847     , -0.525731146     ,
  0.525731146     ,   0.00000000     ,  0.850650847     ,
  0.525731146     ,   0.00000000     , -0.850650847     ,
 -0.525731146     ,   0.00000000     ,  0.850650847     ,
 -0.525731146     ,   0.00000000     , -0.850650847     ,
  0.850650847     ,  0.525731146     ,   0.00000000     ,
  0.850650847     , -0.525731146     ,   0.00000000     ,
 -0.850650847     ,  0.525731146     ,   0.00000000     ,
 -0.850650847     , -0.525731146     ,   0.00000000 };
    Dbl vn[]={   0.00000000     ,  0.850650847     ,  0.525731146     ,
   0.00000000     ,  0.850650847     , -0.525731146     ,
   0.00000000     , -0.850650847     ,  0.525731146     ,
   0.00000000     , -0.850650847     , -0.525731146     ,
  0.525731146     ,   0.00000000     ,  0.850650847     ,
  0.525731146     ,   0.00000000     , -0.850650847     ,
 -0.525731146     ,   0.00000000     ,  0.850650847     ,
 -0.525731146     ,   0.00000000     , -0.850650847     ,
  0.850650847     ,  0.525731146     ,   0.00000000     ,
  0.850650847     , -0.525731146     ,   0.00000000     ,
 -0.850650847     ,  0.525731146     ,   0.00000000     ,
 -0.850650847     , -0.525731146     ,   0.00000000 };
    if(makeTSurf(s,nt,nv)) return -1;
    for(j=0;j<3*nt;j++)    s->tvec[j]=vt[j];
    for(j=0;j<3*nv;j++)  {
	s->vvec[j]=vv[j];
	s->nvec[j]=vn[j];
    }
    return 0;
}

Dbl trgNorm(Vec3Dbl w,TSurf *s,int t) {
/* rerturns area of trg */    
    Vec3Dbl u,v; 
    int i,v0,v1,v2;
    Dbl a;
    v0=s->tvec[3*t+0];
    v1=s->tvec[3*t+1];
    v2=s->tvec[3*t+2];
    for(i=0;i<3;i++) {
	u[i]=s->vvec[3*v1+i]-s->vvec[3*v0+i];
	v[i]=s->vvec[3*v2+i]-s->vvec[3*v0+i];
    }
    cross3Dbl(w,u,v);
    a=norm3Dbl(w);
    if(a>MIN_NORM) {
	for(i=0;i<3;i++) w[i]=w[i]/a;
	return 0.5*a;
    } else {
	for(i=0;i<3;i++) w[i]=0.0f;
	return 0.0f;
    }
} 

int getTrgCon(int *ntc,CList *tcvec,TSurf *s) {
 /* for every veretex makes the list of triangles it is included */    
    int n0,k,n,i;
    if(ntc==0) return -1;
    if(tcvec==0) return -2; 
    if(s==0) return -3;
    for(k=0;k<s->nv;k++) {
	ntc[k]=0;
	for(i=0;i<MAX_CONNECT;i++) tcvec[k][i]=-1;
    }
    for(k=0;k<s->nt;k++) {
	for(i=0;i<3;i++) {
	    n0=s->tvec[3*k+i];
	    n=ntc[n0]++;
	    if(n<MAX_CONNECT)  
		tcvec[n0][n]=k;
	    else 
		return -4;
	}
    } 	
    return 0;
}


int getVrtCon(int *nvc,CList *vcvec,TSurf *s) {
/* for every vertex makes list of vertices (with larger numbers) it connects to */
 
    int n0,n1,n2,k,n,key,i;
    if(nvc==0) return -1;
    if(vcvec==0) return -2; 
    if(s==0) return -3;
    for(k=0;k<s->nv;k++) {
	nvc[k]=0;
	for(i=0;i<MAX_CONNECT;i++) vcvec[k][i]=-1;
    }
    for(k=0;k<s->nt;k++) {
	n0=s->tvec[3*k];
	n1=s->tvec[3*k+1];
	n2=s->tvec[3*k+2];
	if(n1>n0) {
	    key=1;
	    for(i=0;i<nvc[n0];i++) {//vert already in
		if(vcvec[n0][i]==n1) {
		    key=0;
		    break;
		}
	    }
	    if(key) {
		n=nvc[n0]++;
		if(n<MAX_CONNECT)   
		    vcvec[n0][n]=n1;
		else 
		    return -4;
	    }
	} else {
	    key=1;
	    for(i=0;i<nvc[n1];i++) {//vert already in
		if(vcvec[n1][i]==n0) {
		    key=0;
		    break;
		}
	    }
	    if(key) {
		n=nvc[n1]++;
		if(n<MAX_CONNECT)   
		    vcvec[n1][n]=n0;
		else 
		    return -4;
	    }
	} 
	if(n2>n0) {
	    key=1;
	    for(i=0;i<nvc[n0];i++) {//vert already in
		if(vcvec[n0][i]==n2) {
		    key=0;
		    break;
		}
	    }
	    if(key) {
		n=nvc[n0]++;
		if(n<MAX_CONNECT)   
		    vcvec[n0][n]=n2;
		else 
		    return -4;
	    }
	} else {
	    key=1;
	    for(i=0;i<nvc[n2];i++) {//vert already in
		if(vcvec[n2][i]==n0) {
		    key=0;
		    break;
		}
	    }
	    if(key) {
		n=nvc[n2]++;
		if(n<MAX_CONNECT)   
		    vcvec[n2][n]=n0;
		else 
		    return -4;
	    }
	} 

	if(n2>n1) {
	    key=1;
	    for(i=0;i<nvc[n1];i++) {//vert already in
		if(vcvec[n1][i]==n2) {
		    key=0;
		    break;
		}
	    }
	    if(key) {
		n=nvc[n1]++;
		if(n<MAX_CONNECT)   
		    vcvec[n1][n]=n2;
		else 
		    return -4;
	    }
	} else {
	    key=1;
	    for(i=0;i<nvc[n2];i++) {//vert already in
		if(vcvec[n2][i]==n1) {
		    key=0;
		    break;
		}
	    }
	    if(key) {
		n=nvc[n2]++;
		if(n<MAX_CONNECT)   
		    vcvec[n2][n]=n1;
		else 
		    return -4;
	    }
	} 


    } 	
    return 0;
}
int getVrtBetw(int *nvc,CList *vcvec,TSurf *s,int k1,int k2) {
    int nmin,j,k;
    Dbl xx,yy,zz,x,y,z,x1,x2,y1,y2,z1,z2,d;
    x1=s->vvec[3*k1];
    x2=s->vvec[3*k2];	
    y1=s->vvec[3*k1+1];
    y2=s->vvec[3*k2+1];
    z1=s->vvec[3*k1+2];
    z2=s->vvec[3*k2+2];
    xx=0.5*(x1+x2);
    yy=0.5*(y1+y2);
    zz=0.5*(z1+z2);
    nmin=k1>k2?k2:k1;
    for(j=0;j<nvc[nmin];j++) {
	k=vcvec[nmin][j];
	x=s->vvec[3*k];	
	y=s->vvec[3*k+1];
	z=s->vvec[3*k+2];
	d=(x-xx)*(x-xx)+(y-yy)*(y-yy)+(z-zz)*(z-zz);
	if(d<MIN_DIST) return k;
    }		
    return -1;
}




inline static int new_vert(TSurf *s,int n,int m,int k) {
    Dbl a,b;
    int j,ret;
    a=dot3Dbl(s->nvec+3*m,s->nvec+3*k);
    if(a>-1.0+MIN_NORM) {
	b=1.0/sqrt(2.0*(1.0+a));
	ret=0;
    } else {/* normal is set to 0 should not be here*/
	b=0.0;
	ret=1;
    }
    for(j=0;j<3;j++) {
	s->vvec[3*n+j]=0.5*(s->vvec[3*m+j]+s->vvec[3*k+j]);
	s->nvec[3*n+j]=b*(s->nvec[3*m+j]+s->nvec[3*k+j]);
    }
    return ret;
}

int refineTSurf2(TSurf *s,TSurf *sold,int *nvc,CList *vcvec) {
    int i,k,v0,v1,v2,ret,n,u0,u1,u2,j,m,v;
  //  Dbl x1,y1,z1,x2,y2,z2;
    n=sold->nv;
    for(j=0;j<sold->nv;j++) n+=nvc[j];/* new vertex count */
    if(makeTSurf(s,4*sold->nt,n)) return -1;
    for(i=0;i<3*sold->nv;i++) {
	s->vvec[i]=sold->vvec[i];
	s->nvec[i]=sold->nvec[i];
    }
    /* vertices */
    n=sold->nv;
    ret=0;
    for(j=0;j<sold->nv;j++) {
	m=nvc[j];
	for(k=0;k<m;k++) {
	    v=vcvec[j][k];
	    ret+=new_vert(s,n,j,v);
	    vcvec[j][k]=n;
	    n++;
	}
    }
    if(ret>0) return ret;

    /* triangles*/
    for(i=0;i<sold->nt;i++) {
	v0=sold->tvec[3*i+0];
	v1=sold->tvec[3*i+1];
	v2=sold->tvec[3*i+2];
	/* 0->1 */
	u0=getVrtBetw(nvc,vcvec,s,v0,v1);
	if(u0<0) return -2;
	/* 1->2 */
	u1=getVrtBetw(nvc,vcvec,s,v1,v2);
	if(u1<0) return -2;
	/* 2->0 */
	u2=getVrtBetw(nvc,vcvec,s,v2,v0);
	if(u2<0) return -2;
	s->tvec[3*(4*i+0)+0]=v0;
	s->tvec[3*(4*i+0)+1]=u0;
	s->tvec[3*(4*i+0)+2]=u2;	
	s->tvec[3*(4*i+1)+0]=v1;
	s->tvec[3*(4*i+1)+1]=u1;
	s->tvec[3*(4*i+1)+2]=u0;	
	s->tvec[3*(4*i+2)+0]=v2;
	s->tvec[3*(4*i+2)+1]=u2;
	s->tvec[3*(4*i+2)+2]=u1;	
	s->tvec[3*(4*i+3)+0]=u0;
	s->tvec[3*(4*i+3)+1]=u1;
	s->tvec[3*(4*i+3)+2]=u2;		
    }
    return ret;
}

int checkNormals(Dbl * d,TSurf *s) {
    int n0,n1,n2;
    int t,i;
    Vec3Dbl w1,w2;
    Dbl a;
    for(t=0;t<s->nt;t++) {
	a= trgNorm(w1,s,t);
	if(a==0.0f) return -1;
	n0=s->tvec[3*t+0];
	n1=s->tvec[3*t+1];
	n2=s->tvec[3*t+2];
	w2[0]=0.0;w2[1]=0.0;w2[2]=0.0;
	for(i=0;i<3;i++) w2[i]=(s->nvec[3*n0+i]+s->nvec[3*n1+i]+s->nvec[3*n2+i]);
	a=norm3Dbl(w2);
	if(a<MIN_NORM) return t;
	for(i=0;i<3;i++) w2[i]=w2[i]/a;
	d[t]=dot3Dbl(w1,w2);
    }
    return 0;
}
    
void mkUnitSphere(TSurf *s) {
    Vec3Dbl cnt,w1;
    int j,i;
    Dbl a=1.0/s->nv;
    Dbl b;
    for(i=0;i<3;i++) cnt[i]=0.0;
    for(j=0;j<s->nv;j++) {
	for(i=0;i<3;i++) cnt[i]+=s->vvec[3*j+i];
    }
    for(i=0;i<3;i++) cnt[i]*=a;
    for(j=0;j<s->nv;j++) {
	for(i=0;i<3;i++) w1[i]=s->vvec[3*j+i]-cnt[i];
	b=1.0/norm3Dbl(w1);
	for(i=0;i<3;i++) {
	    w1[i]*=b;
	    s->nvec[3*j+i]=w1[i];
	    s->vvec[3*j+i]=w1[i]+cnt[i];
	}
    }
}
Dbl dist(TSurf *s,int n1,int n2) {
    Dbl x,y,z;
    x=s->vvec[3*n1+0]-s->vvec[3*n2+0];
    y=s->vvec[3*n1+1]-s->vvec[3*n2+1];
    z=s->vvec[3*n1+2]-s->vvec[3*n2+2];
    return sqrt(x*x+y*y+z*z);
}

void printStat(TSurf *s) {
    Dbl a,amin,amax,amean;//trg areas
    Dbl d,dmin,dmax,dmean;//distances between vertices
    Dbl b,bmin,bmax,bmean;//diff between trinag normal and mean normal
    int n0,n1,n2;
    int t,i;
    Vec3Dbl w1,w2;
    amin=MATH_LARGE_NUMBER_DBL;amax=0.0;amean=0.0;
    bmin=MATH_LARGE_NUMBER_DBL;bmax=0.0;bmean=0.0;
    dmin=MATH_LARGE_NUMBER_DBL;dmax=0.0;dmean=0.0;
    for(t=0;t<s->nt;t++) {
	a= trgNorm(w1,s,t);
	if(a>amax) amax=a;
	if(a<amin) amin=a;
	amean+=a;
	n0=s->tvec[3*t+0];
	n1=s->tvec[3*t+1];
	n2=s->tvec[3*t+2];
	w2[0]=0.0;w2[1]=0.0;w2[2]=0.0;
	for(i=0;i<3;i++) w2[i]=(s->nvec[3*n0+i]+s->nvec[3*n1+i]+s->nvec[3*n2+i]);
	a=norm3Dbl(w2);
	for(i=0;i<3;i++) w2[i]=w2[i]/a;
	b=1.0-dot3Dbl(w1,w2);
	if(b>bmax) bmax=b;
	if(b<bmin) bmin=b;
	bmean+=b;
	d=dist(s,n0,n1);
	if(d>dmax) dmax=d;
	if(d<dmin) dmin=d;
	dmean+=d;
	d=dist(s,n1,n2);
	if(d>dmax) dmax=d;
	if(d<dmin) dmin=d;
	dmean+=d;
	d=dist(s,n2,n0);
	if(d>dmax) dmax=d;
	if(d<dmin) dmin=d;
	dmean+=d;
    }
    printf("distance(min,max,mean): %e %e %e\n",dmin,dmax,dmean/(3*s->nt));
    printf("area (min,max,mean): %e %e %e\n",amin,amax,amean/s->nt);
    printf("norm diff (min,max,mean): %e %e %e\n",bmin,bmax,bmean/s->nt);
 
}

    

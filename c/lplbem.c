/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "lplbem.h"

int mkQtot1(Dbl *q,int *ntc,CList *tcvec,TSurf *s) {
    Vec3Dbl w;
    Dbl a;
    int i,k;
    if(!ntc) return -1;
    if(!tcvec) return -2;
    for(i=0;i<s->nv;i++) {
	a=0.0;
	for(k=0;k<ntc[i];k++) {
	    a+=trgNorm(w,s,tcvec[i][k]);
	}
	q[i]=a/6.0;
    }
    return 0;
}

int mkQtot0(Dbl *q,TSurf *s) {
    Vec3Dbl w;
    int k;
    for(k=0;k<s->nt;k++) {
	q[k]=0.5*trgNorm(w,s,k);	
    }
    return 0;
}
int mkCenters(Dbl *c,TSurf *s) {
    int j,k,n0,n1,n2;
    if(!c) return -1;
    for(j=0;j<s->nt;j++) {
	n0=s->tvec[3*j+0];
	n1=s->tvec[3*j+1];
	n2=s->tvec[3*j+2];
	for(k=0;k<3;k++) c[3*j+k]=(s->vvec[3*n0+k]+s->vvec[3*n1+k]+s->vvec[3*n2+k])/3.0f;
    }
    return 0;
}

int  mkSAMat1(Dbl *lm,int *ntc,CList *tcvec,TSurf *s,TrgPot1 tp1){
    int i,j,k,v0,v1,v2,n,m;
    Dbl *lmc;
    Dbl q0,q1,q2;
    Dbl *cpvec;
    n=s->nv;
    m=s->nv;
    cpvec=s->vvec;
    for(i=0;i<n*n;i++) lm[i]=0.0;
    for(i=0;i<n;i++) {
	lmc=lm+i*m;
	for(k=0;k<ntc[i];k++) {
	    v0=s->tvec[3*tcvec[i][k]+0];
	    v1=s->tvec[3*tcvec[i][k]+1];
	    v2=s->tvec[3*tcvec[i][k]+2];
	    if(i==v0) {
		q0=1.0;q1=0.0;q2=0.0;
	    } else if( i==v1) {
		q1=1.0;q0=0.0;q2=0.0;
	    } else {
		q2=1.0;q1=0.0;q0=0.0;
	    }
	    for(j=0;j<m;j++) {
		lmc[j]+=tp1(cpvec+3*j,s->vvec+3*v0,s->vvec+3*v1,s->vvec+3*v2,q0,q1,q2);
	    }
	}
    }
    return 0;
}

int  mkSAMat0(Dbl *lm,Dbl *cpvec,TSurf *s,TrgPot0 tp0){
    int i,j,v0,v1,v2,n,m;
    Dbl *lmc;
    n=s->nt;
    m=s->nt;
    for(i=0;i<n*n;i++) lm[i]=0.0;
    for(i=0;i<n;i++) {
	lmc=lm+i*m;
	v0=s->tvec[3*i+0];
	v1=s->tvec[3*i+1];
	v2=s->tvec[3*i+2];
	for(j=0;j<m;j++) {
	    lmc[j]=tp0(cpvec+3*j,s->vvec+3*v0,s->vvec+3*v1,s->vvec+3*v2);
	}
    }
    return 0;
}
int  mkSAMat1Tst(Dbl *lm,int *ntc,CList *tcvec,TSurf *s,TrgPot1 tp1){
    int i,j,k,v0,v1,v2,n,m;
    Dbl *lmc;
    Dbl q0,q1,q2,a;
    Dbl *cpvec;
    n=s->nv;
    m=s->nv;
    cpvec=s->vvec;
    //for(i=0;i<n*n;i++) lm[i]=0.0;
    for(j=0;j<m;j++) {
	for(i=0;i<n;i++) {
	    lmc=lm+i*m;
	    a=0.0;
	    for(k=0;k<ntc[i];k++) {
		v0=s->tvec[3*tcvec[i][k]+0];
		v1=s->tvec[3*tcvec[i][k]+1];
		v2=s->tvec[3*tcvec[i][k]+2];
		if(i==v0) {
		    q0=1.0;q1=0.0;q2=0.0;
		} else if( i==v1) {
		    q1=1.0;q0=0.0;q2=0.0;
		} else {
		    q2=1.0;q1=0.0;q0=0.0;
		}
		a+=tp1(cpvec+3*j,s->vvec+3*v0,s->vvec+3*v1,s->vvec+3*v2,q0,q1,q2);
	    }
	    lmc[j]=a;
	}
    }
    return 0;
}


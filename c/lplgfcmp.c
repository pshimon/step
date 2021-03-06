/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "lplbem.h"
#include "timers.h"
double trgint(double x,double y,double z,
             double x0,double y0,double z0,double q0,
             double x1,double y1,double z1,double q1,
             double x2,double y2,double z2,double q2);

void vt11(Dbl *lm,TSurf *s) {
    int i,j,v0,v1,v2,n,t;
    Dbl *lmc;
    Dbl q0,q1,q2;
    Dbl *cpvec;
    n=s->nv;
    t=s->nt;
    cpvec=s->vvec;
    q0=1.0;q1=0.0;q2=0.0;
    for(i=0;i<n*t;i++) lm[i]=0.0;
    for(i=0;i<t;i++) {
	lmc=lm+i*n;
	v0=s->tvec[3*i+0];
	v1=s->tvec[3*i+1];
	v2=s->tvec[3*i+2];
	for(j=0;j<n;j++) {
	    lmc[j]+=lplGfL1(cpvec+3*j,s->vvec+3*v0,s->vvec+3*v1,s->vvec+3*v2,q0,q1,q2);
	}
    }
}

void vt41(Dbl *lm,TSurf *s) {
    int i,j,v0,v1,v2,n,t;
    Dbl *lmc;
    Dbl q0,q1,q2;
    Dbl *cpvec;
    n=s->nv;
    t=s->nt;
    cpvec=s->vvec;
    q0=1.0;q1=1.0;q2=1.0;
    for(i=0;i<n*t;i++) lm[i]=0.0;
    for(i=0;i<t;i++) {
	lmc=lm+i*n;
	v0=s->tvec[3*i+0];
	v1=s->tvec[3*i+1];
	v2=s->tvec[3*i+2];
	for(j=0;j<n;j++) {
	    lmc[j]+=lplGfL1(cpvec+3*j,s->vvec+3*v0,s->vvec+3*v1,s->vvec+3*v2,q0,q1,q2);
	}
    }
}
void vt43(Dbl *lm,TSurf *s) {
    int i,j,v0,v1,v2,n,t;
    Dbl *lmc;
    Dbl *cpvec;
    n=s->nv;
    t=s->nt;
    cpvec=s->vvec;
    
    for(i=0;i<n*t;i++) lm[i]=0.0;
    for(i=0;i<t;i++) {
	lmc=lm+i*n;
	v0=s->tvec[3*i+0];
	v1=s->tvec[3*i+1];
	v2=s->tvec[3*i+2];
	for(j=0;j<n;j++) {
	    lmc[j]+=lplGfC1(cpvec+3*j,s->vvec+3*v0,s->vvec+3*v1,s->vvec+3*v2);
	}
    }
}
void vt44(Dbl *lm,TSurf *s) {
    int i,j,v0,v1,v2,n,t;
    Dbl *lmc;
    Dbl *cpvec;
    n=s->nv;
    t=s->nt;
    cpvec=s->vvec;

    for(i=0;i<n*t;i++) lm[i]=0.0;
    for(i=0;i<t;i++) {
	lmc=lm+i*n;
	v0=s->tvec[3*i+0];
	v1=s->tvec[3*i+1];
	v2=s->tvec[3*i+2];
	for(j=0;j<n;j++) {
	    lmc[j]+=lplGfC2(cpvec+3*j,s->vvec+3*v0,s->vvec+3*v1,s->vvec+3*v2);
	}
    }
}

void vt12(Dbl *lm,TSurf *s) {
    int i,j,v0,v1,v2,n,t;
    Dbl *lmc;
    Dbl q0,q1,q2;
    Dbl *cpvec;
    n=s->nv;
    t=s->nt;
    cpvec=s->vvec;
    q0=1.0;q1=0.0;q2=0.0;
    for(i=0;i<n*t;i++) lm[i]=0.0;
    for(i=0;i<t;i++) {
	lmc=lm+i*n;
	v0=s->tvec[3*i+0];
	v1=s->tvec[3*i+1];
	v2=s->tvec[3*i+2];
	for(j=0;j<n;j++) {
	    lmc[j]+=lplGfL2(cpvec+3*j,s->vvec+3*v0,s->vvec+3*v1,s->vvec+3*v2,q0,q1,q2);
	}
    }
}

void vt10(Dbl *lm,TSurf *s) {
    int i,j,v0,v1,v2,n,t;
    Dbl *lmc;
    Dbl q0,q1,q2;
    Dbl *cpvec;
    Dbl x,y,z,x0,y0,z0,x1,y1,z1,x2,y2,z2;
    n=s->nv;
    t=s->nt;
    cpvec=s->vvec;
    q0=1.0;q1=0.0;q2=0.0;
    for(i=0;i<n*t;i++) lm[i]=0.0;
    for(i=0;i<t;i++) {
	lmc=lm+i*n;
	v0=s->tvec[3*i+0];
	x0=s->vvec[v0*3+0];
	y0=s->vvec[v0*3+1];
	z0=s->vvec[v0*3+2];

	v1=s->tvec[3*i+1];
	x1=s->vvec[v1*3+0];
	y1=s->vvec[v1*3+1];
	z1=s->vvec[v1*3+2];

	v2=s->tvec[3*i+2];
	x2=s->vvec[v2*3+0];
	y2=s->vvec[v2*3+1];
	z2=s->vvec[v2*3+2];

	for(j=0;j<n;j++) {
	    x=cpvec[3*j+0];
	    y=cpvec[3*j+1];
	    z=cpvec[3*j+2];
	    lmc[j]+=trgint(x,y,z,x0,y0,z0,q0,x1,y1,z1,q1,x2,y2,z2,q2);
	}
    }
}
void vt40(Dbl *lm,TSurf *s) {
    int i,j,v0,v1,v2,n,t;
    Dbl *lmc;
    Dbl q0,q1,q2;
    Dbl *cpvec;
    Dbl x,y,z,x0,y0,z0,x1,y1,z1,x2,y2,z2;
    n=s->nv;
    t=s->nt;
    cpvec=s->vvec;
    q0=1.0;q1=1.0;q2=1.0;
    for(i=0;i<n*t;i++) lm[i]=0.0;
    for(i=0;i<t;i++) {
	lmc=lm+i*n;
	v0=s->tvec[3*i+0];
	x0=s->vvec[v0*3+0];
	y0=s->vvec[v0*3+1];
	z0=s->vvec[v0*3+2];

	v1=s->tvec[3*i+1];
	x1=s->vvec[v1*3+0];
	y1=s->vvec[v1*3+1];
	z1=s->vvec[v1*3+2];

	v2=s->tvec[3*i+2];
	x2=s->vvec[v2*3+0];
	y2=s->vvec[v2*3+1];
	z2=s->vvec[v2*3+2];

	for(j=0;j<n;j++) {
	    x=cpvec[3*j+0];
	    y=cpvec[3*j+1];
	    z=cpvec[3*j+2];
	    lmc[j]+=trgint(x,y,z,x0,y0,z0,q0,x1,y1,z1,q1,x2,y2,z2,q2);
	}
    }
}
int main(int argc,char * argv[]) {
    TSurf s;   
    double time_start,time_stop;
    int ret,n,t;
    Dbl * cnt=0;
    Int * ntc=0;
    CList *tcvec=0;
    Dbl *lm=0; 
    char str[100];
  

    if(argc!=3) {
	fprintf(stderr,"usage: %s surf lbl\n",argv[0]);
	exit(1);
    }
    initTSurf(&s);
    ret=readTSurf(&s,argv[1]);
    if(ret) {
	fprintf(stderr,"readTSurf returns %d\n",ret);
	exit(1);
    }
    n=s.nv;
    t=s.nt;
    printf("nt=%d, nv=%d\n",t,n);
    cnt=ALLOC_MEM(Dbl,3*t);
    ret=mkCenters(cnt,&s);
    if(ret) {
	fprintf(stderr,"mkCenters returns %d\n",ret);
	exit(1);
    }
    lm=ALLOC_MEM(Dbl,n*t);
    time_start=cpuClock();
    vt10(lm,&s);
    time_stop=cpuClock();
    printf("vt10 takes %e s\n",time_stop-time_start);
    sprintf(str,"%s-vt10k.bin",argv[2]);
    ret=write2DataBufDbl(lm,n,t,str);
    FREE_MEM(lm);
    FREE_MEM(cnt);
    FREE_MEM(tcvec);
    FREE_MEM(ntc);
    cleanTSurf(&s);
    return 0;
}



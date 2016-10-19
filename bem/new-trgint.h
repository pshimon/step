/******************************************************************
* Author: Shimon Panfil                                           *
* Copyright (c) Shimon Panfil Industrial Physics and Simulations  *
* http://industrialphys.com                                       *
******************************************************************/ 
#ifndef TRGINT_H
#define TRGINT_H
double trgint(double x,double y,double z,
             double x0,double y0,double z0,double q0,
             double x1,double y1,double z1,double q1,
             double x2,double y2,double z2,double q2);
double int_trg(double x,double y,double z,
             double x0,double y0,double z0,double q0,
             double x1,double y1,double z1,double q1,
             double x2,double y2,double z2,double q2);
double int_trg1(double x,double y,double z,
             double x0,double y0,double z0,double q0,
             double x1,double y1,double z1,double q1,
             double x2,double y2,double z2,double q2);

double int_trg2(double x,double y,double z,
             double x0,double y0,double z0,double q0,
             double x1,double y1,double z1,double q1,
             double x2,double y2,double z2,double q2);

inline static double trgint1(double x,double y,double z,const double * trg) {	
    
    //triangle:
    double x0=trg[0],y0=trg[1],z0=trg[2],q0=trg[3];//vertex 0
    double x1=trg[4],y1=trg[5],z1=trg[6],q1=trg[7];//vertex 1
    double x2=trg[8],y2=trg[9],z2=trg[10],q2=trg[11];//vertex 2
    return trgint(x,y,z,x0,y0,z0,q0,x1,y1,z1,q1,x2,y2,z2,q2);

}
inline static double trgint2(double x,double y,double z,const double * trg) {	
    
    //triangle:
    double x0=trg[0],y0=trg[1],z0=trg[2],q0=trg[3];//vertex 0
    double x1=trg[4],y1=trg[5],z1=trg[6],q1=trg[7];//vertex 1
    double x2=trg[8],y2=trg[9],z2=trg[10],q2=trg[11];//vertex 2
    return int_trg(x,y,z,x0,y0,z0,q0,x1,y1,z1,q1,x2,y2,z2,q2);
}
inline static double trgint3(double x,double y,double z,const double * trg) {	
    
    //triangle:
    double x0=trg[0],y0=trg[1],z0=trg[2],q0=trg[3];//vertex 0
    double x1=trg[4],y1=trg[5],z1=trg[6],q1=trg[7];//vertex 1
    double x2=trg[8],y2=trg[9],z2=trg[10],q2=trg[11];//vertex 2
    return int_trg1(x,y,z,x0,y0,z0,q0,x1,y1,z1,q1,x2,y2,z2,q2);
}
inline static double trgint4(double x,double y,double z,const double * trg) {	
    
    //triangle:
    double x0=trg[0],y0=trg[1],z0=trg[2],q0=trg[3];//vertex 0
    double x1=trg[4],y1=trg[5],z1=trg[6],q1=trg[7];//vertex 1
    double x2=trg[8],y2=trg[9],z2=trg[10],q2=trg[11];//vertex 2
    return int_trg2(x,y,z,x0,y0,z0,q0,x1,y1,z1,q1,x2,y2,z2,q2);
}

#define CHARGE_TRG_LENGTH 12
inline static double potCTV1(double x,double y,double z,const double *ctv,int L) {
    double res=0.0;
    int j;
    for(j=0;j<L;j+=CHARGE_TRG_LENGTH) res+=trgint1(x,y,z,ctv+j);
    return res;
};
inline static double potCTV2(double x,double y,double z,const double *ctv,int L) {
    double res=0.0;
    int j;
    for(j=0;j<L;j+=CHARGE_TRG_LENGTH) res+=trgint2(x,y,z,ctv+j);
    return res;
};
inline static double potCTV3(double x,double y,double z,const double *ctv,int L) {
    double res=0.0;
    int j;
    for(j=0;j<L;j+=CHARGE_TRG_LENGTH) res+=trgint3(x,y,z,ctv+j);
    return res;
};
inline static double potCTV4(double x,double y,double z,const double *ctv,int L) {
    double res=0.0;
    int j;
    for(j=0;j<L;j+=CHARGE_TRG_LENGTH) res+=trgint4(x,y,z,ctv+j);
    return res;
};

#endif

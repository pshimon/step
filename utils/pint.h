/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#ifndef PINT_H
#define PINT_H
#include "func.h"
#include "surf.h"
#include "trgint.h"
#define ZZERO 1.0e-12
//#define ZZERO 1.0e-16
inline static dbl pot_trg0(v3_flt pnt,patch_flt patch) {/* unit charge density */
    dbl x,y,z,x0,y0,z0,x1,y1,z1,x2,y2,z2;
    x=pnt[0];y=pnt[1];z=pnt[2];
    x0=patch[0][0];y0=patch[0][1];z0=patch[0][2];
    x1=patch[1][0]+x0;y1=patch[1][1]+y0;z1=patch[1][2]+z0;
    x2=patch[2][0]+x0;y2=patch[2][1]+y0;z2=patch[2][2]+z0;
    return trgint0(x,y,z,x0,y0,z0,x1,y1,z1,x2,y2,z2);
}
inline static dbl pot_trg1(v3_flt pnt,patch_flt patch) {/* unit charge density */
    dbl x,y,z,x0,y0,z0,x1,y1,z1,x2,y2,z2,q=1.0;
    x=pnt[0];y=pnt[1];z=pnt[2];
    x0=patch[0][0];y0=patch[0][1];z0=patch[0][2];
    x1=patch[1][0]+x0;y1=patch[1][1]+y0;z1=patch[1][2]+z0;
    x2=patch[2][0]+x0;y2=patch[2][1]+y0;z2=patch[2][2]+z0;
    return trgint1(x,y,z,x0,y0,z0,x1,y1,z1,x2,y2,z2,q,q,q);
}
    
dbl pot_patch_t(v3_flt pnt,patch_flt patch);/* trg, unit charge density */
dbl pot_patch_p(v3_flt pnt,patch_flt patch);/* parallelogramm, unit charge density
not used actually*/
dbl pot_patch_r(v3_flt pnt,patch_flt patch);/* rectangular, unit charge density */
inline static dbl rtpot(v3_flt pnt,patch_flt *parr,dbl* carr,int rpnum,int tpnum ) {
    dbl p;
    int j;
    p=0.0;
    for(j=0;j<rpnum;j++) 
	p+=pot_patch_r(pnt,parr[j])*carr[j];
    for(j=0;j<tpnum;j++) 
	p+=pot_patch_t(pnt,parr[j+rpnum])*carr[j+rpnum];
    return p;
};  


#endif /* PINT_H */

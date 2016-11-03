/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#ifndef LPLGF_H
#define LPLGF_H
#include "databuf.h"
/* Green functions for Laplace equation integrated
 * over triangle with linear charge distribution
 * values q0,q1,q2 in vertices
 */

/* linear charge distribution */
double trgInt1(double x,double y,double z,
             double x0,double y0,double z0,
             double x1,double y1,double z1,
             double x2,double y2,double z2,
	     double q0,double q1,double q2);
double intTrg(double x,double y,double z,
             double x0,double y0,double z0,double q0,
             double x1,double y1,double z1,double q1,
             double x2,double y2,double z2,double q2);

/* Green functions for Laplace equation integrated
 * over triangle with constant unit charge distribution
 */
double lplGfC(float dst[3],float vrt0[3],float vrt1[3],float vrt2[3]);
typedef Flt patch_flt[3][3];//temporary
double potPatchT(Vec3Flt pnt,patch_flt patch);/* trg, unit charge density */
double potPatchP(Vec3Flt pnt,patch_flt patch);/* parallelogramm, unit charge density
not used actually*/
double potPatchR(Vec3Flt pnt,patch_flt patch);/* rectangular, unit charge density */

#endif

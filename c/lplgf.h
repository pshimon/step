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

double trgInt1(double x0,double y0,double z0,
             double x1,double y1,double z1,
             double x2,double y2,double z2,
	     double q0,double q1,double q2);

double trgInt(double x,double y,double z,
             double x0,double y0,double z0,double q0,
             double x1,double y1,double z1,double q1,
             double x2,double y2,double z2,double q2);


typedef Dbl (*TrgPot1) (Vec3Dbl,Vec3Dbl,Vec3Dbl,Vec3Dbl,Dbl,Dbl,Dbl); 

Dbl lplGfL1 (Vec3Dbl dst ,Vec3Dbl vrt0,Vec3Dbl vrt1,Vec3Dbl vrt2,Dbl q0,Dbl q1,Dbl q2);
/*double intTrg(double x,double y,double z,
             double x0,double y0,double z0,double q0,
             double x1,double y1,double z1,double q1,
             double x2,double y2,double z2,double q2);*/
Dbl lplGfL2 (Vec3Dbl dst ,Vec3Dbl vrt0,Vec3Dbl vrt1,Vec3Dbl vrt2,Dbl q0,Dbl q1,Dbl q2);
/* Green functions for Laplace equation integrated
 * over triangle with constant unit charge distribution
 */
typedef Dbl (*TrgPot0) (Vec3Dbl,Vec3Dbl,Vec3Dbl,Vec3Dbl);
Dbl lplGfC1(Vec3Dbl dst,Vec3Dbl vrt0,Vec3Dbl vrt1,Vec3Dbl vrt2);
Dbl lplGfC2(Vec3Dbl dst,Vec3Dbl vrt0,Vec3Dbl vrt1,Vec3Dbl vrt2);

typedef Dbl PatchDbl[3][3];//temporary
double potPatchT(Vec3Dbl pnt,PatchDbl patch);/* trg, unit charge density */
double potPatchP(Vec3Dbl pnt,PatchDbl patch);/* parallelogramm, unit charge density
not used actually*/
double potPatchR(Vec3Dbl pnt,PatchDbl patch);/* rectangular, unit charge density */

#endif

/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#ifndef C_LPLGF_CC_H
#define C_LPLGF_CC_H
/* Green functions for Laplace equation integrated
 * over triangle with constant unit charge distribution
 */

double trgint0(double x,double y,double z,
             double x0,double y0,double z0,
             double x1,double y1,double z1,
             double x2,double y2,double z2);
dbl pot_patch_t(v3_flt pnt,patch_flt patch);/* trg, unit charge density */
dbl pot_patch_p(v3_flt pnt,patch_flt patch);/* parallelogramm, unit charge density
not used actually*/
dbl pot_patch_r(v3_flt pnt,patch_flt patch);/* rectangular, unit charge density */

#endif

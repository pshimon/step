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
double c_lplgf_cc1(float dst[3],float vrt0[3],float vrt1[3],float vrt2[3]);
double pot_patch_t(v3_flt pnt,patch_flt patch);/* trg, unit charge density */
double pot_patch_p(v3_flt pnt,patch_flt patch);/* parallelogramm, unit charge density
not used actually*/
double pot_patch_r(v3_flt pnt,patch_flt patch);/* rectangular, unit charge density */

#endif

/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#ifndef C_LPLGF_LC_H
#define C_LPLGF_LC_H
/* Green functions for Laplace equation integrated
 * over triangle with linear charge distribution
 * values q0,q1,q2 in vertices
 */

/* linear charge distribution */
double trgint1(double x,double y,double z,
             double x0,double y0,double z0,
             double x1,double y1,double z1,
             double x2,double y2,double z2,
	     double q0,double q1,double q2);
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

#endif

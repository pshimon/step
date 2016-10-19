/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
/*
 * $Log: trgint.h,v $
 * Revision 1.1  2011/01/10 10:06:45  shimon
 * Initial revision
 *
 */
#ifndef TRGINT_H
#define TRGINT_H

/* linear charge distribution */
double trgint1(double x,double y,double z,
             double x0,double y0,double z0,
             double x1,double y1,double z1,
             double x2,double y2,double z2,
	     double q0,double q1,double q2);
/* constant unit charge density q0=q1=q2=1.0 */
double trgint0(double x,double y,double z,
             double x0,double y0,double z0,
             double x1,double y1,double z1,
             double x2,double y2,double z2);

#endif

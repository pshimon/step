/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef LPLBEM_H
#define LPLBEM_H
#include "tsurf.h"
#include "lplgf.h"
/* ntc, tcvec should be updated before using this function*/
int mkQtot1(Dbl *q,int *ntc,CList *tcvec,TSurf *s);//linear charge distribution
int mkQtot0(Dbl *q,TSurf *s);//constant charge distribution
int mkCenters(Dbl *q,TSurf *s);//centers of triangles


/* self action matrices */
int mkSAMat1(Dbl *lm,int *ntc,CList *tcvec,TSurf *s,TrgPot1 tp);
int mkSAMat0(Dbl *lm,Dbl *cpvec,TSurf *s,TrgPot0 tp);
#endif

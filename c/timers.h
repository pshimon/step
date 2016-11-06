/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef TIMERS_H
#define TIMERS_H
/*High-resolution per-process timer from the CPU 
 * on linux -lrt library should be linked*/

double cpuClock();
void getWallTime(double* wcTime);
void getWallTime_(double* wcTime);

#endif

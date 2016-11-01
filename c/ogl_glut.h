/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef OGL_GLUT_H
#define OGL_GLUT_H
#include <GL/freeglut.h>
/* returns window handle, 0 -failure */
int glut_init_window(int * argc,char * argv[]);
/* interface */
void resize_cb(int w,int h);
void display_cb();
void timer_cb(int v);
void idle_cb();
void clean_cb(void);
void keyboard_cb(unsigned char key,int x,int y);
void special_cb(int key,int x,int y);

#endif

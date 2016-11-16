/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef OGL_H
#define OGL_H
#include  <GL/glew.h>
#define MAX_FILE_NAME_LENGTH 256
#include "databuf.h"

#define BUFFER_OFFSET(offs) ((void *) (offs))

char * readBytes(char* filename);
void oglInit3d();
/* returns program id on success 0 on error */
GLuint  oglMakeProg(const GLchar * vs,const GLchar * fs);
void oglDeleteProg(GLuint pid);
void oglDeleteBuffs(GLuint vaoid[],GLuint vboid[]);
void oglPrintLog(GLuint sid);

void oglSetUniformInt(GLuint pid,const char * name,int val);
/* float */
void oglMakeBuffsFlt(GLuint vaoid[],GLuint vboid[],int nv,int nt,float * vvec,float * nvec,int * tvec);
void oglSetUniformM4Flt(GLuint pid,const char * name,Mat4Flt val);
void oglSetUniformM3Flt(GLuint pid,const char * name,Mat3Flt val);
void oglSetUniformV3Flt(GLuint pid,const char * name,Flt * val);
void oglSetUniformV4Flt(GLuint pid,const char * name,Flt * val);
void oglSetUniformFlt(GLuint pid,const char * name,Flt val);

#endif

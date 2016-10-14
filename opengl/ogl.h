/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#ifndef OGL_H
#define OGL_H
#include  <GL/glew.h>
#define MAX_FILE_NAME_LENGTH 256
#include "my_cdefs.h"

#define BUFFER_OFFSET(offs) ((void *) (offs))

char * read_bytes(char* filename);

/* returns program id on success 0 on error */
GLuint  ogl_make_prog(const GLchar * vs,const GLchar * fs);
void ogl_delete_prog (GLuint pid);
void ogl_init_3d();
void ogl_make_buffs(GLuint vaoid[],GLuint vboid[],int nv,int nt,float * vvec,float * nvec,int * tvec);
void ogl_delete_buffs(GLuint vaoid[],GLuint vboid[]);
void ogl_set_uniform_m4(GLuint pid,const char * name,t_m4 val);
void ogl_print_log(GLuint sid);
void ogl_set_uniform_m3(GLuint pid,const char * name,t_m3 val);
void ogl_set_uniform_v3(GLuint pid,const char * name,float * val);
void ogl_set_uniform_v4(GLuint pid,const char * name,float * val);
void ogl_set_uniform_f32(GLuint pid,const char * name,float val);
void ogl_set_uniform_i32(GLuint pid,const char * name,int val);
#endif

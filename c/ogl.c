/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/

#include "ogl.h"
void oglInit3d() {
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glFrontFace(GL_CCW);
}
/* shaders are described as strings */
GLuint  oglMakeProg(const GLchar * vs,const GLchar * fs) {
    GLuint vid,fid,pid;
    GLint flag;
    GLenum err = glGetError();//reset glGetError()

    vid = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vid, 1, &vs, NULL);
    glCompileShader(vid);
    glGetShaderiv(vid, GL_COMPILE_STATUS, &flag);
    if (!flag ) {
	fprintf(stderr,"vertex shader compilation error \n");
	oglPrintLog(vid);
	return 0;
    }

    fid = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fid, 1, &fs, NULL);
    glCompileShader(fid);
    glGetShaderiv(fid, GL_COMPILE_STATUS, &flag);
    if (!flag ) {
	fprintf(stderr,"fragment shader compilation error ");
	oglPrintLog(fid);
	return 0;
    }

    pid = glCreateProgram();
    glAttachShader(pid,vid);
    glAttachShader(pid,fid);
    glLinkProgram(pid);
    glUseProgram(pid);

    err = glGetError();
    if (err != GL_NO_ERROR) return 0;
    return pid;
}

void oglDeleteProg (GLuint pid) {
    GLint ns = 0; 
    GLuint * sid; 
    int i;
    if(pid == 0) return;
    glUseProgram(0);
    glGetProgramiv(pid, GL_ATTACHED_SHADERS, &ns);
    sid = ALLOC_MEM(GLuint,ns);
    glGetAttachedShaders(pid, ns, NULL, sid);
    for (i = 0; i < ns; i++)
    glDeleteShader(sid[i]);
    glDeleteProgram (pid);
    FREE_MEM(sid);
  
}
void oglDeleteBuffs(GLuint vaoid[],GLuint vboid[]) {
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDeleteBuffers(3, vboid);
    
    glBindVertexArray(0);
    glDeleteVertexArrays(1, vaoid);
}
void oglPrintLog(GLuint sid) {
    GLint maxloglen,loglen;
    GLchar * log;
    glGetShaderiv(sid, GL_INFO_LOG_LENGTH, &maxloglen);
    log = ALLOC_MEM(GLchar,maxloglen);
    glGetShaderInfoLog(sid, maxloglen, &loglen, log);
    fprintf(stderr,"****************log:\n");
     fprintf(stderr,"%s\n",log);
     FREE_MEM(log);
    exit(1);
}
void oglSetUniformInt(GLuint pid,const char * name,int val) { 
    GLuint loc=glGetUniformLocation(pid,name);
    glUniform1i(loc,val);
}
/* float */
void oglMakeBuffsFlt(GLuint vaoid[],GLuint vboid[],int nv,int nt,float * vvec,float * nvec,int * tvec) {
    glGenBuffers(3,vboid);
    /* vertices */
    glBindBuffer(GL_ARRAY_BUFFER, vboid[0]);
    glBufferData(GL_ARRAY_BUFFER, (3 * nv) * sizeof(float), vvec, GL_STATIC_DRAW);
    /* normals */
    glBindBuffer(GL_ARRAY_BUFFER, vboid[1]);
    glBufferData(GL_ARRAY_BUFFER, (3 * nv) * sizeof(float), nvec, GL_STATIC_DRAW);
    /* triangles */
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboid[2]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, (3*nt) * sizeof(int), tvec, GL_STATIC_DRAW);

    // Create the VAO
    glGenVertexArrays( 1, vaoid );
    glBindVertexArray(vaoid[0]);

    glEnableVertexAttribArray(0);  // Vertex position
    glBindBuffer(GL_ARRAY_BUFFER, vboid[0]);
    glVertexAttribPointer( (GLuint)0, 3, GL_FLOAT, GL_FALSE, 0, 0 );

    glEnableVertexAttribArray(1);  // Vertex normal
    glBindBuffer(GL_ARRAY_BUFFER, vboid[1]);
    glVertexAttribPointer( (GLuint)1, 3, GL_FLOAT, GL_FALSE, 0, 0 );

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboid[2]);

    glBindVertexArray(0);
}

void oglSetUniformM4Flt(GLuint pid,const char * name,Mat4Flt val) { 
    GLuint loc=glGetUniformLocation(pid,name);
    glUniformMatrix4fv(loc, 1, GL_FALSE,val);
}


void oglSetUniformM3Flt(GLuint pid,const char * name,Mat3Flt val) { 
    GLuint loc=glGetUniformLocation(pid,name);
    glUniformMatrix3fv(loc, 1, GL_FALSE,val);
}

void oglSetUniformV3Flt(GLuint pid,const char * name,Flt * val) { 
    GLuint loc=glGetUniformLocation(pid,name);
    //glUniform3fv(loc,3,val);does not work
    glUniform3f(loc,val[0],val[1],val[2]);
}

void oglSetUniformV4Flt(GLuint pid,const char * name,Flt* val) { 
    GLuint loc=glGetUniformLocation(pid,name);
    glUniform4f(loc,val[0],val[1],val[2],val[3]);
}
void oglSetUniformFlt(GLuint pid,const char * name,Flt val) { 
    GLuint loc=glGetUniformLocation(pid,name);
    glUniform1f(loc,val);
}



  

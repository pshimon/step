/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "ogl.h"
#include "surf.h"
GLuint vaoid[1]; /* VAO id, no need in more than 1 */
GLuint vboid[3]; /* VBO ids vertices, normals, triangles*/
/* vertex shader source */
GLchar* vs_src = 
    "#version 330\n\
    layout(location=0) in vec3 v;\
    layout(location=1) in vec3 n;\
    out vec3 oclr;\
    uniform vec3 clr;\
    uniform vec3 lc;\
    uniform vec3 ld;\
    uniform mat4 mm;\
    uniform mat4 vm;\
    uniform mat4 pm;\
    uniform mat3 nm;\
    void main(void){\
	vec3 n1=normalize(nm*n);\
	vec3 l1=normalize(ld);\
	float sn=max( dot(l1,n1), 0.0 );\
	gl_Position = (pm * vm * mm) * vec4(v,1.0);\
	oclr = clr*lc*0.5f*(1.0f+sn);\
    }";
/* fragment shader source */
GLchar* fs_src =
    "#version 330\n\
    in vec3 oclr;\
    out vec4 fclr;\
    void main(void){\
    fclr = vec4(oclr,1.0);\
    }";

GLuint pid;/* program id */
t_m4  pm = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
} ;

t_m4  vm = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
} ;

t_m4  mm = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
} ;
t_m3  nm = {
  1, 0, 0,
  0, 1, 0,
  0, 0, 1} ;
int nv,nt;
t_v3 cl0={0.5f,0.5f,0.5f};/* clean color (gray)*/
t_v3 cl1={1.0f,0.0f,0.0f};/* surface color (red)*/
t_v3 cl2={1.0f,1.0f,1.0f};/* wire frame color (white)*/
t_v3 lc={1.0f,1.0f,1.0f}; /* light color */
t_v3 ld={0.0f,0.0f,1.0f};/* light direction */
t_f32 left,right,top,bot,near,far;
t_f32 ax=0.0f;
t_f32 ay=0.0f;
void keyboard_cb(unsigned char key,int x,int y) {
    switch (key) {
	case 27:
	case 'q':
	    exit(0);
	    break;
    }
}
void special_cb(int key,int x,int y) {
    t_f32 dx=0.1;
    t_f32 dy=0.1;
     switch (key) {
            case GLUT_KEY_LEFT:
                ax-=dx;
		if(ax<-M_PI) ax+=2.0*M_PI;
                break;
            case GLUT_KEY_RIGHT:
                ax+=dx;
		if(ax>M_PI) ax-=2.0*M_PI;
                break;
             case GLUT_KEY_UP:
                ay+=dy;
		if(ay>M_PI) ay-=2.0*M_PI;
                break;
            case GLUT_KEY_DOWN:
                ay-=dy;
		if(ay<-M_PI) ay+=2.0*M_PI;
                break;
      
       };
    glutPostRedisplay();

}

void idle_cb(void){
  glutPostRedisplay();
}

void clean_cb(void){
    ogl_delete_prog (pid);
    ogl_delete_buffs(vaoid,vboid);
}

void resize_cb(int w, int h) {
    t_f32 asp;
    glViewport(0, 0,w ,h );
    if(w<h) {
	asp=(h+0.0f)/(w+0.0f);
	ortho(pm,left,right,bot*asp,top*asp,near,far);
    } else {
	asp=(w+0.0f)/(h+0.0f);
	ortho(pm,left*asp,right*asp,bot,top,near,far);
    }

  glUseProgram(pid);
  ogl_set_uniform_m4(pid,"pm",pm);
  glUseProgram(0);

}

void display_cb(void){
    t_v3 v;
    t_m4 m1,m2;
    v[0]=1.0f;v[1]=0.0f;v[2]=0.0f;
    rotm4(m1, v,ax);
    v[0]=0.0f;v[1]=1.0f;v[2]=0.0f;
    rotm4(m2, v,ay);
    mxm4(mm,m1,m2);
    m3_from_m4(nm,mm);    
    glUseProgram(pid);
    ogl_set_uniform_m4(pid,"mm",mm);
    ogl_set_uniform_m4(pid,"vm",vm);
    ogl_set_uniform_m3(pid,"nm",nm);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
     glBindVertexArray(vaoid[0]);
    /* wire frame */
    ogl_set_uniform_v3(pid,"clr",cl2);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDrawElements(GL_TRIANGLES, 3*nt, GL_UNSIGNED_INT, (GLvoid*)0);
    /* surface */
    ogl_set_uniform_v3(pid,"clr",cl1);
    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1,1);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDrawElements(GL_TRIANGLES, 3*nt, GL_UNSIGNED_INT, (GLvoid*)0);
    /* done */
    glBindVertexArray(0);
    glUseProgram(0);
    glutSwapBuffers();
}
/*
int read_sphere(t_surf *s,char * fname) {
    FILE *fp=fopen(fname,"rb");
    int n,t,ret,j;
    size_t l;
    t_f32 x,y,z,r;
    if(fp==0) return -1;
    if(fread(&t,sizeof(int),1,fp)!=1) {ret=-5;goto abend;}
    if(fread(&n,sizeof(int),1,fp)!=1) {ret=-5;goto abend;}
    if(make_surf(s,t,n)) {ret=-2;goto abend;}
    l=3*t;
    if(fread(s->tvec,sizeof(int), l,fp)!= l) {ret=-3;goto abend;}
    l=3*n;
    if(fread(s->vvec,sizeof(t_f32), l,fp)!= l) {ret=-4;goto abend;}
    //if(fread(s->nvec,sizeof(t_f32), l,fp)!= l) {ret=-4;goto abend;}
    for(j=0;j<n;j++) {
	x=s->vvec[3*j+0];
	y=s->vvec[3*j+1];
	z=s->vvec[3*j+2];
	r=sqrt(x*x+y*y+z*z);
	s->nvec[3*j+0]=x/r;
	s->nvec[3*j+1]=y/r;
	s->nvec[3*j+2]=z/r;
    }
    ret=0;
abend:
    fclose(fp);    
    return ret;
}
*/    

int main(int argc, char* argv[]){
    int handle; 
    t_v3 v; 
    int ret,i;
    t_surf s;
    t_f32 rmax,fct=1.1,r;
    handle=glut_init_window(&argc,argv);
    if(handle<1) {
	fprintf(stderr,"failed to create window,handle =%d\n",handle);
	exit(1);
    }

    if(argc!=2) {
	fprintf(stderr,"usage: %s surf \n",argv[0]);
	exit(1);
    }
    init_surf(&s);
    ret=read_surf(&s,argv[1]);
    if(ret) {
	fprintf(stderr,"read_surf returns %d\n",ret);
	exit(1);
    }
    nv=s.nv;
    nt=s.nt;
    rmax=0.0;
    for(i=0;i<nv;i++) {
	r=norm3(s.vvec+3*i);
	if(r>rmax) rmax=r;
    }
    left=-rmax*fct;
    right=rmax*fct;
    bot=-rmax*fct;
    top=rmax*fct;
    
    near=0;
    far=2*rmax*fct;

  
 /*   if(!vs_src) vs_src=read_bytes("v2.gsl");
    if(!vs_src) {
	fprintf(stderr,"failed to read vs\n");
	exit(1);
    }
    if(!fs_src) fs_src=read_bytes("f2.gsl");
    if(!fs_src) {
	fprintf(stderr,"failed to read fs\n");
	exit(1);
    }
*/  glutSetWindowTitle(argv[1]);
    ogl_init_3d();
    pid=ogl_make_prog(vs_src,fs_src);
    ogl_make_buffs(vaoid,vboid,3*nv,3*nt,s.vvec,s.nvec,s.tvec);
    glClearColor(cl0[0],cl0[1],cl0[2], 1.0f);
    ogl_set_uniform_v3(pid,"ld",ld);
    ogl_set_uniform_v3(pid,"lc",lc);
    v[0]=0.0f;v[1]=0.0f;v[2]=-rmax*fct;
    translate4(vm,v);
    glutMainLoop();
    return 0;

}


/***********************************************************
* Shimon Panfil: Industrial Physics and Simulations        *
* http://industrialphys.com                                *
* THE SOFTWARE IS PROVIDED "AS IS",USE IT AT YOUR OWN RISK *
***********************************************************/
#include "ogl.h"
#include "ogl_glfw.h"
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
Mat4F  pm = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
} ;

Mat4F  vm = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
} ;

Mat4F  mm = {
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
} ;
Mat3F  nm = {
  1, 0, 0,
  0, 1, 0,
  0, 0, 1} ;
int nv,nt;
Vec3F cl0={0.5f,0.5f,0.5f};/* clean color (gray)*/
Vec3F cl1={1.0f,0.0f,0.0f};/* surface color (red)*/
Vec3F cl2={1.0f,1.0f,1.0f};/* wire frame color (white)*/
Vec3F lc={1.0f,1.0f,1.0f}; /* light color */
Vec3F ld={0.0f,0.0f,1.0f};/* light direction */
Flt left,right,top,bot,near,far;
Flt ax=0.0f;
Flt ay=0.0f;

void clean_cb(void){
    ogl_delete_prog (pid);
    ogl_delete_buffs(vaoid,vboid);
}

void resize_cb(int w, int h) {
    Flt asp;
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

static void render(){
    Vec3F v;
    Mat4F m1,m2;
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
    //glBindVertexArray(0);
    //glUseProgram(0);
 
 
}

int main(int argc, char* argv[]){
 
    Vec3F v; 
    int ret,i;
    TSurf s;
    Flt rmax,fct=1.1,r;
     int width=800, height=600;
  Flt dx=0.1;
    Flt dy=0.1;

    GLFWwindow *window;

    if(argc!=2) {
	fprintf(stderr,"usage: %s surf \n",argv[0]);
	exit(1);
    }
    window=glfw_init_window(width,height,argv[1]);
    if (!window) {
        fprintf(stderr, "GLFW3: failed to initialize\n");
        exit(EXIT_FAILURE);
    }
    printf ("(glfw) Renderer: %s\n", glGetString (GL_RENDERER));
    printf ("(glfw) OpenGL version supported %s\n",glGetString (GL_VERSION) );
    iniTSurf(&s);
    ret=readTSurf(&s,argv[1]);
    if(ret) {
	fprintf(stderr,"readTSurf returns %d\n",ret);
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
*/  
    ogl_init_3d();
    pid=ogl_make_prog(vs_src,fs_src);
    ogl_make_buffs(vaoid,vboid,3*nv,3*nt,s.vvec,s.nvec,s.tvec);
    glClearColor(cl0[0],cl0[1],cl0[2], 1.0f);
    ogl_set_uniform_v3(pid,"ld",ld);
    ogl_set_uniform_v3(pid,"lc",lc);
    v[0]=0.0f;v[1]=0.0f;v[2]=-rmax*fct;
    translate4(vm,v);
    
 
    while (!glfwWindowShouldClose(window)) {

	glfwGetFramebufferSize(window, &width, &height);
	resize_cb(width,height);
	render();	

	glfwPollEvents();

	if (GLFW_PRESS == glfwGetKey (window, GLFW_KEY_ESCAPE))
	    glfwSetWindowShouldClose(window, GL_TRUE);
	if (glfwGetKey (window, GLFW_KEY_LEFT)) {
	    ax-=dx;
	    if(ax<-M_PI) ax+=2.0*M_PI;
	}
	if (glfwGetKey (window, GLFW_KEY_RIGHT)) {
	    ax+=dx;
	    if(ax>M_PI) ax-=2.0*M_PI;
	}
	if (glfwGetKey (window, GLFW_KEY_UP)) {
	    ay+=dy;
	    if(ay>M_PI) ay-=2.0*M_PI;
	}
	if (glfwGetKey (window, GLFW_KEY_DOWN)) {
	    ay-=dy;
	    if(ay<-M_PI) ay+=2.0*M_PI;
	}
	glfwSwapBuffers(window);
    }
    clean_cb();
    glfwTerminate();
    return 0;

}



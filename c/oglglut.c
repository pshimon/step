/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#include "ogl.h" 
#include "oglglut.h"
int glutInitWindow(int * argc,char * argv[]) {
    GLenum GlewInitResult;
    int w=800,h=600,handle=0;
    char * wt="glut window";
    glutInit(argc, argv);
    glutInitContextVersion(3, 3); 
    glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
    glutInitContextProfile(GLUT_CORE_PROFILE);
    glutSetOption(
	GLUT_ACTION_ON_WINDOW_CLOSE,
	GLUT_ACTION_GLUTMAINLOOP_RETURNS
    );
    glutInitWindowSize(w,h);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    handle = glutCreateWindow(wt);
    if(handle <1) {
	fprintf(stderr,"glutCreateWindow returns %d\n",handle);
	return handle;
    }
    glutReshapeFunc(resize_cb);
    glutDisplayFunc(display_cb);
    glutIdleFunc(idle_cb);
    glutCloseFunc(clean_cb);
    glutKeyboardFunc(keyboard_cb);
    glutSpecialFunc(special_cb);
    //glutTimerFunc(0,timer_cb , 0);
    //extensions
    glewExperimental = GL_TRUE;
    GlewInitResult = glewInit();
    if (GLEW_OK != GlewInitResult) {
	fprintf(stderr,"ERROR: %s\n",glewGetErrorString(GlewInitResult));
	handle=0;
	return handle;
    }
 
    return handle;
}


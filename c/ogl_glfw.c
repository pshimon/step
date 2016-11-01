/***********************************************************************
* Shimon Panfil: Industrial Physics and Simulations                   **
* http://industrialphys.com                                           **
* THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           **
***********************************************************************/
#include "ogl.h"
#include "ogl_glfw.h"
GLFWwindow* glfw_init_window(int width,int height,char * title) {
    GLFWwindow* window;
     GLenum GlewInitResult;
     
    if (!glfwInit ()) {
	fprintf (stderr, "ERROR: could not start GLFW3\n");
	window=0;
	return window;
    }
    glfwWindowHint (GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint (GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint (GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint (GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    window = glfwCreateWindow (width,height,title, NULL, NULL);
    if (!window) {
	fprintf (stderr, "ERROR: could not open window with GLFW3\n");
	glfwTerminate();
	window=0;
	return window;
    }
    glfwMakeContextCurrent (window);
    glfwWindowHint (GLFW_SAMPLES, 4);
    glewExperimental = GL_TRUE;
    GlewInitResult = glewInit();
    if (GLEW_OK != GlewInitResult) {
	fprintf(stderr,"ERROR: %s\n",glewGetErrorString(GlewInitResult));
	window=0;
	return window;
    }
    return window;
}



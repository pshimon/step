########################################################################
# Shimon Panfil: Industrial Physics and Simulations                   ##
# http://industrialphys.com                                           ##
# THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           ##
########################################################################
ROOTDIR:=../step
include $(ROOTDIR)/makefile.def
include $(ROOTDIR)/makefile.gcc
include $(ROOTDIR)/makefile.rules
# fortran library sources
LSTF:=fdefs.f90 data_buf.f90 vec3d.f90 tsurf.f90  
# C library sources
LSTC:=databuf.c ogl.c ogl_glut.c  ogl_glfw.c read_bytes.c surf.c lplgf.c lplbem.c timers.c
# combining
SLST:=$(LSTF:.f90=.$(OEXT)) $(LSTC:.c=.$(OEXT))
#static library name
SLIB:=steplib.$(LEXT)
# fortran program sources
LSTBF:=data_buf_tst.f90 #mk_sphere.f90 cmp_tsurf.f90
# C program sources
LSTBC:=databuf-tst.c draw_surf1.c draw_surf2.c surf-test.c mksphere.c lplbem0-tst.c lplbem1-tst.c 
# combining
TRG:=$(LSTBF:.f90=.$(BEXT)) $(LSTBC:.c=.$(BEXT))
all:$(TRG)  
$(SLIB): $(SLST) 
	$(LBS)$(SLIB) $(SLST) 
$(TRG):$(SLIB)
clean:	
	$(RM) *.$(LEXT) *.$(BEXT) *.$(MEXT) *.$(OEXT) *.$(SHEXT)
#special treatment
GLLIBS:=  -lGL -lglut -lGLEW
GLLIBS2:=  -lGL -lglfw -lGLEW
draw_surf1.$(BEXT): draw_surf1.c 
	$(CC)  $(CLOPT) $(NOUT) $@ $^    $(GLLIBS) $(LM)
draw_surf2.$(BEXT): draw_surf2.c 
	$(CC)  $(CLOPT) $(NOUT) $@ $^    $(GLLIBS2) $(LM)
#lplbem1-tst.$(BEXT): lplbem1-tst.c $(SLIB)
#	$(CC) $(NOUT) $@ $^ $(MKLCOPT) $(MKLCLOPT) 




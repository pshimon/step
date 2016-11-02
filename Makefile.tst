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
LSTF:=fdefs.f90 data_buf.f90 #poly.f90 sdist1.f90 rng.f90 arrio.f90 sde_diff.f90 hsort.f90 stat.f90 
# C library sources
LSTC:=databuf.c ogl.c ogl_glut.c  ogl_glfw.c read_bytes.c surf.c  
# combining
SLST:=$(LSTF:.f90=.$(OEXT)) $(LSTC:.c=.$(OEXT))
#static library name
SLIB:=steplib.$(LEXT)
# fortran program sources
LSTBF:=data_buf_tst.f90
# C program sources
LSTBC:=databuf-tst.c draw_surf1.c draw_surf2.c surf-test.c
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






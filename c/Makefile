########################################################################
# Shimon Panfil: Industrial Physics and Simulations                   ##
# http://industrialphys.com                                           ##
# THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           ##
########################################################################
ROOTDIR:=..
include $(ROOTDIR)/makefile.def
include $(ROOTDIR)/makefile.gcc
include $(ROOTDIR)/makefile.rules
# C library sources
LSTC:= databuf.c # ogl.c  ogl_glut.c  read_bytes.c   surf.c 
SLST:=$(LSTC:.c=.$(OEXT))
#static library name
SLIB:=stepclib.$(LEXT)
# C program sources
LSTBC:=databuf-tst.c# draw_surf1.c surf-test.c
TRG:=$(LSTBC:.c=.$(BEXT))
all:$(TRG)  
$(SLIB): $(SLST) 
	$(LBS)$(SLIB) $(SLST) 
$(TRG):$(SLIB)
clean:	
	$(RM) *.$(LEXT) *.$(BEXT) *.$(MEXT) *.$(OEXT) *.$(SHEXT)
#special treatment
draw_surf1.$(BEXT): draw_surf1.c $(ULIB) $(GLIB)
	$(CC)  $(CLOPT) $(NOUT) $@ $^    $(GLLIBS) $(LM)





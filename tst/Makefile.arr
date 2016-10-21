########################################################################
# Shimon Panfil: Industrial Physics and Simulations                   ##
# http://industrialphys.com                                           ##
# THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           ##
########################################################################
ROOTDIR:=..
include $(ROOTDIR)/makefile.def
include $(ROOTDIR)/makefile.gcc
include $(ROOTDIR)/makefile.rules
LSTC:= arr.c  
SLST:=$(LSTF:.f90=.$(OEXT)) $(LSTC:.c=.$(OEXT))
SLIB:=alib.$(LEXT)
LSTBC:=arr-test.c
TRG:=$(LSTBF:.f90=.$(BEXT)) $(LSTBC:.c=.$(BEXT))
all:$(TRG)  
$(SLIB): $(SLST) 
	$(LBS)$(SLIB) $(SLST) 
$(TRG):$(SLIB)
clean:	
	$(RM) *.$(LEXT) *.$(BEXT) *.$(MEXT) *.$(OEXT) *.$(SHEXT)



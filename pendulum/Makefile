########################################################################
# Shimon Panfil: Industrial Physics and Simulations                   ##
# http://industrialphys.com                                           ##
# THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           ##
########################################################################
OEXT:=o
BEXT:=x
SHEXT:=so
FPIC:=-fpic
RM:=rm -f
LEXT:=a
LBS:=ar crv 
LM:=-lm -lrt
NOUT:= -o 
MEXT:=mod
vpath %.f90 . ../utils
FC:=gfortran
FOPT:=-O3 $(FPIC) -Wall -g  -fopenmp
FPE:= -ffpe-trap=invalid,zero,overflow
FLOPT:=$(FOPT) $(FPE)

%.$(OEXT): %.f90
	$(FC) -c $(FOPT)    $<  $(NOUT) $@
%.$(BEXT): %.f90
	$(FC) $(FLOPT)  $(NOUT) $@ $^ 
	
LSTF:=fdefs.f90  arrio.f90 kpend.f90 
SLST:=$(LSTF:.f90=.$(OEXT)) 
SLIB:=pendlib.$(LEXT)
LSTBF:=xkpend.f90
LSTBC:=data_int.c data_cmb.c
TRG:=$(LSTBF:.f90=.$(BEXT)) 
all:$(TRG)  
$(SLIB): $(SLST) 
	$(LBS)$(SLIB) $(SLST) 
$(TRG):$(SLIB)

clean:	
	$(RM) *.$(LEXT) *.$(BEXT) *.$(MEXT) *.$(OEXT) *.$(SHEXT)



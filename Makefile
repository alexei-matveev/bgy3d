# Include Slepc and Petsc
# on Debian Lenny do
#
#   export PETSC_DIR=/usr/lib/petscdir/2.3.3
#
# Then
#
#   make all
#
include ${PETSC_DIR}/bmake/common/base

# Shell
SHELL = /bin/sh

#USERFLAGS
USERFLAGS += -DCENTRALDIFF
#USERFLAGS += -DMATPRECOND
USERFLAGS += -DINTRA2
USERFLAGS += -DL_BOUNDARY
#USERFLAGS += -DH2O
#USERFLAGS += -DCS2
#USERFLAGS += -DPROPANE
#USERFLAGS += -DN2
USERFLAGS += -DHCl


# Path to sources
srcdir = .

# Compiler and compiler options
CC       = gcc 
CPPFLAGS = 
CFLAGS   = -O2 -DFFT_FFTW $(USERFLAGS) 
#CFLAGS	 =   -g  -Wall  -DFFT_FFTW $(USERFLAGS)
LDFLAGS  = 



INCDIRS = ${PETSC_INCLUDE} -I./fft/ -I/home/opt64/packages/fftw-2.1.5/include/ 
LIBDIRS = 
LIBS    = -lm -lfftw -lfftw_mpi  ${PETSC_LIB} 

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

OBJECTS =  bgy3d.o functions.o hnc3d.o bgy3ddiv.o bgy3dtest.o bgy3dfourier.o bgy3dmolecule.o bgy3dH2O.o bgy3dH2OS.o bgy3dH2ONewton.o bgy3dH2OSNewton.o bgy3dH2O_solutes.o fft/fft_3d.o fft/fft_3d_f.o fft/pack_3d.o fft/remap_3d.o fft/factor.o 

all: bgy3d

bgy3d : $(OBJECTS)
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $(OBJECTS) ${LIBDIRS} ${LIBS}
#	cp $@ test/
#	cp $@ test_H2O/
#	cp $@ test_H2OII/
#	cp $@ test_H2OIII/

myclean:
	rm -f *.o fft/*.o 
	rm -f bgy3d
.PHONY: myclean

distclean:
	rm -f *.o *.d  fft/*.o fft/*.d 
	rm -f bgy3d
.PHONY: distclean

include $(OBJECTS:.o=.d)

# This do-it-all rule is specific to gcc version >= 3.0,
# and updates the object file and the dependency file at the same time.
# The -MP also causes it to output empty targets for each header.
# Comment out this rule if you are using a different compiler.
#%.o %.d: %.c
#	$(CC) -MT "$*.o $*.d " -MD -MP $(CFLAGS) $(LDFLAGS) $(INCDIRS) -c $<
%.o %.d: %.c
	$(CC) -M -MF '$*.d' -MP $(CFLAGS) $(LDFLAGS) $(INCDIRS) -c $<
	$(CC) -o '$*.o' $(CFLAGS) $(LDFLAGS) $(INCDIRS) -c $<

# The next two rules are fairly portable across compilers.
%.o: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCDIRS) -c $<

%.d: %.c
	set -e; $(CC) -M $(CFLAGS) $(LDFLAGS) $(INCDIRS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
	[ -s $@ ] || rm -f $@




#---End of Makefile---

# Include Slepc and Petsc
# on Debian Lenny do
#
#   export PETSC_DIR=/usr/lib/petscdir/2.3.3
#
# Then
#
#   make -s
#
all: bgy3d

include ${PETSC_DIR}/bmake/common/base
# include ${PETSC_DIR}/conf/base # on Ubuntu 12.04 LTS

# Shell
SHELL = /bin/sh

#
# Compile rarely used solvers.   These solvers have not beed converted
# to direct  use of  FFTW MPI API  and rely  on the wrappers  by Steve
# Plimpton,  see  ./fft  directory.    The  fft_3d.h  header  needs  a
# -DFFT_FFTW  in  CFLAGS,  otherwise   it  cannot  be  included.  Keep
# WITH_EXTRA_SOLVERS = 1  by default to minimize bit  rotting of those
# solvers:
#
WITH_EXTRA_SOLVERS = 1

USERFLAGS = -DFFT_FFTW
ifeq ($(WITH_EXTRA_SOLVERS),1)
USERFLAGS += -DWITH_EXTRA_SOLVERS
endif

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
CFLAGS   = -std=c99 -Wall -Wextra -O2 $(USERFLAGS)
LDFLAGS  =



INCDIRS = ${PETSC_INCLUDE} -I./fft
LIBDIRS = 
LIBS    = -lm -lfftw -lfftw_mpi  ${PETSC_LIB} 

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

bgy3d-objs = \
	bgy3d-main.o \
	bgy3d.o \
	bgy3dH2O.o \
	bgy3dH2OS.o \
	bgy3d-solutes.o \
	bgy3d-getopt.o \
	bgy3d-fft.o \

bgy3d-extra-objs = \
	bgy3d-simple.o \
	hnc3d.o \
	bgy3d-molecule.o \
	bgy3ddiv.o \
	bgy3dfourier.o \
	bgy3d-test.o \
	bgy3dH2ONewton.o \
	bgy3dH2OSNewton.o \
	bgy3d_MG.o \

fft3d-objs = \
	fft/fft_3d.o \
	fft/fft_3d_f.o \
	fft/pack_3d.o \
	fft/remap_3d.o \
	fft/factor.o

OBJECTS = $(bgy3d-objs)

ifeq ($(WITH_EXTRA_SOLVERS),1)
OBJECTS += $(bgy3d-extra-objs) $(fft3d-objs)
endif

bgy3d : $(OBJECTS)
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $(OBJECTS) ${LIBDIRS} ${LIBS}

#
# Dont call the target "test" because we have a directory called so:
#
test-all:
	$(MAKE) -C ./test

#
# One of the include files defines a target named clean already, we
# augment its prerequisites here:
#
clean: myclean

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
%.d: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCDIRS) -M -MF $(*).d -MP $(<)

# The next two rules are fairly portable across compilers.
%.o: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCDIRS) -o $(*).o -c $(<)

# node capital D here, this rule has no effect:
%.D: %.c
	set -e; $(CC) -M $(CFLAGS) $(LDFLAGS) $(INCDIRS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
	[ -s $@ ] || rm -f $@

TAGS: $(OBJECTS:.o=.c)
	etags $(OBJECTS:.o=.c)

#---End of Makefile---

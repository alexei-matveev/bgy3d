# Include Slepc and Petsc
# on Debian Lenny do
#
#   export PETSC_DIR=/usr/lib/petscdir/2.3.3
#
# On Ubuntu 12.04 LTS do
#
#   export PETSC_DIR=/usr/lib/petscdir/3.1
#
# Then
#
#   make -s
#
all: bgy3d

include $(PETSC_DIR)/bmake/common/base
# include $(PETSC_DIR)/conf/base # on Ubuntu 12.04 LTS

# Shell
SHELL = /bin/sh

#
# Guile  is  a  Scheme interpreter.  Set  this  to  one to  compile  a
# BGY3d-enabled interpreter:
#
WITH_GUILE = 0

#
# Compile rarely used solvers.   These solvers have not beed converted
# to direct  use of  FFTW MPI API  and rely  on the wrappers  by Steve
# Plimpton,  see  ./fft  directory.    The  fft_3d.h  header  needs  a
# -DFFT_FFTW  in  CFLAGS,  otherwise   it  cannot  be  included.  Keep
# WITH_EXTRA_SOLVERS = 1  by default to minimize bit  rotting of those
# solvers:
#
WITH_EXTRA_SOLVERS = 0

#
# Compile  a   shared  library  libbgy3d.so  and  link   that  to  the
# executable.   May require  setting LD_LIBRARY_PATH  in order  to run
# that. The following flag is used in expressions like $(if $(shared),
# ...) so  in order to turn  the feature off  it is not enough  to set
# shared = 0, one has to comment the line instead:
#
# shared = 1

USERFLAGS = -DFFT_FFTW
ifeq ($(WITH_EXTRA_SOLVERS),1)
USERFLAGS += -DWITH_EXTRA_SOLVERS
endif

USERFLAGS += -DCENTRALDIFF
#USERFLAGS += -DMATPRECOND
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
CFLAGS   = -std=c99 -Wall -Wextra -O3 $(USERFLAGS) $(if $(shared), -fPIC)
LDFLAGS  =



INCDIRS = $(PETSC_INCLUDE) -I./fft
fftw3-libs = -lfftw3_mpi -lfftw3
fftw2-libs = -lfftw_mpi -lfftw
rfftw2-libs = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
LIBS = $(fftw3-libs) -lm $(PETSC_LIB)

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

libbgy3d.a = \
	hnc3d.o \
	hnc3d-sles.o \
	rism-dst.o \
	bgy3d.o \
	bgy3d-force.o \
	bgy3d-pure.o \
	bgy3d-impure.o \
	bgy3d-solvents.o \
	bgy3d-solutes.o \
	bgy3d-getopt.o \
	bgy3d-poisson.o \
	bgy3d-dirichlet.o \
	bgy3d-snes.o \
	bgy3d-vec.o \
	bgy3d-fft.o \
	bgy3d-fftw3.o \
	bgy3d-potential.o

ifeq ($(WITH_GUILE),1)
	libbgy3d.a += bgy3d-guile.o
	LIBS += $(shell guile-config link)
	INCDIRS += $(shell guile-config compile)
	USERFLAGS += -DWITH_GUILE
endif

bgy3d-extra-objs = \
	bgy3d-simple.o \
	hnc3d-newton.o \
	bgy3d-molecule.o \
	bgy3ddiv.o \
	bgy3d-fourier.o \
	bgy3d-test.o \
	bgy3d-newton-pure.o \
	bgy3d-newton-impure.o \
	bgy3d-multigrid.o \

fft3d-objs = \
	fft/fft_3d.o \
	fft/fft_3d_f.o \
	fft/pack_3d.o \
	fft/remap_3d.o \
	fft/factor.o

ifeq ($(WITH_EXTRA_SOLVERS),1)
libbgy3d.a += $(bgy3d-extra-objs) $(fft3d-objs)
endif

OBJECTS = bgy3d-main.o

libbgy3d.a: $(libbgy3d.a)
	$(AR) ruv $@  $(^)
	$(RANLIB) $@

libbgy3d.so: $(libbgy3d.a)
	$(CC) -shared -Wl,-soname,libbgy3d.so.1 -o libbgy3d.so.1.0 $(libbgy3d.a)
	ln -sf libbgy3d.so.1.0 libbgy3d.so.1
	ln -sf libbgy3d.so.1 libbgy3d.so

bgy3d: $(OBJECTS) $(if $(shared), libbgy3d.so, libbgy3d.a)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) -L. -lbgy3d $(LIBS)

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
	rm -f *.a *.so *.o fft/*.o *.bin *.info
	rm -f bgy3d
.PHONY: myclean

distclean:
	rm -f *.o *.so *.d  fft/*.o fft/*.d
	rm -f bgy3d
.PHONY: distclean

include $(OBJECTS:.o=.d)
include $(libbgy3d.a:.o=.d)

# This do-it-all rule  is specific to gcc version  >= 3.0, and updates
# the object file  and the dependency file at the  same time.  The -MP
# also causes it to output empty targets for each header.  Comment out
# this rule if you are using a different compiler.
# %.o %.d: %.c
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

TAGS: $(OBJECTS:.o=.c) $(libbgy3d.a:.o=.c)
	etags $(^)

tags: $(OBJECTS:.o=.c) $(libbgy3d.a:.o=.c)
	ctags $(^)

#---End of Makefile---

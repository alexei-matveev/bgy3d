#
# Include Slepc and  Petsc.  Petsc relies on a  lot of libraries.  The
# easiest way to set compiler flags is to use the variables defined by
# Petsc  makefiles. For  Petsc  2.3 the  file  included below  defines
# PETSC_CC_INCLUDES  and  PETSC_LIB  among  other things.   On  Debian
# derivatives do
#
#   export PETSC_DIR=/usr/lib/petsc
#
# This should  be a  symbolic link to  the actual location  managed by
# update-alternatives.   The  structure  of config  directory  changed
# around Petsc 3.0 --- adapt the path accordingly:
#
include $(PETSC_DIR)/bmake/common/variables # on Lenny
# include $(PETSC_DIR)/conf/variables # on Wheezy, Ubuntu 12.04

#
# Finally, execute
#
#   make -s
#
all: bgy3d

# Shell
SHELL = /bin/sh

#
# Guile  is  a  Scheme interpreter.  Set  this  to  one to  compile  a
# BGY3d-enabled interpreter:
#
WITH_GUILE = 0

#
# Some  code  is  written  in  Fortran using  F2008  features.   Older
# compilers   including   GFortran   4.3   on  Lenny   cannot   handle
# that. Disable those features:
#
WITH_FORTRAN = 0

#
# Compile  a   shared  library  libbgy3d.so  and  link   that  to  the
# executable.   May require  setting LD_LIBRARY_PATH  in order  to run
# that. The following flag is used in expressions like $(if $(shared),
# ...) so  in order to turn  the feature off  it is not enough  to set
# shared = 0, one has to comment the line instead:
#
# shared = 1

USERFLAGS = -DFFT_FFTW
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
FC       = gfortran
CFLAGS   = -std=c99 -Wall -Wextra -O3 $(USERFLAGS) $(if $(shared), -fPIC)
FFLAGS   = -std=f2008 -Wall -O3 $(if $(shared), -fPIC)
LDFLAGS  =



INCDIRS = $(PETSC_CC_INCLUDES) -I./fft
fftw3-libs = -lfftw3_mpi -lfftw3
fftw2-libs = -lfftw_mpi -lfftw
rfftw2-libs = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
LIBS = $(fftw3-libs) -lm $(PETSC_LIB)

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

libbgy3d.a = $(c-objs) $(f-objs)

c-objs = \
	hnc3d.o \
	hnc3d-sles.o \
	rism-dst.o \
	rism-rdf.o \
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
	c-objs += bgy3d-guile.o
	LIBS += $(shell guile-config link)
	INCDIRS += $(shell guile-config compile)
	USERFLAGS += -DWITH_GUILE
endif

ifeq ($(WITH_FORTRAN),1)
	f-objs += rism.o lebed/lebed.o lebed/Lebedev-Laikov.o
	USERFLAGS += -DWITH_FORTRAN
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

clean:
	rm -f *.a *.so *.o *.bin *.info
	rm -f bgy3d

distclean:
	rm -f *.o *.so *.d
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

#
# FIXME:  this effectively  ignores the  dependencies  between fortran
# sources:
#
%.d: %.f90
	touch $(@)

%.d: %.F
	touch $(@)

# The next two rules are fairly portable across compilers.
%.o: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCDIRS) -o $(*).o -c $(<)

%.o: %.f90
	$(FC) $(FFLAGS) -o $(*).o -c $(<)

# node capital D here, this rule has no effect:
%.D: %.c
	set -e; $(CC) -M $(CFLAGS) $(LDFLAGS) $(INCDIRS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
	[ -s $@ ] || rm -f $@

TAGS: $(OBJECTS:.o=.c) $(c-objs:.o=.c) $(f-objs:.o=.f90)
	etags $(^)

tags: $(OBJECTS:.o=.c) $(c-objs:.o=.c) $(f-objs:.o=.f90)
	ctags $(^)

#---End of Makefile---

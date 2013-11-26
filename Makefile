#
# Petsc relies on a lot of libraries.  The easiest way to set compiler
# flags is to  use the variables defined by  Petsc makefiles. The file
# "variables" included  below defines PETSC_CC_INCLUDES  and PETSC_LIB
# among other things. To overwrite the default location do e.g.:
#
#   export PETSC_DIR=/usr/lib/petsc
#
# On Debian derivatives this is a symbolic link to the actual location
# managed by  update-alternatives.  The structure  of config directory
# was  different   before  Petsc  3.0  ---  adapt   the  include  path
# accordingly:
#
PETSC_DIR ?= /usr/lib/petsc
include $(PETSC_DIR)/conf/variables # on Wheezy, Ubuntu 12.04
# include $(PETSC_DIR)/bmake/common/variables # on Lenny

#
# Finally, execute
#
#   make clean
#   make -s
#
# FIXME: Note that generating  dependencies for Fortran files may fail
# in a  parallel build  (make -j).   Delete all the  *.d files  if you
# encounter  the  problem  with  GFortran  being  unable  to  generate
# dependencies.  A  serial "make clean" might  be sometimes sufficient
# to initally populate *.d files.
#
all: bgy3d

# Shell
SHELL = /bin/sh

#
# Guile  is  a  Scheme interpreter.  Set  this  to  one to  compile  a
# BGY3d-enabled interpreter:
#
WITH_GUILE = 1

#
# Some  code  is  written  in  Fortran using  F2008  features.   Older
# compilers   including   GFortran   4.3   on  Lenny   cannot   handle
# that. Disable those features:
#
WITH_FORTRAN = 1

#
# Compile  a   shared  library  libbgy3d.so  and  link   that  to  the
# executable.   May require  setting LD_LIBRARY_PATH  in order  to run
# that. The following flag is used in expressions like $(if $(shared),
# ...) so  in order to turn  the feature off  it is not enough  to set
# shared = 0, one has to comment the line instead:
#
# shared = 1

USR-FLAGS = -DFFT_FFTW
USR-FLAGS += -DL_BOUNDARY
# USR-FLAGS += -DH2O
# USR-FLAGS += -DCS2
# USR-FLAGS += -DPROPANE
# USR-FLAGS += -DN2
USR-FLAGS += -DHCl

# Path to sources
srcdir = .

# Compiler and compiler options
CC       = gcc
FC       = gfortran
CFLAGS = -g -std=c99 -Wall -Wextra -Ofast $(PIC-FLAGS) $(USR-FLAGS)
FFLAGS = -g -std=f2008 -Wall -O3 $(PIC-FLAGS) $(OMP-FLAGS) $(DBG-FFLAGS)
LDFLAGS  = $(OMP-FLAGS)

# Flags to generate position independent code (PIC) for shared objects
# happen to be the same for gcc and gfortran:
PIC-FLAGS = $(if $(shared), -fPIC)

# I was not able to get a real speedup for urany/water with N=8192, so
# this  remains disabled.   Note  that -fopenmp  is  accepted by  both
# gcc/gfortran (not used in C though) and also needs to be supplied at
# link stage:
OMP-FLAGS = # -fopenmp

# Fortran flags to assist debugging and experiments:
DBG-FFLAGS = # -fcheck-array-temporaries -fbounds-check -fexternal-blas -fblas-matmul-limit=1



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
	bgy3d-mat.o \
	bgy3d-fft.o \
	bgy3d-fftw3.o \
	bgy3d-potential.o

ifeq ($(WITH_GUILE),1)
	c-objs += bgy3d-guile.o
	LIBS += $(shell guile-config link)
	INCDIRS += $(shell guile-config compile)
	USR-FLAGS += -DWITH_GUILE
endif


#
# There is a chicken and  egg problem here.  To make GFortran generate
# dependencies for Fortran files  the compiler expects all *.mod files
# the source depends on to be present. The order of the files here was
# manually crafted so that  a *serial* make will generate dependencies
# of Fortran modules in  the proper sequence simultaneousely producing
# the *.mod  files. FIXME:  sill fails for  parallel build as  in make
# -j4.
#
ifeq ($(WITH_FORTRAN),1)
	f-objs += \
		rism.o \
		drism.o \
		snes.o \
		bessel.o \
		fft.o \
		linalg.o \
		options.o \
		lisp.o \
		foreign.o \
		units.o \
		kinds.o \
		lebed/lebed.o \
		lebed/Lebedev-Laikov.o
	USR-FLAGS += -DWITH_FORTRAN
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
	rm -f *.a *.so *.o *.mod *.bin *.info
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
# FIXME: running the command also  produces the *.mod file, but no *.o
# file.  Also GFortran expects the *.mod files of the prerequisites to
# exist when parsing the source:
#
%.d: %.f90
	$(FC) -cpp -M -MF $(@) $(<)

%.d: %.F
	touch $(@)

# The next two rules are fairly portable across compilers.
%.o: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCDIRS) -o $(*).o -c $(<)

%.o %.mod: %.f90
	$(FC) $(FFLAGS) -o $(*).o -c $(<)

# Note capital D here, this rule has no effect:
%.D: %.c
	set -e; $(CC) -M $(CFLAGS) $(LDFLAGS) $(INCDIRS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
	[ -s $@ ] || rm -f $@

TAGS: $(OBJECTS:.o=.c) $(c-objs:.o=.c) $(f-objs:.o=.f90)
	etags $(^)

tags: $(OBJECTS:.o=.c) $(c-objs:.o=.c) $(f-objs:.o=.f90)
	ctags $(^)

#---End of Makefile---

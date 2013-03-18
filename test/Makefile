# Shell
SHELL = /bin/sh

#
# Flags  for  bgy3d executable.  With  current  version  it takes  312
# iterations to converge butanoic acid in HCl.
#
BGY-FLAGS = \
	--bgy \
	--N 32 \
	--rho 0.018 \
	--beta 1.1989 \
	--norm-tol 1.0e-7 \
	--max-iter 1000 \
	--L 10.0 \
	--damp-start 1.0 \
	--lambda 0.02 \
	--snes-solver jager \

HNC-FLAGS = \
	--hnc \
	--N 128 \
	--L 10.0 \
	--beta 0.606060 \
	--max-iter 320 \
	--lambda 0.1 \
	--rho 0.1 \
	--norm-tol 1e-12 \

#
# In case $(EXE) is not defined, assume the executable is one level up
# in directory structure:
#
EXE ?= ../bgy3d

solutes = \
	hydrogen_chloride \
	carbon_disulfide \
	water \
	methanol \
	butanoic_acid \
	hexane

#
# These  are  phony  targets,  they   are  not  created  as  files  to
# re-evaluate the diffs every time Make is executed:
#
bgy-diffs = \
	hydrogen_chloride.g2.diff \
	$(patsubst %,%.g1.diff,$(solutes))

hnc-diffs = \
	LJ.hnc3d.g2.diff \
	LJ.hnc3d.g1.diff \
	$(patsubst %,%.hnc3d.g1.diff,$(solutes))

#
# By default "make"  all the diffs. Unomment the  second entry to also
# run HNC3D tests:
#
all: bgy-tests # hnc-tests
bgy-tests: $(bgy-diffs)
hnc-tests: $(hnc-diffs)

#
# Target specific variables. Four cases (g1, g2) x (BGY, HNC) at most:
#
$(bgy-diffs): base-flags = $(BGY-FLAGS)
$(hnc-diffs): base-flags = $(HNC-FLAGS)

#
# For g1 calculations need to specify additional flags:
#
%.hnc3d.g1 %.g1: solute-flags= --solute "$(subst _, ,$(*))"


#
# Prerequisites for a  diff are the two files, keep  them in the order
# old, new:
#
%.diff: out/% %
	diff $(^)

#
# Dont  delete  distribution  summary  files rebuilt  as  intermediate
# targets.  Running a  test takes some time. FIXME: why  do I need the
# second rule, given the first?
#
.PRECIOUS: %.g2 %.g1
.PRECIOUS: %.hnc3d.g2 %.hnc3d.g1

#
# Nothing to  be done for  reference outputs.  Without this  rule Make
# will try to  re-build the reference outputs too  when the executable
# gets changed:
#
out/%:
	 # nothing to be done, is there a better way to say that?

#
# Command line  to print  the distribution summary  from *.bin  or *.m
# file(s). Darcs does not track file attributes, so that the *.py file
# will not always have an execution bit set in a fresh repository:
#
moments = python ../python/moments.py

#
# This script is also used in QM regression tests, see ../test-qm:
#
make-summary = $(SHELL) ./summary.sh

#
# To produce a *.g2 or *.g1 summary file execute BGY3d and examine the
# resulting  distributions.  Updating  the executable  should  lead to
# re-running  the  test, so  add  it  as  a pre-requisite.   Also  the
# solvent-solvent g2-distribution  is required for  g1-distribution of
# the solvent in the presence of a solute:
#
# Testing HNC3d code is not as elaborate, thus the special case.
#
%.hnc3d.g2: ../bgy3d
	$(EXE) $(base-flags) | tee $(@).out
	$(moments) g00.bin > $(@)
	$(make-summary) $(@).out >> $(@)

%.hnc3d.g1: ../bgy3d LJ.hnc3d.g2
	$(EXE) $(base-flags) $(solute-flags) | tee $(@).out
	$(moments) g0.bin > $(@)
	$(make-summary) $(@).out >> $(@)

%.g1: ../bgy3d hydrogen_chloride.g2
	$(EXE) $(base-flags) $(solute-flags) | tee $(@).out
	$(moments) g0.bin g1.bin > $(@)
	$(make-summary) $(@).out >> $(@)

%.g2: ../bgy3d
	$(EXE) $(base-flags) | tee $(@).out
	$(moments) g00.bin g11.bin g01.bin > $(@)
	$(make-summary) $(@).out >> $(@)


# This only works if ../bgy3d executable was built WITH_GUILE:
%.m: %.bin
	../guile/runbgy.scm dump $(<) > $(@)

# This pops a graphics window for  a second. Is there a way to disable
# that?
%.m.png: %.m
	../octave/bgy3d_plot.m $(<)

# Do  not delte  intermediate  ascii *.m  files  produced when  making
# figures:
.PRECIOUS: %.m

clean:
	rm -f *.g1 *.g2 *.bin *.info *.m *.out

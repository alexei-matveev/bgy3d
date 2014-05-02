# Shell
SHELL = /bin/sh

#
# Flags for  bgy3d executable.  With  current settings it takes  up to
# 1000 iterations some cases.
#
BGY-FLAGS = \
	--bgy \
	--N $(N) \
	--L $(L) \
	--rho 0.018 \
	--beta 1.1989 \
	--norm-tol 1.0e-7 \
	--max-iter 1000 \
	--lambda 0.02 \
	--solvent "hydrogen chloride" \
	--snes-solver jager \

solvent = water
HNC-FLAGS = \
	--hnc \
	--N $(N) \
	--L $(L) \
	--beta 1.6889 \
	--max-iter 1000 \
	--rho 0.033427745 \
	--norm-tol 1e-12 \
	--solvent "$(solvent)" \
	--closure KH

#
# In  case $(cmd) is  not defined,  use the  script that  emulates its
# behaviour.  Here $(exe) is the actual executable, $(cmd) is a script
# with  that executable  in  shebang line  accepting and  interpreting
# command line  flags to tune  settings. Finally, $(run) is  a command
# line prefix that can be used to interprete arbitrary Scheme files:
#
TOP ?= ..
cmd ?= mpirun $(TOP)/guile/runbgy.scm
exe ?= $(TOP)/bgy3d
run = $(exe) -L .. -s

solutes = \
	water \
	hydrogen_chloride \
	carbon_disulfide \
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
	$(solvent).hnc3d.g2.diff \
	$(patsubst %,%.hnc3d.g1.diff,$(solutes))

#
# By default "make" all the  diffs. FIXME: An artificial dependency is
# there  because  both  type   of  tests  generate  g2-files  used  by
# g1-calculations. They cannot be run in parallel.
#
all: rism-tests bgy-tests hnc-tests
bgy-tests: $(bgy-diffs)
hnc-tests: $(hnc-diffs) bgy-tests
rism-tests: run-ions-out.diff

run-ions-out: run-ions.scm $(exe)
	$(run) $(<) --snes-solver trial

#
# Target specific variables. Four cases (g1, g2) x (BGY, HNC) at most:
#
$(bgy-diffs): N = 32
$(hnc-diffs): N = 96
$(bgy-diffs): L = 10.0
$(hnc-diffs): L = 10.0
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
# file(s). The  important flags are  those that specify the  domain (L
# and N). In fact L could  be arbitrary --- then the units the moments
# are printed in are arbitrary too:
#
# moments = python $(TOP)/python/moments.py
moments = $(cmd) moments --L $(L) --N $(N)

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
%.hnc3d.g2: $(exe)
	$(cmd) $(base-flags) | tee $(@).out
	$(moments) g00.bin g01.bin g02.bin g11.bin g12.bin g22.bin > $(@)
	$(make-summary) $(@).out >> $(@)

%.hnc3d.g1: $(exe) $(solvent).hnc3d.g2
	$(cmd) $(base-flags) $(solute-flags) | tee $(@).out
	$(moments) g0.bin g1.bin g2.bin > $(@)
	$(make-summary) $(@).out >> $(@)

%.g1: $(exe) hydrogen_chloride.g2
	$(cmd) $(base-flags) $(solute-flags) | tee $(@).out
	$(moments) g0.bin g1.bin > $(@)
	$(make-summary) $(@).out >> $(@)

%.g2: $(exe)
	$(cmd) $(base-flags) | tee $(@).out
	$(moments) g00.bin g11.bin g01.bin > $(@)
	$(make-summary) $(@).out >> $(@)


# This only works if ../bgy3d executable was built WITH_GUILE:
%.m: %.bin
	$(TOP)/guile/runbgy.scm dump $(<) > $(@)

# This may pop  a graphics window for a second,  unless you use figure
# ("visible", "off").
%.m.png: %.m
	$(TOP)/octave/bgy3d_plot.m $(<)

# Do  not delte  intermediate  ascii *.m  files  produced when  making
# figures:
.PRECIOUS: %.m

all-bin=$(wildcard *.bin)
all-png: $(all-bin:.bin=.m.png)

clean:
	rm -f *.g1 *.g2 *.bin *.info *.m *.out

# Shell
SHELL = /bin/sh

#
# Flags  for  bgy3d executable.  With  current  version  it takes  312
# iterations to converge butanoic acid in HCl.
#
CMDFLAG = \
	--N 32 \
	--rho 0.018 \
	--beta 1.1989 \
	--norm-tol 1.0e-2 \
	--max-iter 320 \
	--L 10.0 \
	--damp-start 1.0 \
	--lambda 0.02

#
# In case $(EXE) is not defined, assume the executable is one level up
# in directory structure:
#
EXE ?= ../bgy3d

#
# These  are  phony  targets,  they   are  not  created  as  files  to
# re-evaluate the diffs every time Make is executed:
#
diffs = \
	hydrogen_chloride.g2.diff \
	hydrogen_chloride.g1.diff \
	carbon_disulfide.g1.diff \
	water.g1.diff \
	methanol.g1.diff \
	butanoic_acid.g1.diff \
	hexane.g1.diff

#
# By default "make" all the diffs:
#
all: $(diffs)

#
# Prerequisites for a  diff are the two files, keep  them in the order
# old, new:
#
%.diff: out/% %
	diff $(^)

#
# Dont  delete  distribution  summary  files rebuilt  as  intermediate
# targets.  Running a test takes some time:
#
.PRECIOUS: %.g2 %.g1

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
# To produce a *.g2 or *.g1 summary file execute BGY3d and examine the
# resulting  distributions.  Updating  the executable  should  lead to
# re-running  the  test, so  add  it  as  a pre-requisite.   Also  the
# solvent-solvent g2-distribution  is required for  g1-distribution of
# the solvent in the presence of a solute:
#
# Testing HNC3d code  is not as elaborate, this is  a special case. At
# the  moment  it  requires  a  "molecule"  file  to  be  loaded  with
# Load_Molecule():
#
hnc3d.g2: ../bgy3d
	../bgy3d --HNC --L 5.0 --N 128 --beta 0.606060 --max-iter 320 --lambda 0.1 --rho 0.1 --norm-tol 1e-12
	$(moments) g00.bin > $(@)

%.g1: ../bgy3d hydrogen_chloride.g2
	$(EXE) --BGYM2Site $(CMDFLAG) --solute "$(subst _, ,$(*))"
	$(moments) g0.bin g1.bin > $(@)

%.g2: ../bgy3d
	$(EXE) --BGY2Site $(CMDFLAG)
	$(moments) g00.bin g11.bin g01.bin > $(@)

clean:
	rm -f *.g1 *.g2 *.bin *.info *.m

# Shell
SHELL = /bin/sh

#
# In case $(EXE) is not defined, assume the executable here. Note that
# the stock  PG script  will not interprete  the BGY3d section  in the
# *.scm input files (yet):
#
EXE ?= ./runqm

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
diffs = \
	$(patsubst %,%.g1.diff,$(solutes))

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
# This script is also used in MM regression tests, see ../test:
#
make-summary = $(SHELL) ../test/summary.sh

#
# To produce a *.g2 or *.g1 summary file execute BGY3d and examine the
# resulting distributions:
#
%.g1: %.scm
	$(EXE) $(*).scm | tee $(*).out
	$(moments) g0.bin g1.bin > $(@)
	$(SHELL) summary.sh $(*).out >> $(@)
	rm -f *.info

clean:
	rm -f *.g1 *.g2 *.bin *.info *.m *.nml *.out
	rm -rf $(patsubst %,o.%,$(solutes))

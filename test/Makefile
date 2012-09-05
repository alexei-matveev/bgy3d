#
# Use  current working  directory by  default.  WARNING:  "make clean"
# will "rm -rf $(WORK_DIR)" so be carefull:
#
WORK_DIR ?= .

# Shell
SHELL = /bin/sh

#
# Flags for bgy3d executable:
#
CMDFLAG = \
	-N 32 \
	-rho 0.018 \
	-beta 1.1989 \
	-norm-tol 1.0e-2 \
	-max-iter 10000 \
	-L 10.0 \
	-zpad 10.0 \
	-damp-start 1.0 \
	-lambda 0.02

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
	HCl_moments.diff \
	HClM_moments.diff

#
# By default "make" all the diffs:
#
all: $(diffs)

#
# Prerequisites for a  diff are the two files, keep  them in the order
# old, new:
#
%.diff: out/%.out %.out
	diff $(^)

#
# Dont delete *.out files  rebuilt as intermediate targets.  Running a
# test takes some time:
#
.PRECIOUS: %_moments.out

#
# Nothing to  be done for  reference outputs.  Without this  rule Make
# will try to  re-build the reference outputs too  when the executable
# gets changed:
#
out/%_moments.out:
	 # nothing to be done, is there a better way to say that?

#
# To produce an XYZ_moments.out  file execute ./runqm.sh XYZ. Updating
# the executable  should lead to re-running  the test, so add  it as a
# pre-requisite:
#
%_moments.out: $(EXE)
	$(SHELL) runtest.sh $(*)

BGY2site:
	(cd ${WORK_DIR}; ${EXE} -BGY2site ${CMDFLAG})

BGYM2site:
	(cd ${WORK_DIR}; ${EXE} -BGYM2site ${CMDFLAG})

clean:
	rm -rf ${WORK_DIR}

#!/bin/bash

#
# FIXME: will not  work if executed from arbitrary  directory. We need
# to locate  the python scripts.  All references are relative  to this
# directory:
#
dir=$PWD

# Set python helper function directory
PYDIR=$dir/../python

if [ ! -d $PYDIR ] ; then
    echo "Check the python helper function directory."
    exit 1
fi

# Directory containing standard outputs
ref=$dir/out

# Create links to the solvents g2 binary files if not exist
solute_init(){
    local svtdir curdir
    svtdir=$1
    curdir=$2
    for f in g00.bin g01.bin g11.bin; do
        ln -sf $svtdir/$f $curdir
    done
}



main(){

    local cmd
    cmd=$1

    # Name test directories with initial "test_"
    workname=test_$cmd

    # Use absolute path
    workdir=$dir/$workname

    # Set file name for computing moments and standard output
    resmoments=$dir/$cmd\_moments.out
    stdmoments=$ref/$cmd\_moments.out

    case "$cmd" in

        # 2-site bgy3d, HCl as solvent
        HCl)
            mkdir -p $workdir
            # FIXME: maybe passing command line parameters from here
            # rather than setting in makefile?
            make HCl WORK_DIR=$workdir

            # Get moments information for each particle pair
            echo | tee $resmoments
            echo "Moments for H-H:" | tee -a $resmoments
            python $PYDIR/moments.py $workdir/g00.bin 2>&1 | tee -a $resmoments
            echo | tee -a $resmoments
            echo "Moments for Cl-Cl:" | tee -a $resmoments
            python $PYDIR/moments.py  $workdir/g11.bin 2>&1 | tee -a $resmoments
            echo | tee -a $resmoments
            echo "Moments for H-Cl:" | tee -a $resmoments
            python $PYDIR/moments.py  $workdir/g01.bin 2>&1 | tee -a $resmoments
            ;;

        # 2-site bgy3dM, HCl as both solvent and solute
        HClM)
            mkdir -p $workdir
            # Check g2 function data first
            solute_init $dir/test_HCl $workdir
            make HClM WORK_DIR=$workdir

            # Get moments information for each particle pair
            echo | tee $resmoments
            echo "Moments for H-HCl:" | tee -a $resmoments
            python ${PYDIR}/moments.py $workdir/g0.bin 2>&1 | tee -a $resmoments
            echo | tee -a $resmoments
            echo "Moments for Cl-HCl:" | tee -a $resmoments
            python ${PYDIR}/moments.py $workdir/g1.bin 2>&1 | tee -a $resmoments
            # mv temp files
            rm -f ./vec*.dat
            ;;

        clean)
            if [[ -n `ls | grep test_` ]]; then
                workdir=$(ls -d test_* | awk -F" " '{print $1}')
                for i in $(seq ${#workdir[@]}); do
                    make clean WORK_DIR=${workdir[$i]}
                done
            else
                echo "No test directory"
            fi
            exit 1
            ;;

        help)
            echo "HCl   :   2-site bgy3d, HCl as solvent"
            echo "HClM  :   2-site bgy3dM, HCl as both solvent and solute"
            echo "clean :   Clean test directories"
            ;;

        *)
            echo "Unknown options"
            echo "type './runtest.sh help' for more information"
            exit 1
            ;;

    esac
}

main $*

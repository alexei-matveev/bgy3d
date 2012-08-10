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

# Set path for bgy3d executable
bgyexe=$dir/../bgy3d

if [ ! -x $bgyexe ] ; then
    echo "Check the bgy3d executable."
    exit 1
fi

# Directory containing standard outputs
STDOUT=$dir/out

# Comparing two files
function diff_out(){
    local target std
    target=$1
    std=$2
    if [ -f $target -o -f $std ]; then
        if (diff $target $std) then
            echo "No differences!"
        else
            echo -e "Possible problems with $target, see diffs above."
        fi
    else
        echo -e "Cannot find $target/$std."
        exit 1
    fi
}

# Create links to the solvents g2 binary files if not exist
function solute_init(){
    local svtdir curdir
    svtdir=$1
    curdir=$2
    if [ ! -f $curdir/g2H.bin -o ! -f $curdir/g2O.bin -o ! -f $curdir/g2HO.bin ]; then
        if [ -f $svtdir/g2H.bin -o -f $svtdir/g2O.bin -o -f $svtdir/g2HO.bin ]; then
            ln -sf $svtdir/g2H.bin $curdir/g2H.bin
            ln -sf $svtdir/g2O.bin $curdir/g2O.bin
            ln -sf $svtdir/g2HO.bin $curdir/g2HO.bin
        else
            echo "Error: g2 function data cannot be found, run test for solvent only first."
            exit 1
        fi
    fi
}



function main(){

    local cmd
    cmd=$1

    # Name test directories with initial "test_"
    workname=test_$cmd

    # Use absolute path
    workdir=$dir/$workname

    # Set output file
    testout=$workdir/$workname.out


    # Set file name for computing moments and standard output
    resmoments=$workdir/$cmd\_moments.out
    stdmoments=$STDOUT/$cmd\_moments.out

    case "$cmd" in

        # 2-site bgy3d, HCl as solvent
        HCl)
            echo -e "Run $workname:"
            mkdir -p $workdir
            # FIXME: maybe passing command line parameters from here rather than setting in makefile?
            make -f Makefile BGY2site WORK_DIR=$workdir EXE=$bgyexe  2>&1 | tee $testout

            # Get moments information for each particle pair
            echo "Caculating moments of distribution functions:"
            echo -e "\nMoments for H-H:" | tee  $resmoments
            python $PYDIR/testmoments.py  $workdir/g2H.bin 2>&1 | tee -a $resmoments
            echo -e "\nMoments for Cl-Cl:" | tee -a $resmoments
            python $PYDIR/testmoments.py  $workdir/g2O.bin 2>&1 | tee -a $resmoments
            echo -e "\nMoments for H-Cl:" | tee -a $resmoments
            python $PYDIR/testmoments.py  $workdir/g2HO.bin 2>&1 | tee -a $resmoments
            ;;

        # 2-site bgy3dM, HCl as both solvent and solute
        HClM)
            echo -e "Run $workname:"
            mkdir -p $workdir
            # Check g2 function data first
            solute_init $dir/test_HCl $workdir
            make -f Makefile BGYM2site WORK_DIR=$workdir EXE=$bgyexe 2>&1 | tee $testout

            # Get moments information for each particle pair
            echo "Caculating moments of distribution functions:"
            echo -e "\nMoments for H-HCl:" | tee $resmoments
            python ${PYDIR}/testmoments.py $workdir/vecH-0.m 2>&1 | tee -a $resmoments
            echo -e "\nMoments for Cl-HCl:" | tee -a $resmoments
            python ${PYDIR}/testmoments.py $workdir/vecO-0.m 2>&1 | tee -a $resmoments
            # mv temp files
            rm -rf ./vec*.dat
            ;;

        clean)
            if [[ -n `ls | grep test_` ]]; then
                workdir=(`ls -d test_* | awk -F" " '{print $1}'`)
                for((i=0;i<${#workdir[@]};i++));do
                    make -f Makefile clean WORK_DIR=${workdir[$i]}
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

    # Check differences
    echo -e "\nComparing with verified outputs:"
    diff_out $resmoments $stdmoments
}

main $*




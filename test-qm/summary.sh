#!/bin/bash

grep_match () {

    # $1 is grep pattern
    local greppattern=$1

    # $2 is options for grep, e.g. -A 1
    local opt=$2
    if [ -n "$greppattern" ]; then

      local tmp=`grep $opt "$greppattern" $file | tail -1`

      if [ -n "$tmp" ]; then

	# text seperated with space, so print the last field after "="
	# FIXME: what if no space between "=" with target value
	echo $tmp | awk '
		    BEGIN {FS="[ ]+"}
			  {
			    for (i = 1; i <= NF; i++) 
			      if ($i~"=") 
				j=i+1 
			  }
		    # only print six digits after the decimal point
		    END {printf "%.6f\n", $j}'
      else
	echo "entry not found"
	exit 1
      fi
    else
      echo "no entry to match"
      exit 1
    fi
}

main () {

    # receive file name from input
    file=$1

    # Get charge
    echo "****charge information****"

    q_u=$(grep_match "integrated charge")
    echo "Solute charge:" $q_u

    q_v=$(grep_match "charge of solvent medium")
    echo "Induced charge:" $q_v

    # Get dipole
    echo "****dipole information****"

    d_u=$(grep_match "dipole moment of solute")
    echo "Solute dipole:" $d_u

    d_v=$(grep_match "dipole moment of solvent")
    echo "Solute dipole:" $d_v

    # Get energy
    echo "****Solvation energy information****"

    # Interaction between point nuclei with solvent potential
    bgy_e_nuc=$(grep_match "with solute point nuclei")
    echo "e_nuc:" $bgy_e_nuc

    # Interaction between diffuse solute charge density with solvent potential
    bgy_e_dif=$(grep_match "diffuse charge density")
    echo "e_dif:" $bgy_e_dif

    # Interaction between solvent charge density with solute electrostatic potential
    bgy_e_slv=$(grep_match "with long-range electrostatic")
    echo "e_slv:" $bgy_e_slv

    # Energy from PG
    e_qm=$(grep_match "Total energy of solute with solvent" "-A 1")
    echo "e_qm :" $e_qm
}

main $*

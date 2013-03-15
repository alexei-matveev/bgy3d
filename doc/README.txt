
QM and MM Solutes in HCl
========================

Options:

 `--N 32 --rho 0.018 --beta 1.1989 --norm-tol 1.0e-7 --max-iter 1000 --L
 10.0 --damp-start 1.0 --lambda 0.02 --snes-solver jager --save-guess
 --solute ...`

This corresponds to

- $V$ = $20^3$ \AA$^3$
- $\beta$ = 1.1989 (420.0 K)
- $\rho$ = 0.018 (144 per cell)

Hydrogen chloride (HCl)
-----------------------

Properties  of hydrogen chloride  in HCl  solvent. Charges  in natural
units e, dipole moments in e\AA, energies and potentials in kcal/mol.

Property               MM, 32     MM, 64     QM, 64          SCF, 64
-----------------   ---------  ---------  ---------  ---------------
Solute dipole        0.251400   0.251400   0.358877         ?
Induced dipole       0.250201   0.229164   0.308945         0.333666
Solute charge        0.0        0.0       -0.021331        -0.021161
Induced charge      -0.040317  -0.051056  -0.024913        -0.023833
Solvation:
point nuclei        -0.842385  -0.833298  40.076264        40.590794
diffuse density     -0.834841  -0.833376  -2.011084        -2.254806
solvent density     -0.840428  -0.837090  -2.016555        -2.261501
QM                          -          -  -1.993208        -2.237503
Electron reorg.             -          -          -         0.258717
Potentials:
H                   -0.939689  -0.898025  -0.348204        -0.651956
Cl                   3.272238   3.268467   5.774924         5.891822
---------------------------------------------------------------------------

There is something wrong with the dipole moments as computed by QM.

The five entries for  "Solvation" are all slightly different estimates
of  the electrostatic  interaction of  the solute  and  solvent charge
distributions:

- The  first  estimate, labeled  with  "point  nuclei"  relies on  the
  interpolated values of the solvent  potential at the location of the
  nuclei.

- The second  estimate, "diffuse density",  is the interaction  of the
  diffuse  solute  charge density  (Gaussian-smeared  cores) with  the
  solvent field

- The  third  entry, "solvent  density",  is  the  interaction of  the
  solvent  charge  density with  the  electrostatic  potential of  the
  diffuse solute density.

- The fourth entry, "QM", is computed using the matrix representations
  of the  electron density and  solvent field as well  as interpolated
  values of the solvent field at locations of nuclei.

- The fifth  entry, "Electron  reorg.", is the  difference of  SCF and
  post-SCF  expectation  values of  the  "gas  phase" hamiltonian.  It
  reflects  the strength  of  the reorganization  of electron  density
  under influence of the solvent field.

Carbon disulfide (CS2)
----------------------

Property               MM, 32     MM, 64     QM, 64          SCF, 64
-----------------   ---------  ---------  ---------  ---------------
Solute dipole        0.0        0.0        0.0
Induced dipole       0.001715   0.000855   0.000857
Solute charge        0.0        0.0        0.018068
Induced charge      -0.084692  -0.125740  -0.142392
Solvation:
point nuclei        -0.372028  -0.373906  14.373086
diffuse density     -0.370351  -0.373662  -0.035914
solvent density     -0.370255  -0.374992  -0.035318
QM                         -          -    0.164666
Electron reorg.            -          -           -
Potentials:
---------------------------------------------------------------------------

Water (H2O)
-----------

Property               MM, 32     MM, 64     QM, 64          SCF, 64
-----------------   ---------  ---------  ---------  ---------------
Solute dipole        0.488557   0.488557   0.523955
Induced dipole       0.482920   0.431496   0.416620
Solute charge        0.0        0.0        0.164691
Induced charge      -0.004372   0.002604  -0.143676
Solvation:
point nuclei        -5.711713  -5.489982  -106.9138
diffuse density     -5.515792  -5.488975  -8.239829
solvent density     -5.539995  -5.505404  -8.091826
QM                          -          -  -8.131681
Electron reorg.             -          -          -
Potentials:
---------------------------------------------------------------------------

Methanol (CH3OH)
----------------

Property               MM, 32     MM, 64     QM, 64          SCF, 64
-----------------   ---------  ---------  ---------  ---------------
Solute dipole        0.501684   0.501684   0.452612
Induced dipole       0.554468   0.502932   0.450425
Solute charge        0.0        0.0        0.019721
Induced charge      -0.115578  -0.162064  -0.179812
Solvation:
point nuclei        -5.544705  -5.524472  17.673523
diffuse density     -5.467964  -5.516670  -5.016932
solvent density     -5.500922  -5.538708  -5.037368
QM                          -          -  -4.810845
Electron reorg.             -          -          -
Potentials:
---------------------------------------------------------------------------

Butanoic acid (C3H7COOH)
------------------------

Property               MM, 32     MM, 64     QM, 64          SCF, 64
-----------------   ---------  ---------  ---------  ---------------
Solute dipole        0.617248   0.617248   0.586525
Induced dipole       0.690702   0.690340   0.730560
Solute charge        0.0        0.0       -0.072836
Induced charge      -0.085109  -0.108007  -0.035466
Solvation:
point nuclei        -4.594817  -4.519687  271.32440
diffuse density     -4.542757  -4.514622  -7.136493
solvent density     -4.593047  -4.550776  -7.158122
QM                          -          -  -6.980013
Electron reorg.             -          -          -
Potentials:
---------------------------------------------------------------------------

Hexane (C6H14)
------------------------

Property               MM, 32     MM, 64     QM, 64          SCF, 64
-----------------   ---------  ---------  ---------  ---------------
Solute dipole        0.001209   0.001209   0.008223
Induced dipole       0.015086   0.009558   0.137926
Solute charge        0.0        0.0       -0.010074
Induced charge      -0.035860  -0.042026  -0.085149
Solvation:
point nuclei        -0.070050  -0.027137  82.000621
diffuse density     -0.031091  -0.027207  -0.343400
solvent density     -0.028981  -0.025604  -0.355654
QM                          -          -  -0.068319
Electron reorg.             -          -          -
Potentials:
---------------------------------------------------------------------------

### Raw data QM, 64



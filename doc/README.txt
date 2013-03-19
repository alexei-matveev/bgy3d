
QM and MM Solutes in HCl
========================

Options:

 `--N 32 --rho 0.018 --beta 1.1989 --norm-tol 1.0e-7 --max-iter 1000 --L
 10.0 --damp-start 1.0 --lambda 0.02 --snes-solver jager --save-guess
 --solute ...`

This corresponds to

- $V$ = 20^3^ \AA^3^
- $\beta$ = 1.1989 (420.0 K)
- $\rho$ = 0.018 (144 per cell)

Hydrogen chloride (HCl)
-----------------------

Cl -- H

Properties  of hydrogen chloride  in HCl  solvent. Charges  in natural
units e, dipole moments in e\AA, energies and potentials in kcal/mol.

Property               MM, 32     MM, 64     QM, 64    SCF, 64
-----------------   ---------  ---------  ---------  ---------
Solute dipole        0.251400   0.251400   0.358877   0.384143
Induced dipole       0.250201   0.229164   0.308945   0.333666
Solute charge        0.0        0.0       -0.021331  -0.021161
Induced charge      -0.040317  -0.051056  -0.024913  -0.023833
Solvation:
point nuclei        -0.842385  -0.833298  40.076264  40.590794
diffuse density     -0.834841  -0.833376  -2.011084  -2.254806
solvent density     -0.840428  -0.837090  -2.016555  -2.261501
QM                          -          -  -1.993208  -2.237503
Electron reorg.             -          -          -   0.061836
Potentials:
H                   -0.939689  -0.898025  -0.348204  -0.651956
Cl                   3.272238   3.268467   5.774924   5.891822
--------------------------------------------------------------------

Electronic energies:

* -15.536046511749  (post-scf)
* -15.535947971512  (scf)

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
  post-SCF  expectation  values of  the  "gas  phase" hamiltonian  (or
  rather the  difference of the DFT counterparts  of these expectation
  values).  It reflects the strength of the reorganization of electron
  density under influence of the solvent field.

Carbon disulfide (CS~2~)
------------------------

S1 = C = S2

Property               MM, 32     MM, 64     QM, 64    SCF, 64
-----------------   ---------  ---------  ---------  ---------
Solute dipole        0.0        0.0        0.0        0.000008
Induced dipole       0.001715   0.000855   0.000857   0.000857
Solute charge        0.0        0.0        0.018068   0.018054
Induced charge      -0.084692  -0.125740  -0.142392  -0.142495
Solvation:
point nuclei        -0.372028  -0.373906  14.373086  14.100849
diffuse density     -0.370351  -0.373662  -0.035914  -0.038270
solvent density     -0.370255  -0.374992  -0.035318  -0.037677
QM                         -          -    0.164666   0.162444
Electron reorg.            -          -           -   0.000744
Potentials:
C		     4.627730
S1		     3.420115
S2		     3.419578
--------------------------------------------------------------

Electronic energies:

* -25.969418157956 (post-scf)
* -25.969416971661 (scf)

Water (H~2~O)
-------------

  H1
 /
O
 \
  H2

Property               MM, 32     MM, 64     QM, 64    SCF, 64
-----------------   ---------  ---------  ---------  ---------
Solute dipole        0.488557   0.488557   0.523955   0.555038
Induced dipole       0.482920   0.431496   0.416620   0.442361
Solute charge        0.0        0.0        0.164691   0.165793
Induced charge      -0.004372   0.002604  -0.143676  -0.139845
Solvation:
point nuclei        -5.711713  -5.489982  -106.9138  -109.3312
diffuse density     -5.515792  -5.488975  -8.239829  -9.184369
solvent density     -5.539995  -5.505404  -8.091826  -9.033639
QM                          -          -  -8.131681  -9.078939
Electron reorg.             -          -          -   0.198669
Potentials:
O		     4.263426
H1		    -2.585132
H2		    -2.585170
--------------------------------------------------------------

Electronic energies:

* -17.154353092917 (post-scf)
* -17.154036496839 (scf)

Methanol (CH~3~OH)
------------------

  H2 H3    H4
    \|    /
     C---O
    /
   H1

Property               MM, 32     MM, 64     QM, 64    SCF, 64
-----------------   ---------  ---------  ---------  ---------
Solute dipole        0.501684   0.501684   0.452612   0.503349
Induced dipole       0.554468   0.502932   0.450425   0.505802
Solute charge        0.0        0.0        0.019721   0.017881
Induced charge      -0.115578  -0.162064  -0.179812  -0.192100
Solvation:
point nuclei        -5.544705  -5.524472  17.673523  17.097857
diffuse density     -5.467964  -5.516670  -5.016932  -6.381589
solvent density     -5.500922  -5.538708  -5.037368  -6.406807
QM                          -          -  -4.810845  -6.158992
Electron reorg.             -          -          -   0.277062
Potentials:
C		     2.278130
H1		     3.190031
H2		    -3.398918
H3		     3.206860
O		     6.828317
H4		    -3.184719
--------------------------------------------------------------

Electronic energies:

* -23.983214454853 (post-scf)
* -23.982772933117 (scf)

Butanoic acid (C~3~H~7~COOH)
----------------------------

H6 H5  H2 H1
  \|    \|
  C4    C2
  / \   / \   O2--H8
H7   \ /   \ /
      C3    C1
      |\    ||
     H3 H4  O1

Property               MM, 32     MM, 64     QM, 64    SCF, 64
-----------------   ---------  ---------  ---------  ---------
Solute dipole        0.617248   0.617248   0.586525   0.694934
Induced dipole       0.690702   0.690340   0.730560   0.850978
Solute charge        0.0        0.0       -0.072836  -0.075419
Induced charge      -0.085109  -0.108007  -0.035466  -0.032293
Solvation:
point nuclei        -4.594817  -4.519687  271.32440  288.78482
diffuse density     -4.542757  -4.514622  -7.136493  -8.709857
solvent density     -4.593047  -4.550776  -7.158122  -8.743076
QM                          -          -  -6.980013  -8.565302
Electron reorg.             -          -          -   0.386236
Potentials:
C1		     4.880232
O1		    11.008162
O2		     2.791229
C2		     2.235826
C3		     4.282062
C4		     1.943947
OH		    -1.490686
H2		    -1.273514
H3		     2.052024
H4		     6.580057
H5		     4.709021
H6		     0.541590
H7		     1.494241
H8		     2.451006
--------------------------------------------------------------

Electronic energies:

* -59.326272328142 (post-scf)
* -59.325656827000 (scf)

Hexane (C~6~H~14~)
------------------
 H11 H10 H7 H6 H3
    \ |   \ |   \
     C5    C3    C1--H1
H14 / \   / \   / \
 \ /   \ /   \ /   H2
  C6    C4    C2
  | \   | \   | \
H12 H13 H8 H9 H4 H5

Property               MM, 32     MM, 64     QM, 64    SCF, 64
-----------------   ---------  ---------  ---------  ---------
Solute dipole        0.001209   0.001209   0.008223   0.006491
Induced dipole       0.015086   0.009558   0.137926   0.153552
Solute charge        0.0        0.0       -0.010074  -0.010583
Induced charge      -0.035860  -0.042026  -0.085149  -0.087341
Solvation:
point nuclei        -0.070050  -0.027137  82.000621  81.392287
diffuse density     -0.031091  -0.027207  -0.343400  -0.372157
solvent density     -0.028981  -0.025604  -0.355654  -0.385512
QM                          -          -  -0.068319  -0.089780
Electron reorg.             -          -          -   0.007266
Potentials:
C1		     2.358558
C2		     2.443140
C3		     2.489955
C4		     2.462004
C5		     2.421349
C6		     2.346722
H1		     2.243433
H2		     2.245605
H3		     2.266385
H4		     2.360914
H5		     2.359280
H6		     2.449220
H7		     2.452113
H8		     2.424263
H9		     2.422642
H10		     2.364183
H11		     2.365350
H12		     2.229591
H13		     2.228886
H14		     2.169376
--------------------------------------------------------------

Electronic energies:

* -42.289934379628 (post-scf)
* -42.289922799953 (scf)


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

Hydrogen chloride
-----------------

Properties  of hydrogen chloride  in HCl  solvent. Charges  in natural
units e, dipole moments in e\AA, energies and potentials in kcal/mol.

Property               MM, 32     MM, 64           QM, 64          SCF, 64
-----------------   ---------  ---------  ---------------  ---------------
Solute dipole        0.251400   -                 ?         5.900857?
Induced dipole       0.250201   0.229164   0.308945         0.333666
Solute charge        0.0        -         -0.021331        -0.021161
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

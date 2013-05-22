
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

~~~~~~~~~~~~~~~~~~~~~~~
Cl -- H
~~~~~~~~~~~~~~~~~~~~~~~

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

~~~~~~~~~~~~~~~~~~~~~~~
S1 = C = S2
~~~~~~~~~~~~~~~~~~~~~~~

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
C                    4.627730   4.646586   0.661360   0.632495
S1                   3.420115   3.432765   0.977444   0.964449
S2                   3.419578   3.432447   0.977164   0.964030
--------------------------------------------------------------

Electronic energies:

* -25.969418157956 (post-scf)
* -25.969416971661 (scf)

Water (H~2~O)
-------------

~~~~~~~~~~~~~~~~~~~~~~~
  H1
 /
O
 \
  H2
~~~~~~~~~~~~~~~~~~~~~~~

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
O                    4.263426   4.064380  -11.61762  -11.77347
H1                  -2.585132  -2.518386  -18.60392  -19.34492
H2                  -2.585170  -2.518279  -18.60420  -19.34538
--------------------------------------------------------------

Electronic energies:

* -17.154353092917 (post-scf)
* -17.154036496839 (scf)

Methanol (CH~3~OH)
------------------

~~~~~~~~~~~~~~~~~~~~~~~
  H2 H3    H4
    \|    /
     C---O
    /
   H1
~~~~~~~~~~~~~~~~~~~~~~~

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
C                    2.278130   2.105761   0.584761   0.649005
H1                   3.190031   3.133311   1.686232   2.103915
H2                  -3.398918  -3.477980  -4.257197  -4.986875
H3                   3.206860   3.130518   1.905240   2.346881
O                    6.828317   6.437756   3.724602   3.861846
H4                  -3.184719  -3.694387  -6.347412  -8.133162
--------------------------------------------------------------

Electronic energies:

* -23.983214454853 (post-scf)
* -23.982772933117 (scf)

Butanoic acid (C~3~H~7~COOH)
----------------------------

~~~~~~~~~~~~~~~~~~~~~~~
H6 H5  H2 H1
  \|    \|
  C4    C2
  / \   / \   O2--H8
H7   \ /   \ /
      C3    C1
      |\    ||
     H3 H4  O1
~~~~~~~~~~~~~~~~~~~~~~~

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
C1                   4.880232   4.723773   9.625472  10.434688
O1                  11.008162  10.820006  17.212581  19.333331
O2                   2.791229   2.707382   9.175892   9.578772
C2                   2.235826   2.032847   4.354700   4.308431
C3                   4.282062   4.027545   5.327668   5.698965
C4                   1.943947   1.873328   2.057880   1.689155
OH                  -1.490686  -1.504560   5.103405   5.352918
H2                  -1.273514  -1.400071   0.756783   0.016382
H3                   2.052024   1.933830   3.224199   2.935339
H4                   6.580057   6.548505   7.706036   8.692858
H5                   4.709021   4.633701   6.519890   7.141003
H6                   0.541590   0.440565   0.684162  -0.107680
H7                   1.494241   1.407662   1.233492   0.670217
H8                   2.451006   2.455203   2.302711   2.086211
--------------------------------------------------------------

Electronic energies:

* -59.326272328142 (post-scf)
* -59.325656827000 (scf)

Hexane (C~6~H~14~)
------------------

~~~~~~~~~~~~~~~~~~~~~~~
 H11 H10 H7 H6 H3
    \ |   \ |   \
     C5    C3    C1--H1
H14 / \   / \   / \
 \ /   \ /   \ /   H2
  C6    C4    C2
  | \   | \   | \
H12 H13 H8 H9 H4 H5
~~~~~~~~~~~~~~~~~~~~~~~

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
C1                   2.358558   2.263026   3.750902   3.884198
C2                   2.443140   2.326814   2.914540   2.976435
C3                   2.489955   2.313950   1.997110   1.963126
C4                   2.462004   2.293129   1.481896   1.410374
C5                   2.421349   2.275606   1.177873   1.064181
C6                   2.346722   2.210454   1.474596   1.386596
H1                   2.243433   2.207931   3.994543   4.165565
H2                   2.245605   2.209487   4.005858   4.178340
H3                   2.266385   2.212073   3.694467   3.798112
H4                   2.360914   2.303916   2.947476   3.019432
H5                   2.359280   2.302221   2.944283   3.014983
H6                   2.449220   2.292926   1.870753   1.817494
H7                   2.452113   2.294786   1.879893   1.828096
H8                   2.424263   2.258943   1.414856   1.348360
H9                   2.422642   2.257781   1.413369   1.346068
H10                  2.364183   2.266989   1.040165   0.905699
H11                  2.365350   2.267791   1.044197   0.910570
H12                  2.229591   2.189278   1.489301   1.410292
H13                  2.228886   2.188902   1.488162   1.408692
H14                  2.169376   2.134118   1.585625   1.500951
--------------------------------------------------------------

Electronic energies:

* -42.289934379628 (post-scf)
* -42.289922799953 (scf)


1d- and 3d-RISM with HNC closure
================================

Analytical vs. FFT expression for real space long-range Coulomb
---------------------------------------------------------------

Using analytical  expression for the real space  representation of the
long-range  Coulomb  potential and  the  FFT  transform  of the  (also
analytical)  Fourier space  representation of  the  long-range Coulomb
potential  will  not   necessarily  produce  identical  results.   The
real-space representation of the  long-range Coulomb potential is used
to calculate  the long-range part  of the direct  correlation function
and thus contributes to the  excess chemical potential.  In 1d-RISM we
compute the  finite grid real space representation  from an analytical
expression. In 3d-RISM we use the FFT transform instead.

The excess  chemical potential is  calculated by the  volume integral,

$$
\beta\mu =
         \rho\int \left[\frac{1}{2}h^2(\vec{r})
          - c(\vec{r})
          - \frac{1}{2}h(\vec{r})c(\vec{r})\right] d^3r,
$$

to  get the  long-range component  of $c(\vec{r})$  term, we  need the
long-range tail of Coulomb potential, which could be obtained from two
sources: analytical expression and real-space represention of FFT.

Here we  list the results of  two solvents: simple  charged LJ solvent
and modified  TIP3P water.  It could be  obserbed that  the difference
between the results of anaylical and FFT coulomb potential is tiny.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ("LJC" (("LJ" (0.0 0.0 0.0) 1.0 1.0 -0.1)
         ("LJ" (0.0 0.0 0.5) 0.5 1.0 +0.1)))

 ("water"
   (("O" (-0.2929 0.0 0.0) 3.1506 0.1521 -0.834)
    ("OH" (0.2929 0.757 0.0) 0.4 0.046 0.417)
    ("OH" (0.2929 -0.757 0.0) 0.4 0.046 0.417)))
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Excess chemical potential (kcal):

Solvent            Analytical                  FFT
-------           -----------           ----------
LJC                16.292235             16.292285
Water              -6.320048             -6.309774
---------------------------------------------------


Excess chemical potential of 1D- and 3D-RISM for pure solvents
--------------------------------------------------------------

This section  contains the excess chemical  potential results obtained
from our  implementation of 1D-  and 3D-RISM model. Besides  "LJC" and
modified TIP3P water solvent used in last section, we also compare the
values of several simple solvents:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
("LJ" (("LJ" (0.0 0.0 0.0) 1.0 1.0 0.0)))

("LJ2" (("LJ" (0.0 0.0 0.0) 1.0 1.0 0.0)
        ("LJ" (0.0 0.0 0.5) 0.5 1.0 0.0)))

("OW" (("OW" (0.0 0.0 0.0) 3.16 0.1549 0.0)))
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From the tables  below, we could see that  numbers solved from 3D-RISM
model are close to their counterparts in 1D-RISM model.

Excess chemical potential (kcal):

Solvent              3D-RISM               1D-RISM
-------           ----------            ----------
LJ                 15.580806             15.854951
LJ2                16.476489             16.542005
LJC                16.292285             16.385851
OW                  6.910124              6.940968
Water              -6.309774             -6.303708
---------------------------------------------------


Excess chemical potential in 1D- and 3D-RISM for solute/solvent systems
-----------------------------------------------------------------------

Excess chemical potential  $\mu$ (in kcal/mol) solute/solvent mixtures
at infinite  dilutions computed by  1d- and 3d-RISM.   Excess chemical
potential for pure  solvens computed with the same  settings is quoted
for  reference. Most  of  the calculations  involving charged  solvent
sites in realisitc water model diverge when employing HNC closure with
the default Newton solver (indicated by "-").

Species               in       1D-HNC         3D-HNC
-------            ----- ------------   ------------
Pure solvent:
Ow                           6.922977       6.930328
Nw                           6.981084       6.973021
Water                       -6.262289      -6.281097
Solute/Solvent:
Ow                  Ow       6.922977       6.930328
                    Nw       6.550974       6.542975
                    Water    5.980016       5.943212
Nw                  Ow       7.087025       6.604694
                    Nw       6.981084       6.234882
                    Water    6.627866       5.657242
Water               Ow       7.087025       6.604694
                    Nw       6.981084       6.234882
                    Water   -6.262280              -
Hydrogen chloride   Ow       8.057284       7.801370
                    Nw       7.588731       7.001046
                    Water    6.009940              -
Methanol            Ow      13.621287      11.399124
                    Nw      14.078878      10.526617
                    Water    1.432886              -
Butanoic acid       Ow      33.758116      27.205191
                    Nw      36.250502      24.535412
                    Water   25.450272              -
Carbon disulfide    Ow      15.511894      14.193590
                    Nw      14.597500      12.476879
                    Water   10.910771              -
Hexane              Ow      42.384342      32.110759
                    Nw      46.609190      29.196929
                    Water   44.770328              -
----------------------------------------------------

Here "Ow"  is a  single-site LJ water  model derived from  SPC/E water
model  by removing  the hydrogens  and  setting the  oxygen charge  to
zero. "Nw" is a three-site water model derived from the modified TIP3P
water  model by  setting all  site charges  to zero.   "Water"  is the
modified  TIP3P  water model.   See  `guile/solutes.scm` for  details.
FIXME: As it stands "Ow" does not directly derive from "Nw".

The numbers in this section were  obtained with 1d-RISM (N = 1024) and
3d-RISM (N = 128) and the rest of the settings like this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--N $N --L 10.0 --beta 1.6889
--max-iter 1000 --rho 0.033427745 --norm-tol 1e-12
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

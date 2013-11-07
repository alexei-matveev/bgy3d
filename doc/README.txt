
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
LJ                 15.725374             15.808186
LJ2                16.593576             16.579971
LJC                16.473840             16.435115
OW                  6.930328              6.922977
Water              -6.281097             -6.262288
---------------------------------------------------

Except setting "--rho 0.75" and "--beta 0.25" for "LJ", "LJ2" and "LJC"
solvents, all the other calculation parameters for the comparison above
could be found in the next section in which excess chemical potentials
for solute/solvent systems are calcualted.


Excess chemical potential in 1D- and 3D-RISM for solute/solvent systems
-----------------------------------------------------------------------

Excess chemical potential  $\mu$ (in kcal/mol) solute/solvent mixtures
at infinite  dilutions computed by  1d- and 3d-RISM.   Excess chemical
potential for pure  solvens computed with the same  settings is quoted
for  reference. Most  of  the calculations  involving charged  solvent
sites in realisitc water model diverge when employing HNC closure with
the default Newton solver (indicated by "-").

Table 1. <!-- tab:1d3d -->

Species                HNC1      (GF)     HNC3     (GF)      KH1      (GF)       KH3      (GF)
-------      -----  ------- ---------  -------  -------  -------  --------  --------  --------
Solvent
O~w~                  6.923    4.550     6.930    4.610    7.511     6.489     7.498     6.484
X~2~Y                 6.981    2.458     6.973    2.487    7.521     4.556     7.502     4.540
H~2~O                -6.262  -12.097    -6.281  -12.126   -5.767    -9.562    -5.787    -9.587
Solute        in
O~w~          O~w~    6.923    4.550     6.930    4.610    7.511     6.489     7.498     6.484
              X~2~Y   6.551    3.535     6.543    3.563    7.093     5.380     7.075     5.366
              H~2~O   5.980    3.106     5.943    3.060    6.495     4.571     6.457     4.530
X~2~Y         O~w~    7.087    3.824     6.605    4.267    7.671     5.963     7.167     6.159
              X~2~Y   6.981    2.458     6.235    3.245    7.521     4.556     6.761     5.067
              H~2~O   6.628    1.889     5.657    2.767    7.131     3.670     6.162     4.253
H~2~O         O~w~    7.087    3.824     6.605    4.267    7.671     5.963     7.167     6.159
              X~2~Y   6.981    2.458     6.235    3.245    7.521     4.556     6.761     5.067
              H~2~O  -6.262  -12.097         -        -   -5.767    -9.562    -3.860    -5.835
^a)^          H~2~O                                                           -5.632    -7.667
HCl           O~w~    8.057    4.471     7.801    4.389    8.898     6.977     8.756     7.364
              X~2~Y   7.589    2.511     7.001    2.601    8.373     5.005     7.902     5.519
              H~2~O   6.010    0.692         -        -    6.779     2.874     6.006     3.382
^a)^          H~2~O                                                            4.554     1.931
CH~3~OH       O~w~   13.621    8.259    11.399    6.991   14.579    10.516    12.742    11.010
              X~2~Y  14.079    5.015    10.527    4.870   14.975     7.336    11.806     8.758
              H~2~O   1.433   -9.721         -        -    2.351    -6.942     0.997    -2.515
^a)^          H~2~O                                                            1.808    -1.734
C~3~H~7~COOH  O~w~   33.758   19.064    27.205   16.371   36.190    22.893    31.239    27.548
              X~2~Y  36.251    8.699    24.535   10.573   38.502    12.858    28.462    21.525
              H~2~O  25.450   -6.318         -        -   27.818    -1.944    17.466     9.823
^a)^          H~2~O                                                           16.244     8.508
CS~2~         O~w~   15.512    9.915    14.194    8.471   17.026    13.438    16.115    13.852
              X~2~Y  14.598    6.008    12.477    4.945   16.030     9.411    14.323    10.268
              H~2~O  10.911    1.938         -        -   12.365     4.854    11.499     6.988
^a)^          H~2~O                                                           12.366     7.765
C~6~H~14~     O~w~   42.384   23.272    32.111   19.496   45.177    27.220    36.854    32.673
              X~2~Y  46.609    9.636    29.197   13.132   49.167    14.167    33.814    25.907
              H~2~O  44.770    2.118         -        -   47.536     7.115    29.083    20.433
^a)^          H~2~O                                                           29.367    20.562
-----------------------------------------------------------------------------------------------

Here O~w~  is a  single-site LJ water  model derived from  SPC/E water
model  by removing  the hydrogens  and  setting the  oxygen charge  to
zero.  X~2~Y is  a three-site  water model  derived from  the modified
TIP3P water model  by setting all site charges to  zero.  H~2~O is the
modified  TIP3P  water model.   See  `guile/solutes.scm` for  details.
FIXME: As it stands O~w~ does not directly derive from X~2~Y.

The numbers in this section were  obtained with 1d-RISM (N = 1024) and
3d-RISM (N = 128) and the rest of the settings like this:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--N $N --L 10.0 --beta 1.6889
--max-iter 1000 --rho 0.033427745 --norm-tol 1e-12
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rows  labeled with  ^a)^ correspond  to  the water  solvent around  an
unpolarized  (gas-phase)  QM  solute.   These calculations  have  been
performed with N=96 and correspond to what we called post-SCF above.

Influence of grid resolution and box size on RISM results for water
-------------------------------------------------------------------

This is not  really a 3D model, rather 1D-RISM  with 3D machinery. All
distribution functions  here are spherically  symmetric as far  as the
cubic domain shape allows that.  Also all the local singularities that
are  often  present  in  3D  models are  less  pronounced  here  since
effectively the  method operates with  distribution functions averaged
over angular  coordinates.  Compare to the  third row of  the table in
the previous section.

Here we listed the value changes of excess chemical potential of water
solvent  from  3D-RISM  results  with different  grid  numbers(N)  and
lengths(L).   The  default "--snes-solver  newton  "  was switched  to
"--snes-solver trial"  for HNC closure to  guarantee convergence. Note
that grid resolutions are identical for diagonal  pairs (h = 2 * L / N
= 0.15625).

Due to the limit of  computation resources, the maximum grid number we
used in this  test is 224 for each  dimension. Calculations using more
192  grid numbers had  already begun  to consume  swap on  single node
(neh13 and neh19).

HNC closure:

 N      L=10.0    L=12.5    L=15.0    L=17.5
----    ------    ------    ------    ------
128     -6.281    -6.309    -6.334    -6.328
160     -6.279    -6.317    -6.330    -6.327
192     -6.280    -6.316    -6.332    -6.344
224     -6.280    -6.316    -6.331    -6.340
--------------------------------------------

KH closure:

water(solvent)/water(solvent)

 N      L=10.0    L=12.5    L=15.0    L=17.5
----    ------    ------    ------    ------
128     -5.787    -5.808    -5.835    -5.828
160     -5.783    -5.823    -5.836    -5.828
192     -5.784    -5.821    -5.839    -5.849
224     -5.785    -5.821    -5.837    -5.847
--------------------------------------------

water(solvent)/water(solute)

 N      L=10.0    L=12.5    L=15.0    L=17.5
----    ------    ------    ------    ------
128     -3.860    -3.880    -3.909    -3.920
160     -3.860    -3.895    -3.926    -3.913
192     -3.861    -3.897    -3.911    -3.923
224     -3.861    -3.898    -3.913    -3.919
--------------------------------------------

Solvent electrostatic field with and without boundary condition
---------------------------------------------------------------

When solving Poisson equation  to get the electrostatic potential from
solvent  charge density, we  supply an  option "--no-cage"  to control
whether  not  or zeroing  the  potential  values  at the  boundary  of
simulation box. Potential corrected by boundary condition is in effect
a  superposition of  the  solvent electrostatic  field  and a  surface
charge on that metallic cage.

Here we listed the results of the integration of solvent electrostatic
field with  solute point nuclei $\langle U_v|\rho_N  \rangle$ and with
diffuse charge  density of solute $\langle  U_v|\rho_u \rangle$, which
are the  only two  output entries that  changed when running  with and
without boundary condition.

Calculation parameters used are the same with those in section "Excess
chemical potential in 1D- and 3D-RISM for solute/solvent systems".


Solute        cage    $<U_v|\rho_N>$    $<U_v|\rho_u>$
------        ----    --------------    --------------
H~2~O           no        -22.830917        -22.152216
               yes        -22.816307        -22.137614
HCl             no         -1.894681         -1.892238
               yes         -1.891148         -1.888711
CH~3~OH         no        -21.772003        -21.176797
               yes        -21.756322        -21.161100
C~3~H~7~COOH    no        -15.937807        -15.547307
               yes        -15.896105        -15.505732
CS~2~           no         -2.153173         -2.120207
               yes         -2.149663         -2.116716
C~6~H~14~       no         -0.240059         -0.272374
               yes         -0.239890         -0.272191
------------------------------------------------------


MM and QM solutes in water
--------------------------

### Water in PCM-water

This section discusses COSMO/PCM reference results for water using the
same minimal  CRENBL ECP for  reference purposes. For this  reason the
geometry of the water solute is kept fixed to that of the TIP3P model.
The difference of the two energies,  PCM and gas phase, amounts to the
total  solvation energy  estimate of  -9.58 kcal.   The  solute dipole
moments are also shown:


Method             E~tot~, au    D~u~, eA     D~v~, eA
----------  -----------------  ----------  -----------
qm/gp        -17.154353092917   0.5239553
qm/pcm       -17.169614518847   0.6021832    0.5426852


The next table shows the terms  in the PCM solvation energy as printed
by  a  recent version  of  PG.  The  convention  is  that the  Coulomb
interaction  energy includes the  (Born) factor  1/2.  Note  that this
decomposition  is  missing  a  term  that we  dubbed  "the  electronic
reorganization", a  positive term, here of  10.86 - 9.58  = 1.28 kcal,
which can only be deduced by comparison with the gas phase results.


Term            E, kcal
-------------  --------
Electrostatic    -11.97
Cavitation         5.05
Dispersion        -4.88
Repulsion          0.94
$\Sigma$         -10.86

The  dipole  of  a  gas  phase  molecule  computed  usong  CRENBL  ECP
corresponds to 2.52  D which is significantly more  than 1.85 D quoted
by some  sources [wikipedia].  The dipole  increases to 2.89  D due to
PCM  polarization. However  it is  mostly compensated  by  the induced
medium dipole of 2.61 D of opposite direction.

### MM and QM solutes in KH-water

<!-- Table 1 refers to tab:1d3d -->

The chemical potentials for MM-rows  copied from Table 1 were obtained
with N=128.   FIXME: However  the data  on QM solutes  and all  of the
dipole moments,  including those for  MM-row, were obtained  with N=96
but otherwise the same parameters.  Thus the numbers in MM-row are not
mutually consistent --- the  quality of chemical potential is "better"
than that of the dipole moments.


Solute       Model        KH3      (GF)         E~el~, au    D~u~, eA   D~v~, eA    E~c~    E~c~
-------      -------- -------   -------  ----------------   ---------  ---------  ------  ------
H~2~O        mm        -3.860    -5.835                      0.489         0.455   -11.1
             qm/gp     -5.632    -7.667  -17.154353092917    0.522         0.488   -15.1   -13.0
             qm/aq    -11.576   -13.695  -17.149256062555    0.644         0.604   -22.3   -19.4
             qm/aq/nc -11.582   -13.701  -17.149246364375    0.644         0.604
                                          (3.198)
HCl          mm         6.006     3.382                      0.251         0.225
             qm/gp      4.554     1.931  -15.536046511749    0.358         0.323
             qm/aq      3.966     1.353  -15.535551571194    0.414         0.377
             qm/aq/nc   3.964     1.351  -15.535548519827    0.414         0.377
                                          (0.311)
CH~3~OH      mm         0.997    -2.515                      0.502         0.517
             qm/gp      1.808    -1.734  -23.983214454853    0.454         0.470
             qm/aq     -4.245    -7.883  -23.977789344961    0.615         0.624
             qm/aq/nc  -4.253    -7.891  -23.977775709821    0.615         0.624
                                          (3.404)
C~3~H~7~COOH mm        17.466     9.823                      0.617         0.732
             qm/gp     16.244     8.508  -59.326272328142    0.581         0.641
             qm/aq     ??.??      ??.??  -59.32??????????    ?.??          ?.??
             qm/aq/nc  ??.??      ??.??  -59.32??????????    ?.??          ?.??
                                          (?.???)
CS~2~        mm        11.499     6.988                      0.000         0.035
             qm/gp     12.366     7.765  -25.969418157956    0.000         0.033
             qm/aq     12.347     7.748  -25.969402573454    0.000         0.033
             qm/aq/nc  12.347     7.748  -25.969402441004    0.000         0.033
                                          (0.010)
C~6~H~14~    mm        29.083    20.433                      0.001         0.082
             qm/gp     29.367    20.562  -42.289934379624    0.008         0.077
             qm/aq     29.321    20.518  -42.289896194950    0.009         0.076
             qm/aq/nc  29.320    20.518  -42.289895181470    0.009         0.076
                                          (0.024)
---------------------------------------------------------------------------------

a) The dipoles of the QM solutes have been taken from the medium model
output.  QM  program offers  an analytically integrated  dipole moment
which   is  0.5239553  and   0.6454237  eA   for  gp-   and  aq-water,
respectively. Thus only two digits can be trusted.

b) Note  that 1 D  = 0.20819434  eA.  So that  the MM water  dipole is
about 2.35 D, gas-phase QM water is 2.51 D (???), and aqueous phase QM
water is 3.09 D.

c)  Butanoic acid  in aqueous  medium  did not  converge. Notably  the
dipoles "oscillated" irregularly between ~0.5 and ~1.5 eA.

d)  For  qm/aq/nc calculations  without  boundary condition  (no-cage)
applied to solvent field, all the  values are very close to those with
boundary condition  (qm/aq).  The  boundary condition does  not affect
qm/gp calculations as the action of  the solvent field on QM solute is
ignored.

e) Numbers in parens in the  E~el~ column is the difference of the two
electronic energies  (qm/aq and qm/gp)  in kcal. This is  the electron
reorganization energy due to electric field of the solvent.  The SPC/E
paper quotes  an estimate  for water models:  "With a value  [...]  of
0.001445  nm^3^ for  the polarizability  $\alpha/4\pi  \epsilon_0$, or
1.608  $\times$ 10^-40^ F  m for  $\alpha$,the average  induced dipole
moment  in  the SPC  model  corresponds  to  3.74 kJ/mol  polarization
energy." FIXME: Hm, this is not even close to 3 kcals that we get.

f) The two E~c~ values, both in kcals, correspond to the electrostatic
energy  and include  the factor  1/2. One  is obtained  by integrating
"solvent   electrostatic  field   with  diffuse   charge   density  of
solute". The  other is computed by QM  code by forming a  trace of the
density matrix with the matrix representatio of the potential.


Reference values for execess chemical potential
-----------------------------------------------

For comparison, we collected values for excess chemical potential
reported in others' work in which results for solutes that we used for
regression tests are available. There's no guarantee that our numbers
must be the same with other sources because force field and geometry for
solutes in our code are different with theirs, but good agreement might
indicate that we are not doing things totally wrong.

In source [1], "UC" is short for "universal correction" and in source
[2] and [3], "SDC" is for "structural descriptors correction".  They are
actually different names for the same method which they developed to
"improve the accuracy of hydration free energies".



source [1]:

Solute/FF Solvent/FF Method     KH3  KH3/UC    GF3  GF3/UC $\Delta G_{exp}$
--------- ---------- ------- ------ ------- ------ ------- ----------------
CH~3~OH/  water/     3D-RISM  1.041  -5.971 -2.399  -6.086  -5.100
OPLS-AA   MSPC-E     -KH
---------------------------------------------------------------------------


source [2]:

Solute/FF Solvent/FF Method     GF1    GF3   SDC1   SDC3 $\Delta G_{exp}$
--------- ---------- ------- ------ ------ ------ ------ ----------------
CH~3~OH/   water/    1D/3D-   -6.34  12.24  -7.49  -5.16  -5.11
Amber      Amber     RISM-KH
C~6~H14/   water/    1D/3D-    2.85  20.41   2.49   2.55   2.50
Amber      Amber     RISM-KH
-------------------------------------------------------------------------

source [3]:

Solute/FF Solvent/FF Method    PW1   GF1  SDC/PW  SDC/GF $\Delta G_{exp}$
--------- ---------- ------- ----- ----- ------- ------- ----------------
CH~3~OH/  water/     1D-RISM -0.07 -5.46   -5.76   -6.02  -5.11
OPLS-AA   MSPC-E     -KH
C~6~H14/  water/     1D-RISM 13.33  5.23    2.89    3.03   2.50
OPLS-AA   MSPC-E     -KH
-------------------------------------------------------------------------

[1] "3DRISM Multigrid Algorithm for Fast Solvation Free Energy
    Calculations", Volodymyr P. Sergiievskyi, and Maxim V. Fedorov,
    J. Chem. Theory Comput., 2012, 8 (6), pp 2062-2070,
    http://dx.doi.org/10.1021/ct200815v

[2] "Hydration Thermodynamics Using the Reference Interaction Site
    Model: Speed or Accuracy?",  A. I. Frolov, E. L. Ratkova,
    D. S. Palmer and M. V. Fedorov, J. Phys. Chem. B, 2011, 115 (19),
    pp 6011-6022, http://dx.doi.org/10.1021/jp111271c

[3] "An Accurate Prediction of Hydration Free Energies by Combination of
    Molecular Integral Equations Theory with Structural Descriptors",
    E. L. Ratkova, G. N. Chuev, V. P. Sergiievskyi and M. V. Fedorov,
    J. Phys. Chem. B, 2010, 114 (37), pp 12068-12079,
    http://dx.doi.org/10.1021/jp103955r


Influence of modified hydrogen parameter
----------------------------------------

In our regression test, solute "methanol" and "butanoic acid" are
directly taken from Jager's code. When examining the force field
parameters listed in these two molecules, we found that most of them are
from OPLSAA force field with the exception that vdw interaction
parameters are modified for hydrogen atom on "Alcohol -OH" and
"Carboxylate acid -COOH":

Site                    sigma(A)    epsilon(kcal/mol)
---------------------  ---------   ------------------
Carboxylic Acid -COOH     0.0             0.0
modified                  3.4             0.046
Alcohol -OH               0.0             0.0
modified                  0.4             0.04
-----------------------------------------------------

These modifications bring different changes in 1D- and 3D-RISM
calculations(units in kcal/mol):

Solute           KH1      KH1^*^      KH3   KH3^*^
------------  ------   ---------   ------  -------
CH~3~OH        2.351   -1045.833    0.997    0.973
C~3~H~7~COOH  27.818   -1206.856   17.466   11.195
--------------------------------------------------

^*^Hydrogen vdw parameters are zero


Alcohols, carboxylate acids and aldehydes in water
--------------------------------------------------

1D and 3D-RISM results of neutral solutes, molecule coordinates are from
NIST database. As already seen in last section, the zero vdw parameters
of hydrogen atom in OPLSAA result in unexpected values in 1D model.

Calculation parameter settting is the same as previous sections, energy
units is "kcal/mol"

Solute                   KH1       (GF)      KH3    (GF)   exp^1^
-----------------  ---------  ---------  -------  ------  -------
Alcohols
CH~3~OH            -1012.245  -1026.191    2.842  -0.644   -5.115
C~2~H~5~OH          -992.651  -1012.942    6.225   1.504   -5.014
C~3~H~7~OH          -974.739   -999.246   11.409   5.574   -4.826
C~4~H~9~OH          -982.240  -1012.521   14.829   7.915   -4.716
C~5~H~11~OH         -974.985  -1011.930   19.347  11.221   -4.474
Carboxylate acids
HCOOH              -1168.869  -1183.655    3.535   0.013   -6.99^2^
CH~3~COOH          -1151.262  -1172.668    7.305   2.566   -6.704
C~2~H~5~COOH       -1156.799  -1186.053   10.960   5.080   -6.475
C~3~H~7~COOH       -1130.024  -1166.344   15.869   8.867   -6.379
C~4~H~9~COOH       -1124.953  -1169.595   20.294  12.107  -6.2^2^
Aldehydes
CH~2~O                 5.961     -1.395    4.145   1.030   -1.7
C~2~H~4~O             13.080      0.175    8.150   3.778   -3.504
C~3~H~6~O             20.547      1.240   12.411   6.882   -3.442
C~4~H~8~O             27.879      1.428   16.463   9.842   -3.176
C~5~H~10~O            34.337      0.864   20.642  12.837   -3.031
-----------------------------------------------------------------

[1] Group contributions to the thermodynamic properties of non-ionic
    organic solutes in dilute aqueous solution, S. Cabani, P. Gianni,
    V. Mollica and L. Lepori, J. Solution Chem. 10, 563 (1981),
    http://dx.doi.org/10.1007/BF00646936

[2] Model for Aqueous Solvation Based on Class IV Atomic Charges and
    First Solvation Shell Effects, C. C. Chambers, G. D. Hawkins, C. J.
    Cramer and D. G. Truhlar, The Journal of Physical Chemistry 100,
    16385 (1996). http://dx.doi.org/10.1021/jp9610776


Choice of hydrogen vdw parameters for alcohols and carboxylate acids
--------------------------------------------------------------------

The zero vdw parameters of hydrogen atom in alcohols and carboxylate
acids will cause the overlap between hyrogen atom and solvent particle
with opposite charge sign. At least we can observe unphysical chemical
potential values from converged 1D-RISM calculations from the previous
section. (FIXME: what about 3D?)

One workaround is assigning non-zero vdw parameters for "polar
hydrogens", with which the site-overlap would be avoided. In this
section we list two options to choose the reasonable non-zero vdw
parameters.

1. use the parameters for hdyrogen atom "Hw" in modified TIP3P water
model [1], which acts as solvent in our study. However, it might not be
a good idea to apply vdw parameter of "Hw" in other molecules. In [1],
Reiher also argues that "...the H-H van der Waals $r_{min}$ for water
hydrogens interacting with each other, 0.45 \AA, is too small when used
in protein-water and protein-protein interactions. ...The H-H $r_{min}$
value used in CHARMM for all polar hydrogens, 1.60 \AA, was found to be
a suitable value for describing H-H protein-water van der Waals
interactions." (page 4-19)

 $\sigma$(\AA)    $\epsilon$(kcal/mol)
--------------   ---------------------
   0.4             0.046
--------------------------------------

Excess chemical potential:

Solute                   KH1       (GF)      KH3    (GF)
-----------------  ---------  ---------  -------  ------
Alcohols
CH~3~OH                4.217     -5.122    2.737  -0.747
C~2~H~5~OH            10.054     -4.981    6.123   1.404
C~3~H~7~OH            18.710     -3.413   11.318   5.485
C~4~H~9~OH            24.899     -4.101   14.739   7.826
C~5~H~11~OH           32.443     -3.875   19.258  11.132
Carboxylate acids
HCOOH                  3.894     -4.934    3.478  -0.043
CH~3~COOH             10.826     -3.775    7.247   2.509
C~2~H~5~COOH          16.855     -4.145   10.943   5.063
C~3~H~7~COOH          24.878     -3.180   15.812   8.811
C~4~H~9~COOH          32.666     -2.990   20.243  12.056
--------------------------------------------------------

2. use the rule introduced in [2] for "any model with embedded sites",
which was first applied to model water when solving RISM equation in
Amber code:

$$
\frac{\sigma_e}{2} = \frac{\sigma_h}{2} - b_{he} $$
$$
\epsilon_e = 0.1\epsilon_h $$

where $\sigma_e$ is the radius of the embedded site, $\sigma_h$ is the
radius of the host site, and $b_{he}$ is the bond length between the
two. Following this rule, we build the new vdw parameters for H in "-OH"
and "-COOH":

Site    $\sigma_h$  $\epsilon_h$  $b_{he}^*$  $\sigma_e$  $\epsilon_e$
------ ----------- ------------- ----------- ----------- -------------
-OH           3.12          0.17      0.9488      1.2224         0.017
-COOH         3.00          0.17      0.9525       1.095         0.017
----------------------------------------------------------------------

^*^ $b_{he}$ is calculated by averaging the O-H bond length relatively
in alcohols and carboxylate acids which we currently use. Note that
number from pronoic acid is excluded because H atom sites there are
derived by ViewerPro automatically with bond length equal to 1.0\AA
(FIXME: even in optimized strucure, O-H bond in -COOH is still 0.98 \AA
long).

Solute                   KH1       (GF)      KH3    (GF)
-----------------  ---------  ---------  -------  ------
Alcohols
CH~3~OH                5.483     -3.818    2.889  -0.607
C~2~H~5~OH            11.327     -3.647    6.277   1.546
C~3~H~7~OH            19.887     -2.334   11.401   5.558
C~4~H~9~OH            26.214     -2.884   14.903   7.978
C~5~H~11~OH           33.779     -2.660   19.422  11.285
Carboxylate acids
HCOOH                  5.789     -2.931    3.758   0.229
CH~3~COOH             12.706     -1.739    7.515   2.769
C~2~H~5~COOH          18.971     -1.858   11.305   5.418
C~3~H~7~COOH          26.804     -1.084   16.079   9.070
C~4~H~9~COOH          34.503     -0.957   20.517  12.321
--------------------------------------------------------

[1] Theoretical Studies of Hydrogen Bonding, Reiher, III., W.E., Ph.D.
    Thesis, Department of Chemistry, Harvard University, Cambridge,
    MA, USA, 1985

[2] Three-Dimensional Molecular Theory of Solvation Coupled with
    Molecular Dynamics in Amber, Tyler Luchko, Sergey Gusarov, Daniel R.
    Roe, Carlos Simmerling, David A. Case , Jack Tuszynski and Andriy
    Kovalenko, J. Chem. Theory Comput., 2010, 6 (3), pp 607–624
    http://dx.doi.org/10.1021/ct900460m


Effects of molecule conformation, coordinates and location
----------------------------------------------------------

In the original code, methanol model used by Jager has the eclipsed
conformation but it has been confirmed both by experimental [1] and
computational [2] results that staggered conformation is preferred in
equilibrium structure.

Here we compared the different chemical potential results obtained from
eclipsed and staggered conformations, the latter of which is built by
rotating the "-OH" bond 180 degrees along "C-O" axis. Here we used
option 1 in last section for the force field parameters of H in hydroxyl
group, the rests are from OPLS-AA.

Conformation   KH1   (GF)    KH3    (GF)
------------ ----- ------ ------  ------
eclipsed     2.412 -6.877  0.995  -2.517
staggered    2.316 -6.975  0.955  -2.558
----------------------------------------

In order to compare the results between the same type of solute but with
different atomic coordinates, we first need to examine the changes
caused by molecule translation, which is actually the consequance of
boundary condition introduced to obtain solute electrostatic field when
solving 3D-RISM equation (no such problem in 1D-RISM). In the following
table, atomic coordinates in methanol(staggered) solute are shifted
different distances along x axis.

shift($\AA$)    KH3    (GF)
------------ ------  ------
no shift      0.955  -2.558
1.0           0.956  -2.558
2.0           0.966  -2.547
5.0           1.118  -2.398
8.0          -2.259  -5.645
---------------------------

It could be seen that results don't change much if solute molecule is
shifted a tiny distance (~ 1.0 \AA). So we could safely translate
methanol molecules to make C atom at the center of simulation box and
compare the values of two different sets of atomic coordinates

solute             KH1   (GF)   KH3   (GF)   $\angle$COH   L~CO~   L~OH~
---------------- ----- ------ ----- ------ ------------- ------- -------
Jager(staggered) 2.317 -6.974 0.955 -2.559 109.5$^\circ$ 1.41\AA 1.08\AA
NIST             4.223 -5.119 2.736 -0.748 107.4$^\circ$ 1.39\AA 0.95\AA
diff             1.906  1.855 1.781  1.811
------------------------------------------------------------------------

We observed another case of value changes when using different molecule
geometry for propanoic acid. Since hydrogen site information is missing in the
original strucure obtained from NIST database, we used ViewerPro to add
hydrogen atom with C-H and O-H bond length equal to 1.0 \AA, while the
determined values in other carboxylate acids are about 1.10 \AA and 0.95 \AA
relatively. The unphysical bond length didn't result in the failure of
molecular mechanic RISM calculation because all the solutes and solvents are
treated as rigid molecules, however this is not true in QM calculation. SCF
calculation in PG diverged (gp) or converged to a unreasonable value (aq) due
to the unpysical C-H and O-H bond length. Here we also compare the results for
both molecule structures, "raw" for the one auto-completed by ViewerPro and
"opt" for that optimized by PG

solute     KH1    (GF)     KH3    (GF)    $\angle$HOC    L~HO~    L~HC~
------  ------  ------  ------  ------  -------------  -------  -------
raw     16.065  -4.682   9.849   4.125  109.6$^\circ$  1.00\AA  1.00\AA
opt     16.855  -4.145  10.943   5.063  107.8$^\circ$  0.98\AA  1.10\AA
-----------------------------------------------------------------------

[1] Serrallach, Alberto. Conformation of methanol. matrix and gas
    studies, and normal coordinate analysis of methyl alcohol and seven
    deuterated species Alberto Serrallach. (1973).
    http://dx.doi.org/10.3929/ethz-a-000086069

[2] Pophristic, V. and L. Goodman, Origin of Staggered Conformational
    Preference in Methanol. The Journal of Physical Chemistry A, 2002.
    106(8): p. 1642-1646. http://dx.doi.org/10.1021/jp014287d


Dielectric effect of water solvent
----------------------------------
The dielectric constant calculated from RISM equation is very low (~ 20)
when being compared to the real value (78.4) at the given condition in
our study. Here we compared the excess chemical potential values of
alcohols, carboxylate acids and aldehyes in aqueous solution with and
without consistant dielectric correction.

Again vdw parameters for H in hydroxyl and carboxyl group are set as
$\sigma$ = 0.4 \AA  and $\epsilon$ = 0.046 kcal/mol. Note that since the
dielectric correction is only available in 1D-RISM code, comparison
calculations of 3D both use the 1D solvent susceptibility, which result
in the slight differences (numbers differ from the second digit after
decimal point) between unmodified 3D-RISM results here with those using
solvent susceptibility calculated from 3D machinery.

Solute               KH1    (GF)  KH1^d^  (GF)^d^     KH3    (GF)   KH3^d^  (GF)^d^
----------------- ------  ------ ------- --------  ------  ------  ------- --------
Alcohols
CH~3~OH            4.217  -5.122   3.661   -5.639   2.771  -0.713    2.411   -1.079
C~2~H~5~OH        10.054  -4.981   9.448   -5.512   6.167   1.447    5.799    1.076
C~3~H~7~OH        18.710  -3.413  18.160   -3.869  11.364   5.530   11.104    5.270
C~4~H~9~OH        24.899  -4.101  24.304   -4.565  14.791   7.877   14.539    7.629
C~5~H~11~OH       32.443  -3.875  31.836   -4.323  19.316  11.189   19.104   10.984
Carboxylate acids
HCOOH              3.894  -4.934   3.409   -5.384   3.510  -0.012    3.216   -0.309
CH~3~COOH         10.826  -3.775  10.341   -4.206   7.286   2.548    7.025    2.285
C~2~H~5~COOH      16.855  -4.145  16.315   -4.604  10.988   5.109   10.738    4.860
C~3~H~7~COOH      24.878  -3.180  24.366   -3.584  15.861   8.859   15.685    8.689
C~4~H~9~COOH      32.666  -2.990  32.153   -3.370  20.298  12.109   20.152   11.973
Aldehydes
CH~2~O             5.961  -1.395   5.665   -1.663   4.174   1.059    3.950    0.833
C~2~H~4~O         13.080   0.175  12.796   -0.006   8.187   3.814    8.010    3.637
C~3~H~6~O         20.547   1.240  20.272    1.041  12.452   6.923   12.328    6.802
C~4~H~8~O         27.879   1.428  27.602    1.251  16.509   9.887   16.425    9.809
C~5~H~10~O        34.337   0.864  34.009    0.676  20.694  12.888   20.628   12.832
-----------------------------------------------------------------------------------

^d^ marks for dielectric correction

3D-RISM coupling QM calculation
-------------------------------

Here we list results of 3D-RISM coupling PG calcuations. "gp" is short
for "gas phase" in which we only add the solvation contribution to total
hamiltonian after SCF is converged, on the other hand, "aq" is short for
"aqueous", in which total hamiltonian is updated with the contribution
of solvation during SCF.

Solute               KH3-gp  (GF)-gp    KH3-aq   (GF)-aq
-----------------  -------- --------  --------  --------
Alcohols
CH~3~OH               2.664   -0.871    -1.757    -5.362
C~2~H~5~OH            7.117    2.356     2.403    -2.431
C~3~H~7~OH           11.945    6.047     7.826     1.863
C~4~H~9~OH           15.717    8.730    11.352     4.298
C~5~H~11~OH          20.104   11.874    15.649     7.352
Carboxylate acids
HCOOH                 0.872   -2.718    -3.063    -6.714
CH~3~COOH             4.139   -0.704    -0.578    -5.497
C~2~H~5~COOH          7.613    1.784     2.569    -3.342
C~3~H~7~COOH         12.734    5.591     7.932     0.710
C~4~H~9~COOH         17.259    8.903    12.489     4.052
Aldehydes
CH~2~O                4.973    1.828     1.029    -2.155
C~2~H~4~O             8.402    3.981     3.517    -0.965
C~3~H~6~O            12.525    6.934     7.868     2.214
C~4~H~8~O            16.640    9.941    12.032     5.270
C~5~H~10~O           21.194   13.287    16.307     8.339
Others
NH~3                  1.001   -1.300    -0.348    -2.672
HCN                   3.346    0.551    -0.883    -3.725
--------------------------------------------------------


PCM results for comparison
--------------------------

COSMO PCM calculations for neutral solutes were performed to help us to
evaluate the results from RISM/QM. Relevant settings in PG:

~~~~~~~~~~~~~~~~~~~~~~~~
((tasks (solvation #t))
...
 (solvation
   (smoothing "fixpva")
   (sol_start_cycle 10))
~~~~~~~~~~~~~~~~~~~~~~~~


Solute                       gp, au           pcm, au  $\Delta$G, kcal/mol
-----------------  ----------------  ---------------- --------------------
Alcohols
CH~3~OH            -23.989683122773  -23.998319182508    -5.41904112311515
C~2~H~5~OH         -30.849679034526  -30.856471978985    -4.26250471857791
C~3~H~7~OH         -37.700008107295  -37.704165508302    -2.60872755788243
C~4~H~9~OH         -44.551459664039  -44.553794563805    -1.46512625416734
C~5~H~11~OH        -51.401870031843  -51.402318967560    -0.28170267306033
Carboxylate acids
HCOOH              -38.791144066439  -38.802950998127    -7.40873156427563
CH~3~COOH          -45.658145206524  -45.668683191517    -6.61248020325757
C~2~H~5~COOH       -52.472298627589  -52.481897083463    -6.02293507637626
C~3~H~7~COOH       -59.359945155032  -59.366946909083    -4.39353064946199
C~4~H~9~COOH       -66.211399426322  -66.216041225270    -2.91268242188052
Aldehydes
CH~2~O             -22.770614800323  -22.776602359488    -3.75713350044585
C~2~H~4~O          -29.655554129757  -29.660873088262    -3.33759327230245
C~3~H~6~O          -36.505028147915  -36.508490300566    -2.17246616697599
C~4~H~8~O          -43.355501463930  -43.357098924263    -1.00239038435417
C~5~H~10~O         -50.211029197881  -50.210801248069     0.14303622753188
Others
H~2~O              -17.154353092917  -17.169639196257    -9.59187698481660
NH~3               -11.679500598278  -11.685839458397    -3.97757133607131
HCN                -15.983500809327  -15.990357517679    -4.30251592379648
--------------------------------------------------------------------------

1D-RISM for ions
----------------

Joung et al. showed comparion of DRISM and MD results for akali halide
solutions [JLC13], with a newly developed force field for different
water models [JC08]. Using the same force field, our 1D-RISM code
produced quite close values. Settings used in the following table:

~~~
--solvent "water, SPC/E" --rho 0.0333295 --beta 1.6889 --L 160 --N 4096 \
--solute "Ion, JC08, SPC/E" --norm-tol 1e-14 --dielectric 78.4
~~~

Ion          KH    KH^*^       HNC   HNC^*^   MD^*^    Exp^*^
-----  -------- --------  -------- --------  -------  -------
Li+     -108.63  -108.30   -111.42  -111.87   -113.3   -113.8
Na+      -85.15   -84.95    -86.24   -86.48    -88.4    -88.7
K+       -68.83   -68.57    -69.35   -69.64    -71.0    -71.2
Cl-      -77.03   -78.09    -79.51   -78.48    -89.3    -89.1
Br-      -71.53   -72.20    -73.72   -73.09    -82.7    -82.7
I-       -62.41   -62.65    -64.46   -64.28    -74.4    -74.3
-------------------------------------------------------------
^*^ values are from [JLC13]

[JC08] Determination of Alkali and Halide Monovalent Ion Parameters for
  Use in Explicitly Solvated Biomolecular Simulations, I. S.  Joung and
  T.  E. Cheatham, The Journal of Physical Chemistry B, 2008, 112 (30),
  pp 9020-9041 http://dx.doi.org/10.1021/jp8001614

[JLC13] Simple electrolyte solutions: Comparison of DRISM and molecular
  dynamics results for alkali halide solutions, In Suk Joung, Tyler
  Luchko, and David A. Case, J. Chem. Phys. 138, 044103 (2013),
  http://dx.doi.org/10.1063/1.4775743


Free energy surface of uranyl
-----------------------------

Gas-phase and PCM energies:

d(UO), au    E(gp), au         E(pcm), au
---------    -----------------  ---------------------
2.8          -28107.5416236806  -28108.0349729507
2.85         -28107.5915625210  -28108.0870568919
2.9          -28107.6316805462  -28108.1293724874
2.95         -28107.6635140115  -28108.1634036738
3.0          -28107.6884378135  -28108.1904915288
3.05         -28107.7074880759  -28108.2116606585
3.1          -28107.7213799436  -28108.2276315498
3.15         -28107.7306623692  -28108.2389684036
3.2          -28107.7357913169  -28108.2461423378
3.25         -28107.7372311851  -28108.2496023818
3.3          -28107.7355376646  -28108.2498799759
3.35         -28107.7312878842  -28108.2475555122
3.4          -28107.7249650674  -28108.2431471236
3.45         -28107.7169052478  -28108.2370303650
3.5          -28107.7073109006  -28108.2294489411
3.55         -28107.6963179679  -28108.2205366155
3.6          -28107.6840838658  -28108.2104315647
3.65         -28107.6707831001  -28108.1950774705
3.7          -28107.6565888808  -28108.2193280230 (a)
3.75         -28107.6416991295  -28108.2061146728 (a)

(a) The two PCM energies at the largest d(UO) separation
diverged. Convergence quality degrades towards the larger distances.

There  are  a  few  water  models  in  use  and  also  several  uranyl
parametrization. This  is the matrix  of excess chemical  potential in
kcal (kJ) as an overview:

                     SPC/E             TIP3P
-------------   ---------------   ---------------
Uranyl, SPC     -322.19 (-1349)   -305.18 (-1278)
Uranyl, TIP3P   -320.15 (-1340)   -303.39 (-1270)
Uranyl, GW96    -312.39 (-1308)   -296.51 (-1241)
Uranyl, PM13    -309.46 (-1296)   -293.80 (-1230)

Here $\rho =  0.0333295$, $\beta = 1.6889$, $L = 160$,  $N = 4096$ and
dielectric constant is 78.4 (see example command line flags below).

To  estimate the effect  of the  finite LJ  potential on  the hydrogen
atoms of a modified SPC/E water  we vary slightly the $\epsilon$(H) and
observe the  excess chemical  potential of  uranyl. The  last row
contains an estimate of  a (dimensionless) slope that reflects  the
sensitivity of the results with respect to the input $\epsilon$(H):

d$\epsilon$, kcal  dG(hnc), kcal     dG(gf), kcal      ddG(rbc), kcal
-----------------  ----------------- ----------------- -----------------
-0.0015            -322.066056282789 -355.601006217798 -15.1127353301548
-0.001             -322.108485948297 -355.664512331336 -15.1241598092300
-5.0e-4            -322.134819602101 -355.711740869483 -15.1354954986111
 0.0               -322.191083707032 -355.788724192364 -15.1467436581335
+5.0e-4            -322.229503240141 -355.847684712192 -15.1579054513841
+0.001             -322.292239631960 -355.930798268651 -15.1689839525631
+0.0015            -322.312339116943 -355.971098261618 -15.1799785845785
Slope              -85.79            -127.06           -22.41

Given the  factual $\epsilon$(H)  = 0.0545  kcal instead of  plain zero
in our modified SPC/E  water an extrapolation to  $\epsilon$(H) = 0
would result in a dG(hnc) by 4.6 kcal lower in  an absolute value. For
GF functional the corresponding number would be 6.9 kcal.

Uranyl in TIP3P water using these command line settings:

~~~
--solvent "water" --rho 0.0333295 --beta 1.6889 --L 160 --N 4096 \
--solute "uranyl, TIP3P" --norm-tol 1e-14 --dielectric 78.4
~~~

The grid stepping for L=160 and N=4096 corresponds to dr = 0.039 A and
dk = 0.0196 A^-1^.  Convergence with respect to grid settings:

N\L         40 A               80 A               160 A              320 A
----  -----------------  -----------------  -----------------  -----------------
1024  -303.385029937429  -303.393707484138  -303.451851082673
2046  -303.385028774009  -303.388402173280  -303.388022515336
4096  -303.385070018802  -303.388927301084  -303.385429099049  -303.364092529682
8192  -303.385050215534  -303.388841814498  -303.384446878473  -303.397269051485
16384                                       -303.387139815038  -303.362823504947

This table uses L=160 and N=4096:

d(UO), au   dG(hnc), kcal      dG(gf), kcal       dG(rbc tpt), kcal
---------   -----------------  -----------------  -----------------
2.8         -288.410456707592  -311.119446512732  -10.2460673451379
2.85        -289.795727148196  -312.729114968202  -10.3229597433474
2.9         -291.191983450945  -314.347637528114  -10.3995859574929
2.95        -292.597884745523  -315.974009275080  -10.4756257092612
3.0         -294.012251278381  -317.607531783730  -10.5508007802468
3.05        -295.433996388366  -319.247516647830  -10.6248666677732
3.1         -296.862104758453  -320.892979208202  -10.6976130741169
3.15        -298.295767580739  -322.542120747843  -10.7688864507462
3.2         -299.734755451194  -324.193234239472  -10.8386590408780
3.25        -301.179556168908  -325.848289129347  -10.9071002456897
3.3         -302.630717726626  -327.512706172906  -10.9745164270380
3.35        -304.087980797894  -329.189333644200  -11.0411737645637
3.4         -305.550109896663  -330.875516165704  -11.1071886795370
3.45        -307.015366510514  -332.565823545769  -11.1725593431499
3.5         -308.482051126470  -334.255417376508  -11.2372555348247
3.55        -309.948803459862  -335.941306473582  -11.3012830477001
3.6         -311.414673815431  -337.622127967608  -11.3647028697802
3.65        -312.879075865794  -339.297562424746  -11.4276205378216
3.7         -314.341701699427  -340.967863737705  -11.4901642515472
3.75        -315.802435722865  -342.633560189514  -11.5524629912869

Now  almost  the  same  with  SPC/E  water  and  corresponding  uranyl
parametrization. Command line settings:

~~~
--solvent "water, SPC/E" --rho 0.0333295 --beta 1.6889 --L 160 --N 4096 \
--solute "uranyl, SPC" --norm-tol 1e-14 --dielectric 78.4
~~~

d(UO), au   dG(hnc), kcal      dG(gf), kcal       dG(rbc tpt), kcal
---------   -----------------  -----------------  -----------------
2.8         -304.864370033662  -335.189869413676  -14.3255579935222
2.85        -306.485402711989  -337.151379829175  -14.4099986645670
2.9         -308.114600499118  -339.111627262736  -14.4936426860625
2.95        -309.750585559458  -341.069999583972  -14.5761550920681
3.0         -311.392233376546  -343.026815765087  -14.6572731013174
3.05        -313.038571601192  -344.982749232697  -14.7367862975397
3.1         -314.688727441137  -346.938286437927  -14.8145316296692
3.15        -316.342044362600  -348.892850964879  -14.8904254167647
3.2         -317.998517724657  -350.845322563302  -14.9645765168533
3.25        -319.659111239045  -352.798961561739  -15.0374010438737
3.3         -321.325123353192  -354.763297198121  -15.1094961063650
3.35        -322.996998829052  -356.746294172095  -15.1813133624776
3.4         -324.673817545879  -358.747645935985  -15.2529810217050
3.45        -326.353716220755  -360.760847461796  -15.3243883375494
3.5         -328.034597349590  -362.778475151548  -15.3953598276761
3.55        -329.714614589889  -364.795107482580  -15.4657703302418
3.6         -331.392353150121  -366.807648151674  -15.5355783955610
3.65        -333.066820378392  -368.814634488579  -15.6048119030658
3.7         -334.737361226935  -370.815516986285  -15.6735381407419
3.75        -336.403560961928  -372.810177067083  -15.7418356723166


A linear regression  between 1.48-1.91 A (for larger  distances PCM is
unstable) gives  the following  slope and intercept.   For convenience
intercept is computed at 1.76 A:

              Intercept, kcal    Slope, kcal/A
----------   ----------------   --------------
PCM             -323.32            -48.48
HNC/SPC         -322.23            -62.71
HNC/T3P         -303.43            -54.42
(GF)/SPC        -355.87            -74.54
(GF)/T3P        -328.44            -62.65
(RBC)/SPC        -15.14             -2.84
(RBC)/T3P        -11.00             -2.64

The free energy  profile along the assymetric stretching  mode is only
slightly convex upward. Here we move  the U along the axis keeping the
uranyl oxygens in place. Again this is uranyl (SPC) in water (SPC/E):

z(U), au    dG(hnc), kcal     dG(gf), kcal      ddG(rbc), kcal
--------  ----------------- ----------------- -----------------
-0.5      -322.925715135858 -356.845648757879 -15.1821632906528
-0.45     -322.793975753692 -356.644958828025 -15.1759846238581
-0.4      -322.673057828675 -356.463816279959 -15.1702584270575
-0.35     -322.563913109193 -356.303290159633 -15.1650278812779
-0.3      -322.467426009134 -356.164249231148 -15.1603487855479
-0.25     -322.384396286184 -356.047278182368 -15.1562821276392
-0.2      -322.315521626981 -355.952586206166 -15.1528862378043
-0.15     -322.261380882580 -355.879976022087 -15.1502096110067
-0.1      -322.222417442073 -355.828919163098 -15.1482852934068
-0.05     -322.198930097758 -355.798710552585 -15.1471289447639
 0.0      -322.191083707032 -355.788724192364 -15.1467436581335
 0.05     -322.198930108592 -355.798710563459 -15.1471289447915
 0.1      -322.222417440885 -355.828919161898 -15.1482852934106
 0.15     -322.261380882610 -355.879976022125 -15.1502096110092
 0.2      -322.315521625962 -355.952586205147 -15.1528862378043
 0.25     -322.384396285832 -356.047278182012 -15.1562821276380
 0.3      -322.467426012274 -356.164249234285 -15.1603487855473
 0.35     -322.563913110512 -356.303290160980 -15.1650278812795
 0.4      -322.673057828630 -356.463816279918 -15.1702584270580
 0.45     -322.793975751603 -356.644958825941 -15.1759846238602


First peak and number integral of uranyl
----------------------------------------

The numbers in the following table need to be used carefllly, because
the location and values for first maximum and minimum and obtained from
the linear interpolated function of RDF. Values for number integral are
also from interpolation. Calulation settings are the same as those for
free energy except L=16, since interpolation in a small range [0.0, 3.0]
can not be accurate for L=160, especially it's hard to locate the peaks.
The energy difference is less than 0.1 kcal/mol even L is shortened to
1/10 (see last section).

SPC/E water:

                r~max~   g(r~max~)   r~min~   g(r~min~)   ni(r~min~)
-------------  -------  ----------  -------  ----------  -----------
Uranyl, GW96    2.4902     12.2846   3.2402      0.0136       6.2542
Uranyl, SPC     2.4238     13.1930   3.1699      0.0102       6.1388
Uranyl, TIP3P   2.4316     12.9512   3.1816      0.0110       6.1241
Uranyl, PM13    2.4238     12.7383   3.1778      0.0069       5.7993
--------------------------------------------------------------------


TIP3P water:

                r~max~   g(r~max~)   r~min~   g(r~min~)   ni(r~min~)
-------------  -------  ----------  -------  ----------  -----------
Uranyl, GW96    2.4746     10.4921   3.2011      0.0124       5.2053
Uranyl, SPC     2.4082     11.2197   3.1348      0.0092       5.0920
Uranyl, TIP3P   2.4199     11.0207   3.1426      0.0099       5.0823
Uranyl, PM13    2.4121     10.8627   3.1387      0.0066       4.8429
--------------------------------------------------------------------

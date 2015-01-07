                           4(a)               5             5 + 1
-----------  ------    --------------  ---------------  -------------
d(U-OW), pm  MM         239.2*4        243.5*5            243.7*3
                                                          242.3*2
                                                          (243.1)
                                                          410.6
             MM+KH      241.0*2        239.9*5            239.9*3
                        241.4*2                           240.0*2
                                                          (239.9)
                                                          426.3

             MM+PSE2    241.7*2        240.0*5            239.9*3
                        241.6*2                           240.0*2
                                                          427.0

             MM+PSE3    NA(c)          240.0*5            239.9*3
                                                          240.0*2
                                                          427.4

E(MM/MM)               -232.5         -271.6             -295.2
E(KH/MM)               -154.2         -144.8             -130.7
E(PSE2/MM)             -159.6         -150.1             -137.6
E(PSE3/MM)             -160.8         -151.1             -138.8
E(MM+KH/MM)            -386.7         -416.4             -425.9
E(MM+PSE2/MM)          -392.0         -421.7             -432.4
E(MM+PSE3/MM)          -393.3         -422.7             -433.5+

E(MM/MM+KH)            -219.1         -270.7             -293.8
E(MM/MM+PSE2)          -217.6         -271.6             -294.8
E(MM/MM+PSE3)          NA             -270.7             -293.8
E(KH/MM+KH)            -183.3         -146.7             -133.5
E(PSE2/MM+PSE2)        -195.6         -151.8             -139.5
E(PSE3/MM+PSE3)        NA             -152.8             -140.7
E(MM+KH/MM+KH)         -402.4         -417.4             -427.4
E(MM+PSE2/MM+PSE2)     -413.2         -422.5+            -433.3
E(MM+PSE3/MM+PSE3)     NA             -423.6             -434.5+

dG           MM+KH     -382.0         -391.9             -396.8
             MM+PSE2   -389.1         -392.5-            -397.2
             MM+PSE3   NA             -392.3             -397.0


                           4(b)
                      --------------
             MM+KH      235.8*4
             MM+PSE2    235.9*4
             MM+PSE3    236.0*4

E(MM/MM+KH)            -231.8
E(MM/MM+PSE2)          -232.5-
E(MM/MM+PSE3)          -232.5-
E(KH/MM+KH)            -155.6
E(PSE2/MM+PSE2)        -160.8
E(PSE3/MM+PSE3)        -162.1
E(MM+KH/MM+KH)         -387.4
E(MM+PSE2/MM+PSE2)     -392.7
E(MM+PSE3/MM+PSE3)     -393.9

dG           MM+KH     -367.0
             MM+PSE2   -368.6
             MM+PSE3   -368.9

                             0
                     ----------------
dG           KH        -363.1
             PSE2      -378.7
             PSE3      -383.4

E(KH/water)            -5.10197173104614
E(PSE2/water)          -6.02261426875449
E(PSE3/water)          -6.2484516865421
E(HNC/water)           -6.34905163542536


a) The water planes in the under-coordinated complex with 4 waters
   with MM+RISM are tilted by 26 and 50 degrees, those water close to
   the missing explicit ligand are tilted more. This non-symmetric
   structure is obtained by starting from the same non-symmetric
   geometry as in MM input. In a re-optimization using PSE2 closure
   and starting from KH geometries the tilting angles became close to
   20 and 77 degrees. The waters close to the missing ligand lie
   almost in the equatorial plane.

b) D4H stationary structure obtained in MM+KH optimization starting
   from the D4H MM geometry.

c) PSE3 closure with under-coordinated uranyl in the non-symmetric
   conformation diverges.



QM water in SPC/E geometry and  uranyl with 172 pm, C1, ADKH, BP.  See
1w,c1.scm and uranyl,c1.scm.

         Water                 Uranyl
30/131   -76.502720816987      -28110.765076538421
50/291   -76.502943687787      -28110.790576816220
50/291   -76.502960020395      -28110.790792414213 (older)


Therefore, to make the energies comparable to MM one needs to subtract
the energy of GP uranyl + 6 waters. Make sure to use

        ase.units.Hartree  = 27.211395655517308 eV/H

(used in ASE PG calculator) when converting Hartrees to eV.

20/91    -28110.931719705168 -76.504226371401 * 6 =
         -28569.957077933574 H =
         -777428.40590879761927829263 eV
30/131   -28110.765076538421 -76.502720816987 * 6 =
         -28569.781401440343 H =
         -777423.62550623293672912596 eV
50/291   -28110.790576816220 -76.502943687787 * 6 =
         -28569.808238942942 H =
         -777424.35579213456409218810 eV

Two  "minima"   obtained  by  unconstrained   minimization  with  bias
potential applied to 5 waters:

Grid 30/131/64:

Struct    G, eV              dG, kcal
5+1       -777442.672597     -439.24 (-0.20, bent uranyl)
5+1       -777442.663931     -439.04 (0.0, linear uranyl)
(340pm)   -777442.696429             (-0.75-)
(320pm)   -777442.293271             (+8.55-)
(310pm)   -777442.23707              (+9.84)
(300pm)   -777442.250265             (+9.54)
6         -777442.445809     -434.01 (+5.03)

Grid 50/291/96:
5+1       -777444.015779     -453.37   (uranyl bent)
(320pm)   -777442.837054     -426.19
(310pm)   -777442.795784     -425.24
(300pm)   -777442.784967     -424.99
6         -777442.947949     -428.75-

The  5+1 model  is apparently  affected by  the bias,  as the  two U-O
distances  are longer  than  250 pm  (unless  uranyl is  forced to  be
linear). Quadratic  bias potential has the  effect of 2.5  meV or 0.06
kcal for an elongation by 1  pm beyond 250 pm.  In the six-coordinated
complex  the 6th  water  with U-O  of  253 pm  has  no bias  potential
applied.  FIXME: with 50/291/96 grid the 5+1 structure features uranyl
bent by 137 deg (reoptimize with linear?).

Distances for 5+1:

30/131/64     50/291/96
0.179         0.181
0.180         0.184
0.395         0.383
0.245         0.247
0.245         0.251
0.242         0.251
0.251         0.252
0.251         0.252

Linear Uranyl:

30/131/64     50/291/96
0.179
0.180
0.401
0.245
0.244
0.242
0.250
0.250

Distances for 6 (nm):

30/131/64     50/291/96
0.180         0.180
0.180         0.180
0.253         0.254
0.250         0.250
0.250         0.250
0.250         0.250
0.250         0.250
0.250         0.250

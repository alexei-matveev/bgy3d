
Here   uranyl-aqua  complexes  have   been  optimized   with  symmetry
constraints and  rigid water molecules in SPC  geometry. Gas-phase and
PCM  geometries thus  differ.  Assuming  Hartree  is 627.5154624673408
kcal here.

NOTE: geomeries in *.gx files and in the input are taken from the last
prediction during  geometry optimization.  The energies,  on the other
hand,  correspond to  the  previous geometry  optimization step.   The
second  number  (where  available)  is  obtained  in  a  single  point
calculation using that  predicted geometry. The difference illustrates
the effect of such a geometry mismatch on the energies.

n    sym   method    E(tot), au            dE, kcal  dG, kcal  U-OU, pm  U-OW, pm
---  ---   -------   -------------------  ---------  --------  --------  --------
0    D8H   gp/gp     -28110.790792414213                            172
           pcm/pcm   -28111.304172875658    -322.15   -322.15       174
                     -28111.304172891396
           gp/pcm    -28110.790109421156     [0.43]
                                             (0.43)
           rs/rs     -28111.394094040123    -378.58                 180
           gp/rs     -28110.780086885854
4    D4H   gp/gp     -28417.190831100434                            177       242
           pcm/pcm   -28417.510911499532    -200.85   -416.23       179       237
                     -28417.510911501795
           gp/pcm    -28417.189028477842  [-242.47]
                                             (1.13)
           rs/rs     -28417.471645579732    -176.22   -388.34       181       246
           gp/rs     -28417.187306573804
5    D5H   gp/gp     -28493.736226891073                            177       249
           pcm/pcm   -28494.041471559878    -191.55   -426.48       179       243
                     -28494.041471561650
           gp/pcm    -28493.734208162346  [-268.96]
                                             (1.27)
6    D3D   gp/gp     -28570.268746168604                            178       253
           pcm/pcm   -28570.554049949551    -179.03   -425.45       180       249
                     -28570.554049948907
           gp/pcm    -28570.265528340671  [-286.76]
                                             (2.02)
wat  C2V   gp/spc       -76.502960020395
           pcm/spc      -76.514210290228      -7.06


1) dE = PCM(mol/pcm) - GP(mol/gp), except for water, where both
   geometries are SPC.

2) dG = PCM(complex/pcm) - n * PCM(water/spc) - GP(uranyl/gp)

3) In square brackets, the "internal" energy dE = GP(complex/pcm) - n
   * GP(water/spc) - GP(uranyl/gp).

4) In round parens, dE  = GP(complex/pcm) - GP(complex/gp).  The value
   is always positive since GP/GP is supposed to be a minimum. This is
   also a  measure of  internal energy difference  between PCM  and GP
   geometries of the complex.

5) The GP/RISM row abbreviated as gp/rs in the table is the SCF energy
   at the  end of combined  QM+RISM geometry optimization.  It happens
   not to contain the RISM term which is added afte SCF.


Using optimized PCM  geometries one could estimate the  dG by the RISM
method. Using The following settings

~~~
runbgy.scm energy --norm-tol 1e-14 --dielectric 78.4 --rho 0.0333295 --beta 1.6889
 --L 160 --N 4096 --solvent "water, PR-SPC/E" --solute "uranyl, 0w, pcm"
~~~

one obtains the two major contributions:

 n    E(self), kcal      E(RISM), kcal      dG, kcal
 ---  -----------------  -----------------  --------
 0       0.0             -376.027350192114   -376.03
 4    -232.074402281419  -174.594711875714   -375.19
 5    -271.654767148776  -159.157414678584   -391.47
 6    -282.556908726866  -143.250721366546   -378.59

The  free solvation  energy  of  the complex  F(complex)  = E(self)  +
E(RISM).    Using  the   free  energy   of  the   water,   F(water)  =
-7.86929437074687 kcal,  one obtains four different  estimates for the
uranyl solvation energy, dG = F(complex) - n * F(water).

After a  moment of doubts, I  verified if analytical  gradients of the
RISM  term are  correct.  Here  is  the comparison  of analytical  and
numerical gradients for  a breathing mode of D4H  structure. Yep, they
are  correct. Also note  how repulsive  the RISM  contribution becomes
beyond ~2.25 A for U-OW distance:

d, A  E(RISM), kcal     g(d), kcal/A         g(d), num, kcal/A
----  ----------------- -----------------    -----------------
1.0   -223.989586837144  51.1111337020776    51.1109874627467
2.0   -175.174420287068  24.784343435638     24.7841831723102
2.25  -172.007010238754  -0.861532764046105  -0.861642707694627
2.5   -176.322534336077 -34.3303279112493   -34.3303526218852
2.75  -188.80506314766  -64.6490259071153   -64.6490680168614
3.0   -209.046655343837 -97.7056202265256   -97.7056442888162
4.0   -313.802313300318 -76.6694172030191   -76.6694266021152
5.0   -359.502910150843 -22.9967182245357   -22.9966927484558

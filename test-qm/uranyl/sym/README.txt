
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
4    D4H   gp/gp     -28417.189498272211                            177       247
           pcm/pcm   -28417.508816589663    -200.38   -414.90       178       243
                     -28417.508816445374
           gp/pcm    -28417.190504525548  [-243.40]
                                            (-0.63)
5    D5H   gp/gp     -28493.736050294861                            177       250
           pcm/pcm   -28494.041286317046    -191.54   -426.36       179       245
                     -28494.041286400123
           gp/pcm    -28493.735197376816  [-269.58]
                                             (0.54)
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

4) In round parens, dE = GP(complex/pcm) - GP(complex/gp).  For some
   reason in one case the value is negative even though GP/GP is
   supposed to be a minimum. This is also a measure of internal energy
   difference between PCM and GP geometries of the complex.


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
 4    -231.098378481216  -175.547503782298   -375.17
 5    -271.373959858618  -158.999570912819   -391.03
 6    -282.556908726866  -143.250721366546   -378.59

The  free solvation  energy  of  the complex  F(complex)  = E(self)  +
E(RISM).    Using  the   free  energy   of  the   water,   F(water)  =
-7.86929437074687 kcal,  one obtains four different  estimates for the
uranyl solvation energy, dG = F(complex) - n * F(water).



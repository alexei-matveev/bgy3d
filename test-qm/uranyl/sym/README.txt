
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

n    sym   method    E(tot), au           dE, kcal  dG, kcal  U-OU, pm  U-OW, pm
---  ---   -------   -------------------  --------  --------  --------  --------
0    D8H   gp/gp     -28110.790792414213                           172
           pcm/pcm   -28111.304172875658   -322.15   -322.15       174
                     -28111.304172891396
           gp/pcm    -28110.790109421156    (0.43)
4    D4H   gp/gp     -28417.189498272211                           177       247
           pcm/pcm   -28417.508816589663   -200.38   -414.90       178       243
                     -28417.508816445374
           gp/pcm    -28417.190504525548   (-0.63)
5    D5H   gp/gp     -28493.736050294861                           177       250
           pcm/pcm   -28494.041286317046   -191.54   -426.36       179       245
                     -28494.041286400123
           gp/pcm    -28493.735197376816    (0.54)
6    D3D   gp/gp     -28570.268746168604                           178       253
           pcm/pcm   -28570.554049949551   -179.03   -425.45       180       249
                     -28570.554049948907
           gp/pcm    -28570.265528340671    (2.02)
wat  C2V   gp/spc       -76.502960020395
           pcm/spc      -76.514210290228     -7.06


*) dE = PCM(mol) - GP(mol)
   dG = PCM(complex) - n * PCM(water) - GP(uranyl)

**) In parens, dE = GP(complex/pcm) - GP(complex/gp). For some reason
    in one case the value is negative even though GP/GP is supposed to
    be a minimum.

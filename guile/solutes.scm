;;; vim: set tw=72:
;;;
;;; Copyright (c) 2007 Lukas Jager
;;; Copyright (c) 2013 Alexei Matveev
;;; Copyright (c) 2013 Bo Li
;;;
;;; Number density is measured in A^-3, you would need a molar mass to
;;; convert that to/from g/l:
;;;
;;;      -3     30  -3
;;;   1 A   = 10   m   = 1660.5 mol/l
;;;
;;; E.g. the molar mass of water  is 18.0153 g/mol, that of the HCl is
;;; 36.46094 g/mol.
;;;
;;; Simplest  species,  single  neutral   site  LJ  particle.   To  be
;;; comparable  to  water  one  might   set  β  =  0.261610  and  ρ  =
;;; 1.054796. These  numbers have been  derived from LJ  parameters of
;;; oxygen  center  of the  two-site  water  model  (see below),  room
;;; temperature  and water density  (see comments  on the  TIP3P water
;;; model).
;;;
("LJ" (("LJ" (0.0 0.0 0.0) 1.0 1.0 0.0)))


;;;
;;;
;;; Table I. Structural data and force field parameters [1].
;;;
;;; Atom        x (A)   y (A)   z (A)   ε (kJ/mol)   σ (A)    q (e)
;;; ----------------------------------------------------------------
;;; K+                                  0.3640       3.14265   1
;;; Cl-                                 0.0422*      4.04468  -1
;;; O(water)    0.000   0.000   0.000   0.6359       3.15080  -0.834
;;; H(water) +/-0.757   0.000   0.586   0.1921       0.40000   0.417
;;;
;;; *) ε probably wrong, see comments below.
;;;
;;; Force field for water is a TIP3P parametrization according to Ref.
;;; [1] with 0.4  A hydrogens (see "water" entry  below).  In Ref. [1]
;;; the  authors used  the "logarithmically  spaced grid  ranging from
;;; 0.0059 to  164.02 A",  T = 298.15  K, ρ  = 0.0333295 A^-3  and the
;;; dielectric constant 78.4 for DRISM equations.
;;;
;;; Note that  the Cl- parameters appear  wrong as quoted  in Ref. [1]
;;; and the table above. The Ref. [1] redirects to [2] which says
;;;
;;;   "Lennard-Jones parameters of the ions were ε = 0.087/0.15 kcal
;;;    mol-1 and Rmin/2 = 1.76375/2.27 A for K+/Cl- ..."
;;;
;;; From which I derive the following table:
;;;
;;;        σ (A)         ε (kJ/mol)
;;; ----------------------------------
;;; K+    3.14265      0.3643  0.3640*
;;; Cl-   4.04468      0.6280  0.6276*
;;;
;;; *) Using thermochemical kcal of 4.184 kJ.
;;;
;;; [1] Treatment of charged solutes in three-dimensional integral
;;;     equation theory, Thomas Kloss and Stefan M. Kast
;;;     J. Chem. Phys. 128, 134505 (2008)
;;;     http://dx.doi.org/10.1063/1.2841967
;;;
;;; [2] Molecular Dynamics Simulation of the Cytosolic Mouth in
;;;     Kcv-Type Potassium Channels, Sascha Tayefeh, Thomas Kloss,
;;;     Gerhard Thiel, Brigitte Hertel, Anna Moroni, and Stefan
;;;     M. Kast, Biochemistry, 2007, 46 (16), pp 4826–4839,
;;;     http://dx.doi.org/10.1021/bi602468r
;;;
;; ("K+" (("K+" (0.0 0.0 0.0) (A 3.14265) (kJ 0.3640) 1.0)))
;; ("Cl-" (("Cl-" (0.0 0.0 0.0) (A 4.04468) (kJ 0.0422) -1.0))) ; ε wrong?
("K+, Kloss & Kast" (("K+" (0.0 0.0 0.0) (A 3.14265) (kcal 0.087) 1.0)))
("Cl-, Kloss & Kast" (("Cl-" (0.0 0.0 0.0) (A 4.04468) (kcal 0.15) -1.0)))

;;;
;;; Almost the same as "water":
;;;
("water, Kloss & Kast"
 (("OW" (0.000 0.000 0.000) (A 3.15080) (kJ 0.6359) -0.834)
  ("HW" (+0.757 0.000 0.586) (A 0.40000) (kJ 0.1921) 0.417)
  ("HW" (-0.757 0.000 0.586) (A 0.40000) (kJ 0.1921) 0.417)))


;;
;; Application of an extended RISM equation to dipolar and quadrupolar
;; fluids, Fumio  Hirata, B. Montgomery Pettitt, and  Peter J. Rossky,
;; J. Chem. Phys. 77, 509 (1982), http://dx.doi.org/10.1063/1.443606
;;
;; TABLE 1. Models for dipolar diatomic fluids.
;;
;; Model       1^a        2^b      3^b
;; -------------------------------------
;; σ++ (A)     3.341     2.735    0.4
;; σ-- (A)     3.341     3.353    3.353
;; ε++ (K)    44.0      20.0     20.0
;; ε-- (K)    44.0     259.0    259.0
;; L+- (A)     1.1       1.257    1.3
;; Z+- (e)^c   0.2       0.2      0.2
;; T   (K)^d  72.0     210.0     210.0
;; ρ (A^-3)^d  0.01867   0.018     0.018
;;
;; a)  N2-like  model.   b)  HCl-like  model;  H  corresponds  to  the
;; positively  charged site,  and Cl  to the  negatively  charged.  c)
;; Magnitude of site charge, if not zero, in units of the magnitude of
;; the electronic  charge.  d) Thermodynamic conditions  used here for
;; fluid.
;;
;; Density: 0.01867
;; T =  72 K => beta = 6.9938
;; T = 125 K => beta = 4.0284
;;
("dipolar nitrogen"
 (("N+" (0.0 0.0 +0.55) 3.341 0.08738 +0.2)   ; ε = 44 K
  ("N-" (0.0 0.0 -0.55) 3.341 0.08738 -0.2))) ; ε = 44 K

;;
;; Application of an extended RISM equation to dipolar and quadrupolar
;; fluids, Fumio  Hirata, B. Montgomery Pettitt, and  Peter J. Rossky,
;; J. Chem. Phys. 77, 509 (1982), http://dx.doi.org/10.1063/1.443606
;;
;; Density: 0.018
;; T= 210 K => beta = 2.39788
;; T= 420 K => beta = 1.1989
;; T= 315 K => beta = 1.5985
;;
("hydrogen chloride"
 (("H" (0.6285 0.0 0.0) 2.735 0.03971 0.2)      ; ε = 20 K
  ("Cl" (-0.6285 0.0 0.0) 3.353 0.51434 -0.2))) ; ε = 259 K

("hydrogen chloride z-aligned"
 (("H" (0.0 0.0 0.6285) 2.735 0.03971 0.2)
  ("Cl" (0.0 0.0 -0.6285) 3.353 0.51434 -0.2)))

;;
;; Parametrisierung CS2 nach Zhu:  H=S O=C
;; sigma_H = 3.520   epsilon_H = 0.39500  q_H = 0.154
;; sigma_O = 3.200  epsilon_O = 0.10128 q_O =-0.308
;; r_OH = 1.56
;; r_HH = 2*1.56
;; theta_HOH = 180.0
;; mass CS2 =  75.999432 u
;; density = 1.263 g cm^-3 =  0.76059653 u/A^3 => 0.010007924 / A^3
;; temperature : T= 298.15 K (25 C) => 0.5921, => beta =1.6889
;; temperature : T= 360 K (? C) => 0.7149, => beta =1.3988
;;
("carbon disulfide"
 (("C" (0.0 0.0 0.0) 3.2 0.10128 -0.308)
  ("S1" (-1.56 0.0 0.0) 3.52 0.395 0.154)
  ("S2" (1.56 0.0 0.0) 3.52 0.395 0.154)))

;;
;; Parametrisierung CS2 nach Tildesley:  H=S O=C
;; sigma_H = 3.520   epsilon_H = 0.36341  q_H = 0
;; sigma_O = 3.35  epsilon_O = 0.10168 q_O =0
;; r_OH = 1.57
;; r_HH = 2*1.57
;; theta_HOH = 180.0
;; mass CS2 =  75.999432 u
;; density = 1.263 g cm^-3 =  0.76059653 u/A^3 => 0.010007924 / A^3
;; temperature : T= 298.15 K (25 C) => 0.5921, => beta =1.6889
;;
("carbon disulfide, Tildesley"
 (("C" (0.0 0.0 0.0) 3.35 0.10168 0.0)
  ("S1" (-1.57 0.0 0.0) 3.52 0.36341 0.0)
  ("S2" (+1.57 0.0 0.0) 3.52 0.36341 0.0)))

;;
;; Parametrisierung TIP3P:
;; sigma_H = 0.400   epsilon_H = 0.046  q_H = 0.417
;; sigma_O = 3.1506  epsilon_O = 0.1521 q_O =-0.834
;; r_OH = 0.9572
;; r_HH = 1.5139
;; theta_HOH = 104.52
;; mass H2O = 18.0154 u
;; density = 1 kg/l = 0.6022142 u/A^3 => 0.033427745 / A^3
;; temperature: T = 298.15 K (25 C) => 0.5921, => beta = 1.6889
;;
;; Here are two  sources for the non-zero sigma  and epsilon parameter
;; for hydrogen:
;;
;; I.  The  one  introduced  by  Pettitt  and  Rossky  [1]  to  "avoid
;; catastrophic  overlap  of the  corresponding  site  charges in  the
;; calculations".   But  it  be  noted  that in  [1],  these  non-zero
;; parameter only contribute to the C_12 term of O-H pair and "present
;; in  a 6-12 potential  with a  depth of  0.2 kcal/mol  and effective
;; diameters of 2.8 Å for oxygen  and 0.4 Å for hydrogen", but in this
;; model ε(H) is quite different from what we used here:
;;
;; In [1], C_12(OH) = 225.180 kcal Å^12 / mol, and LJ potential is
;;
;;   u = C_12/r^12 + C_6/r^6
;;
;; considering the other form of LJ potential
;;
;;   u = 4ε[(σ/r)^12 - (σ/r)^6]
;;
;; since C_6 is zero:
;;
;;   C_12(OH) = 4ε(OH)σ(OH)^12
;;
;; with
;;
;;   σ(OH) = (2.8 + 0.4) / 2 = 1.6
;;
;; so
;;
;;   ε(OH) = C_12(OH) / σ(OH)^12 / 4 = 225.180 / 281.47497671065616 / 4
;;         = 0.2 kcal / mol
;;
;; now we have
;;
;;   ε(H) = ε(OH)^2 / ε(O) = 0.04 / 0.1521 = 0.2630
;;
;; II. By  checking some papers  using this modified TIP3P,  found the
;; paper  introducing the  protein  parameters used  in CHARMM22  [2],
;; which is hopefully the original  source of modified TIP3P. From the
;; discussion in CHARMM official  forum [3], some people metioned that
;; this  modified TIP3P, naming  CHARMM TIP3P,  is only  for practical
;; purposes (to  avoid electrostatic catastrophe).  Further details of
;; the refinement of  CHARMM TIP3P could be found  in the Ph.D. thesis
;; of Reiner  [4] (check the path of  local copy on wiki  page). In my
;; opinion we should cite this or  the CHARMM paper instead of [1] for
;; future usage.  (in the PhD  thesis, ε(H) =  0.04598 kcal /  mol and
;; σ(H) = 0.4490 Å)
;;
;; Reference:
;;
;; [1] Integral equation predictions of liquid state structure for
;;     waterlike intermolecular potentials, B.  M.  Pettitt and P.  J.
;;     Rossky, The Journal of Chemical Physics, 1982, 77 (3),
;;     1451-1457, http://dx.doi.org/10.1063/1.443972
;;
;; [2] All-Atom Empirical Potential for Molecular Modeling and
;;     Dynamics Studies of Proteins, A. D. MacKerell, D. Bashford,
;;     Bellott, R. L.  Dunbrack, J. D. Evanseck, M. J. Field,
;;     S. Fischer, J. Gao, H. Guo , S. Ha, D. Joseph-McCarthy,
;;     L. Kuchnir, K. Kuczera, F. T. K. Lau , C. Mattos, S. Michnick,
;;     T. Ngo, D. T. Nguyen, B. Prodhom, W. E.  Reiher, B. Roux,
;;     M. Schlenkrich, J. C. Smith, R. Stote, J. Straub , M. Watanabe,
;;     J. Wiórkiewicz-Kuczera, D. Yin and M. Karplus, The Journal of
;;     Physical Chemistry B 1998, 102 (18), 3586-3616.
;;     http://dx.doi.org/10.1021/jp973084f
;;
;; [3] http://goo.gl/ldmC8
;;
;; [4] Theoretical Studies of Hydrogen Bonding, Reiher, III., W.E.,
;;     Ph.D. Thesis, Department of Chemistry, Harvard University,
;;     Cambridge, MA, USA, 1985
;;
;;
("water"
 (("O" (-0.2929 0.0 0.0) 3.1506 0.1521 -0.834)
  ("OH" (0.2929 0.757 0.0) 0.4 0.046 0.417)
  ("OH" (0.2929 -0.757 0.0) 0.4 0.046 0.417)))

;;
;; SPC/E (extended simple point charge model):
;;
;;   r(OH) = 0.1 nm
;;   θ(HOH) = 109.47 deg [2 atan (√2), thetrahedral]
;;   q(O) = -0.8476 e
;;   q(H) = 0.4238 e
;;   A = 0.37122 (kJ/mol)^(1/6) * nm
;;   B = 0.3428 (kJ/mol)^(1/12) * nm
;;   d = 2.35 D [dipole moment]
;;
;; Here  A and  B  parametrize LJ  as  -(A/r)^6 +  (B/r)^12, cf.   the
;; dimensions.  Hydrogen bond of a pair is -30.0 kJ/mol which slightly
;; more than -27.6 kJ/mol of the original SPC model [BGS97].
;;
;; * AB-form to σε-form:
;;
;;   σ = B^2/A = 0.316555789019988 nm
;;   ε = (A/B)^12 / 4 = 0.650169580818749 kJ/mol
;;                    = 0.155290336490577 kcal/mol
;;
;;   (assuming IT calorie 4.1868 J)
;;
;; For use with RISM "[t]he original SPC solvent potential is modified
;; to include core repulsion  for the interactions associated with the
;; hydrogen site in  a way that does not alter  the physical nature of
;; intermolecular interactions.   The van der Waals parameters  of σ =
;; 1.0 A  and ε  = 0.0545 kcal/mol  are assigned  to the H  site". See
;; footnote in Ref. [2].
;;
;; In this form thus modified SPC/E potential was used in Ref. [3].
;;
;; References:
;;
;; [BGS97] "The Missing Term in Effective Pair Potentials", H.  J.  C.
;;   Berendsen, J.  R.  Grigera, and T.  P.  Straatsma, J.  Phys.
;;   Chem 1987, 91, 6269-6271. http://dx.doi.org/10.1021/j100308a038
;;
;; [2] "Theoretical study for the basicities of methylamines in
;;     aqueous solution: A RISM-SCF calculation of solvation
;;     thermodynamics", Masaaki Kawata, Seiichiro Ten-no, Shigeki
;;     Kato, Fumio Hirata, Chemical Physics, 203, 1996, 53–67,
;;     http://dx.doi.org/10.1016/0301-0104(95)00352-5
;;
;; [3] "Comparative Study on Solvation Free Energy Expressions in
;;     Reference Interaction Site Model Integral Equation Theory",
;;     Kazuto Sato, Hiroshi Chuman, and Seiichiro Ten-no,
;;     J. Phys. Chem. B, 2005, 109 (36), pp 17290–17295,
;;     http://dx.doi.org/10.1021/jp053259i
;;
("water, SPC/E"
 (("O" (-0.288675134594813  0.000000000000000 0.0) 3.1656 0.1553 -0.8476)
  ("OH" (0.288675134594813  0.816496580927726 0.0) 1.0    0.0545  0.4238)
  ("OH" (0.288675134594813 -0.816496580927726 0.0) 1.0    0.0545  0.4238)))

("water, cSPC/E"
 (("O" (-0.288675134594813  0.000000000000000 0.0) 3.1656 0.1553 -0.8476)
  ("OH" (0.288675134594813  0.816496580927726 0.0) 1.1658 0.01553 0.4238)
  ("OH" (0.288675134594813 -0.816496580927726 0.0) 1.1658 0.01553 0.4238)))

("water, PR-SPC/E"
 (("O" (-0.288675134594813  0.000000000000000 0.0) 3.1656 0.1553 -0.8476)
  ("OH" (0.288675134594813  0.816496580927726 0.0) 0.4    0.046  0.4238)
  ("OH" (0.288675134594813 -0.816496580927726 0.0) 0.4    0.046  0.4238)))


;;
;; Two-site model water :
;;
;; Dyer, K. M., Perkyns, J. S., Stell, G. & Montgomery Pettitt, B.
;; Site-renormalised molecular fluid theory: on the utility of a
;; two-site model of water. Molecular Physics 107, 423-431, 2009.
;; http://dx.doi.org/10.1080/00268970902845313
;;
;; sigma_H = 0.0   epsilon_H = 0.0   q_H = 0.38
;; sigma_O = 3.16  epsilon_O = 78 K  q_O =-0.38
;; r_OH = 1.0
;;
("two-site water"
 (("O" (0.0 0.0 0.0) 3.16 0.1549 -0.38)
  ("H" (1.0 0.0 0.0) 0.0 0.0 0.38)))

;;;
;;; Fake single-site  LJ water model. With  β = 1.6889  the reduced LJ
;;; temperature,  T* =  T/ε, is  ~3.8 (the  corresponding β*  =  βε is
;;; ~0.26) The reduced density, ρ* =  ρσ³, at water number density ρ =
;;; 0.033427745 is ~1.05. For  reference, the Wigner-Seitz radius rs =
;;; 1.92575678516970 A for this density.
;;;
("OW"
 (("OW" (0.0 0.0 0.0) 3.16 0.1549 0.0)))
("OW2"
 (("OW" (0.0 0.0  2.0) 3.16 0.1549 0.0)
  ("OW" (0.0 0.0 -2.0) 3.16 0.1549 0.0)))
("OW2-UA"
 (("UA" (0.0 0.0 0.0) 3.16 0.6196 0.0))) ; united atom 4x epsilon

;;
;; Transferable  Potentials  for  Phase Equilibria.   1.   United-Atom
;; Description of  n-Alkanes. Marcus G. Martin and  J.  Ilja Siepmann,
;; J.   Phys.    Chem.   B,  1998,  102  (14),   pp  2569-2577.   DOI:
;; 10.1021/jp972543+
;;
;; TABLE 1: Comparison of the Lennard-Jones Parameters for the OPLS
;;          [7] SKS [9, 10] and TraPPE Force Fields
;;
;;                     OPLS        SKS            TraPPE
;;               -------------- -------------- -------------
;; pseudoatom    ε/kB [K] σ [Å] ε/kB [K] σ [Å] ε/kB [K] σ [Å]
;; ----------------------------------------------------------
;; CH4           147.9    3.73  N/A      N/A   148      3.73
;; CH3 (ethane)  104.1    3.775 114      3.93   98      3.75
;; CH3 (n-alkane) 88.1    3.905 114      3.93   98      3.75
;; CH2            59.4    3.905 47       3.93   46      3.95
;;
;; Pseudoatoms are  connected by bonds with  a fixed length  of 1.53 Å
;; for  the  OPLS model  and  1.54  Å for  the  SKS  and TraPPE  force
;; fields. [...] The equilibrium angle θ0  is set to 112° for OPLS and
;; 114° for SKS and TraPPE.
;;
;;
;; Parametrisierung C3H8 nach Martin:  H=CH3 O=CH2
;; sigma_H = 3.75   epsilon_H = 0.194616  q_H = 0
;; sigma_O = 3.95  epsilon_O = 0.091350   q_O =0
;; r_OH = 1.54
;; theta_HOH = 114.0 => r_HH = 2.5831
;; mass C3H8 =  44.0099 u
;; density = 0.5077 kg/L =  0.3011071 u/A^3 => 0.0068418038 /A^3
;; temperature : T= 200 K  => 0.3971755, => beta =2.5178
;;
;;
("propane"
 (("CH2" (+0.41937205696157085365 0.0 0.0) 3.95 0.091350 0.0)
  ("CH3" (-0.41937205696157085365 +1.29155267463595300518 0.0) 3.75 0.194616 0.0)
  ("CH3" (-0.41937205696157085365 -1.29155267463595300518 0.0) 3.75 0.194616 0.0)))

;;
;; By comparing  atom by  atom, we  found that FFs  for atom  sites in
;; methanol,  butanoic  acid and  hexane  which  are  recorded in  the
;; original code  are the same  with those in OPLS  All-Atom (OPLS-AA)
;; force field, see the supporting  information of [1]. We added place
;; description  of the  matching  entry in  the  tables of  supporting
;; information of [1] for each type  of atom in solutes, hope that one
;; can  easily  build  new  solutes (alcohols,  carboxylic  acids  and
;; alkanes) by re-using these entries
;;
;; [1] "Development and Testing of the OPLS All-Atom Force Field on
;;     Conformational Energetics and Properties of Organic Liquids",
;;     W.  L. Jorgensen, D. S. Maxwell and J. Tirado-Rives, Journal of
;;     the American Chemical Society, 1996, 118 (45), 11225-11236
;;     http://dx.doi.org/10.1021/ja9621760
("methanol"
 (("C" (-0.748 -0.015 0.024) 3.5 0.066 0.145)  ; 18th entry of table 1
  ("HC1" (-1.293 -0.202 -0.901) 2.5 0.03 0.04) ; 17th entry of table 1
  ("HC2" (-1.263 0.754 0.6) 2.5 0.03 0.04)
  ("HC3" (-0.699 -0.934 0.609) 2.5 0.03 0.04)
  ("O" (0.558 0.42 -0.278) 3.12 0.17 -0.683)   ; 15th entry of table 1
  ("OH" (0.716 1.404 0.137) 0.4 0.04 0.418)))  ; 16th entry of table 1
;; FIXME:  note that LJ  parameters for  "H" are  actually 0.0  in the
;; table,  here  small  values  are assigned  possibly  for  numerical
;; reason, similar to what has been done for "H" in "modified TIP3P"

;; H1 sigma and epsilon adopted:
("butanoic acid"
 (("C1" (1.422 -0.017 0.0) 3.75 0.105 0.52) ; 19th entry of table 5
  ("O1" (1.422 1.353 0.0) 2.96 0.21 -0.44)  ; 20th entry of table 5
  ("O2" (2.643 -0.722 0.0) 3.0 0.17 -0.53)  ; 21th entry of table 5
  ("C2" (0.1 -0.78 0.0) 3.5 0.066 -0.12)    ; 3rd entry of table 1
  ("C3" (-1.06 0.212 0.0) 3.5 0.066 -0.12)
  ("C4" (-2.381 -0.551 0.0) 3.5 0.066 -0.18) ; 2nd entry of table 1
  ;; FIXME: expected to match the last entry of table 5 but shouldn't
  ;; they be "0.4 0.046" even for numerical reason?
  ("OH" (3.21 -0.461 0.882) 3.4 0.046 0.45)
  ("H2" (0.043 -1.407 0.89) 2.5 0.03 0.06)  ; 6th entry of table 1
  ("H3" (0.043 -1.407 -0.89) 2.5 0.03 0.06)
  ("H4" (-1.002 0.838 -0.89) 2.5 0.03 0.06)
  ("H5" (-1.002 0.838 0.89) 2.5 0.03 0.06)
  ("H6" (-2.439 -1.178 0.89) 2.5 0.03 0.06)
  ("H7" (-2.439 -1.178 -0.89) 2.5 0.03 0.06)
  ("H8" (-3.21 0.157 0.0) 2.5 0.03 0.06)))

("hexane"
 (("C" (1.709 -2.812 0.0) 3.5 0.066 -0.18) ;; 2nd entry of table 1
  ("C" (1.684 -1.278 0.0) 3.5 0.066 -0.12) ;; 3rd entry of table 1
  ("C" (0.245 -0.753 0.0) 3.5 0.066 -0.12)
  ("C" (0.241 0.779 0.0) 3.5 0.066 -0.12)
  ("C" (-1.198 1.304 0.0) 3.5 0.066 -0.12)
  ("C" (-1.206 2.834 0.0) 3.5 0.066 -0.18)
  ("H" (2.236 -3.164 0.887) 2.5 0.03 0.06) ;; 6th entry of table 1
  ("H" (2.232 -3.164 -0.89) 2.5 0.03 0.06)
  ("H" (0.691 -3.204 0.003) 2.5 0.03 0.06)
  ("H" (2.202 -0.914 -0.888) 2.5 0.03 0.06)
  ("H" (2.201 -0.914 0.89) 2.5 0.03 0.06)
  ("H" (-0.273 -1.115 0.889) 2.5 0.03 0.06)
  ("H" (-0.272 -1.115 -0.89) 2.5 0.03 0.06)
  ("H" (0.757 1.142 -0.89) 2.5 0.03 0.06)
  ("H" (0.757 1.141 0.89) 2.5 0.03 0.06)
  ("H" (-1.716 0.944 0.89) 2.5 0.03 0.06)
  ("H" (-1.716 0.944 -0.89) 2.5 0.03 0.06)
  ("H" (-0.696 3.204 -0.89) 2.5 0.03 0.06)
  ("H" (-0.696 3.204 0.89) 2.5 0.03 0.06)
  ("H" (-2.236 3.19 0.0) 2.5 0.03 0.06)))

;;;
;;; FIXME: I am still not quite sure  about the relation of R = 1.58 A
;;; as a  property of  Uranium atom [GW96]  and the form  of site-site
;;; pair interaction potentials.  Every single paper I checked implies
;;; that  such a radius  (or rather  an *average*  of two  such radii)
;;; corresponds to a minumum of the pair interaction potential:
;;;
;;;   - Most  of them quote the shape of the  potential as given below
;;;     with r = R being the minumum.
;;;
;;;   - Some  of them either imply, or  explicitly state that combined
;;;     radius  is an  (arithmetic)  average of  the  two radii,  e.g.
;;;     Refs. [GW96] and also a much more recent [KC13].
;;;
;;; For instance, according to the the notes to Table 1 in Ref. [KC13]
;;; the Ow in SPC water model is  characterized by R ~ 1.78 A.  A pair
;;; U-Ow would then be characterized by R ~ 1.68 A and a pair of Ow-Ow
;;; by R ~ 1.78  A. I think it is WRONG to  interpret these numbers as
;;; minimuma  of  respective pair  potential  terms  --- because  both
;;; numbers are too small:
;;;
;;;   - The Ow-Ow water peak in RDF between two electrostatically
;;;     repulsive centers is beyond 3 A.
;;;
;;;   - The U-Ow bonds between electrostatically attractive centers
;;;     are of the order of 2.4 A e.g. according to the very same
;;;     [KC13].
;;;
;;;   - The value of σ(Ow) sometimes used to characterize the SPC
;;;     model of water is about 3.16 A, the corresponding σ for a
;;;     Ow-Ow pair is thus again 3.16 A.  So, by the trivial math the
;;;     corresponding minumum is at 2^(1/6) σ which is even more.
;;;
;;; The text accompanying this  particular choice of atomic parameters
;;; for pair  interactions is  sometimes unclear, misleading,  or even
;;; wrong.   FIXME:  My guess  is  that  the  primary reason  for  the
;;; confusion  is operating  with the  atomic radii  R =  σ  / 2^(5/6)
;;; defined so that the minima of pair term for the like centers is at
;;; 2^(1/6) σ =  2 R.  The combination rule to  obtain the location of
;;; the minimum  of the  pair term is  thus never an  *averaging*, but
;;; rather  *adding*  two  atomic   radii  which  is  not  being  made
;;; sufficiently explicit in the text.
;;;
;;; Here a few  argument that make me think so:
;;;
;;;   - The radius R = 1.7412 A as quoted for Sr^2+ ion in Ref. [GW96]
;;;     as derived from the 6-12 parameters of Aaqvist, C6 = 20.54^2
;;;     kcal * A^6 and C12 = 613.5^2 kcal * A^12, Table II Ref.
;;;     [Aaqvist90] (but note the squares), happens to compare very
;;;     well with 1.7413 A obtained for a minumum of the non-bonding
;;;     VdW term in Sr^2+ - Sr^2+ interaction.
;;;
;;;   - With σ(Ow) of SPC water being 3.166 A (which may be derived
;;;     from the original Ref.  [BGS97]) the corresponding atomic
;;;     parameter would be R = 3.166 / 2^(5/6) = 1.777 A exactly as
;;;     quoted in Ref. [KC13].
;;;
;;; Original Guilbaud-Wipff force field parametrization uses this from
;;; of LJ potential between two species:
;;;
;;;                    12          6
;;;   v  (r) = ε [(R/r)   - 2 (R/r) ]
;;;    LJ
;;;
;;; Parameter R is the location of the minimum of that function and -ε
;;; is the  value at that  point.  Since R  and σ are  proportional to
;;; each other we may (and do) treat R-parameters as just another type
;;; of length units specified by (r0 ...)  form.
;;;
;;; To add to  the confusion the atoms are  characterized by an atomic
;;; (VdW) radius R to be combined into the pair interaction parameters
;;; by adding them (not by averaging as some of the texts erroneousely
;;; imply or state):
;;;
;;;   R   =  R  + R
;;;    ab     a    b
;;;
;;; Hence,  the atomic  radii  together with  the overall  interaction
;;; strengths are related to the atomic parameters σ and ε by
;;;
;;;            5/6
;;;   σ = R * 2
;;;   ε = ε
;;;
;;; The minimum of the  pair-interaction potential of the like-species
;;; is at
;;;
;;;    1/6          1/6 + 5/6
;;;   2    σ = R * 2          = R + R
;;;
;;;   2+
;;; Sr  :  R = 1.7412 A, ε = 0.1182 kcal / mol (σ ~ 3.1025 A)
;;;
;;; [GW96] Force field representation of the UO22+ cation from free
;;;   energy MD simulations in water. Tests on its 18-crown-6 and NO3−
;;;   adducts, and on its calix[6]arene6− and CMPO complexes,
;;;   P. Guilbaud, G. Wipff, Journal of Molecular Structure: THEOCHEM,
;;;   Volume 366, Issues 1–2, 31 July 1996, Pages
;;;   55–63. http://dx.doi.org/10.1016/0166-1280(96)04496-X
;;;
;;; [KC13] Structure, Kinetics, and Thermodynamics of the Aqueous
;;;   Uranyl(VI) Cation, Sebastien Kerisit, Chongxuan Liu,
;;;   J. Phys. Chem. A, 2013, 117 (30), pp 6421–6432,
;;;   http://dx.doi.org/10.1021/jp404594p
;;;
;;; [Aaqvist90] Ion-water interaction potentials derived from free
;;;   energy perturbation simulations, Johan. Aaqvist, J. Phys. Chem.,
;;;   1990, 94 (21), pp 8021–8024,
;;;   http://dx.doi.org/10.1021/j100384a009
;;;
("Sr2+, GW96"                           ; use "Sr2+" instead
 (("Sr2+" (0.0 0.0 0.0) (r0 (A 1.7412)) (kcal 0.1182) 2.0)))

("uranyl, GW96"
 (("U" (0.0 0.0 0.0) (r0 (A 1.58)) (kcal 0.4) 2.5)
  ("O" (-1.80 0.0 0.0) (r0 (A 1.75)) (kcal 0.2) -0.25)
  ("O" (1.80 0.0 0.0) (r0 (A 1.75)) (kcal 0.2) -0.25)))

;;;
;;; Lennard-Jones  parameters (R₀  and ε)  for Alkali  and  halide ion
;;; developed in [JC08],  which were later used in  [JLC13] to compare
;;; DRISM and MD results. There are three sets of parameters in [JC08]
;;; for use with  TIP3P, TIP4P and SPC/E respectively,  but only those
;;; with SPC/E  were later applied  in [JLC13]. Results  obtained with
;;; our  1D-RISM  code  (with   dielectric  correction)  are  in  good
;;; agreement with [JLC13].
;;;
;;; [JC08] Determination of Alkali and Halide Monovalent Ion
;;;   Parameters for Use in Explicitly Solvated Biomolecular
;;;   Simulations, I. S.  Joung and T. E. Cheatham, The Journal of
;;;   Physical Chemistry B, 2008, 112 (30), pp 9020-9041
;;;   http://dx.doi.org/10.1021/jp8001614
;;;
;;; [JLC13] Simple electrolyte solutions: Comparison of DRISM and
;;;   molecular dynamics results for alkali halide solutions, In Suk
;;;   Joung, Tyler Luchko, and David A. Case, J. Chem. Phys. 138,
;;;   044103 (2013), http://dx.doi.org/10.1063/1.4775743
;;;
("Li+, JC08, SPC/E" (("Li" (0.0 0.0 0.0) (r0 (A 0.791)) (kcal 0.3367344) 1.0)))
("Na+, JC08, SPC/E" (("Na" (0.0 0.0 0.0) (r0 (A 1.212)) (kcal 0.3526418) 1.0)))
("K+, JC08, SPC/E" (("K" (0.0 0.0 0.0) (r0 (A 1.593)) (kcal 0.4297054) 1.0)))
("Rb+, JC08, SPC/E" (("Rb" (0.0 0.0 0.0) (r0 (A 1.737)) (kcal 0.4451036) 1.0)))
("Cs+, JC08, SPC/E" (("Cs" (0.0 0.0 0.0) (r0 (A 2.021)) (kcal 0.0898565) 1.0)))
("F-, JC08, SPC/E" (("F" (0.0 0.0 0.0) (r0 (A 2.257)) (kcal 0.0074005) -1.0)))
("Cl-, JC08, SPC/E" (("Cl" (0.0 0.0 0.0) (r0 (A 2.711)) (kcal 0.0127850) -1.0)))
("Br-, JC08, SPC/E" (("Br" (0.0 0.0 0.0) (r0 (A 2.751)) (kcal 0.0269586) -1.0)))
("I-, JC08, SPC/E" (("I" (0.0 0.0 0.0) (r0 (A 2.919)) (kcal 0.0427845) -1.0)))

;;;
;;; Lennard-Jones parameters and calculated hydration free energies of
;;; alkali-metal- and alkaline-earth-metal ions [Aaqvist90].  Here x =
;;; √C6 and  y = √C12  --- Aaqvist used the  following parametrization
;;; for pair VdW interactions (albeit a different notation):
;;;
;;;                     6            12
;;;   v  (r) = - x x / r  +  y y  / r
;;;    LJ         i j         i j
;;;
;;; Here the numbers from OCRed  PDF (Table 11 Ref. [Aaqvist90]), note
;;; the order of y and x:
;;;
;;;   Ion       y     x     dG
;;;   ------------------------
;;;   Li+    25.0  2.60 -122.2
;;;   Na+   143.7  3.89  -98.5
;;;   K+    522.7  4.35  -80.9
;;;   Rb+   824.4  4.64  -75.5
;;;   Cs+  1647.9  5.44  -67.7
;;;   ------------------------
;;;   Mg2+   37.0  8.32 -455.9
;;;   Ca2+  264.1 18.82 -380.6
;;;   Sr2+  613.5 20.54 -345.9
;;;   Ba2+ 1341.5 24.13 -314.6
;;;
;;; The basic units are kcal and angstrom so that the actual units for
;;; x is (kcal  * A^6)^(1/2) and for y (kcal *  A^12)^(1/2).  dG is in
;;; kcals, naturally.  The actual C6  and C12 coefficients may are the
;;; products (or squares) of the table entries.
;;;
;;; FIXME:  does   anyone  dare  to  introduce   (kcal*A^6  ...)   and
;;; (kcal*A^12 ...) forms?
;;;
;;; Conversion from either AB-, or 6-12-, or this weired "3-6" (FIXME:
;;; does it  have a  name?)  parametrization to  σε-parametrization is
;;; not just two separate unit  conversion of two numbers. Here a form
;;; with two  arguments appears  in place of  two LJ  parameters. This
;;; form   is   interpreted   by   tabulate-ff  function   called   by
;;; find-molecule. Note that  (^2 ...)  form is not  a unit conversion
;;; either --- it squares the argument.
;;;
("Li+"  (("Li+"  (0.0 0.0 0.0) (c6/c12 (^2  2.60) (^2   25.0)) 1.0)))
("Na+"  (("Na+"  (0.0 0.0 0.0) (c6/c12 (^2  3.89) (^2  143.7)) 1.0)))
("K+"   (("K+"   (0.0 0.0 0.0) (c6/c12 (^2  4.35) (^2  522.7)) 1.0)))
("Rb+"  (("Rb+"  (0.0 0.0 0.0) (c6/c12 (^2  4.64) (^2  824.4)) 1.0)))
("Cs+"  (("Cs+"  (0.0 0.0 0.0) (c6/c12 (^2  5.44) (^2 1647.9)) 1.0)))
("Mg2+" (("Mg2+" (0.0 0.0 0.0) (c6/c12 (^2  8.32) (^2   37.0)) 2.0)))
("Ca2+" (("Ca2+" (0.0 0.0 0.0) (c6/c12 (^2 18.82) (^2  264.1)) 2.0)))
("Sr2+" (("Sr2+" (0.0 0.0 0.0) (c6/c12 (^2 20.54) (^2  613.5)) 2.0)))
("Ba2+" (("Ba2+" (0.0 0.0 0.0) (c6/c12 (^2 24.13) (^2 1341.5)) 2.0)))

;;;                    2+
;;; Force fields for UO  ion with  the SPC/Fw, TIP3P,  TIP4p and TIP5P
;;;                    2
;;; water models. See [1] for more details. For the geometry, since it
;;; only  mentions that  U-O  bond length  is  1.76 Å,  we choose  the
;;; orientation along x-axis in consistent with other solutes
;;;
;;; BL: Made stupid mistakes when tabulating parameters from literature,
;;; so it's better to copy the original table here to avoid any
;;; human-made mistake in future. Note that these are pair parameters
;;;
;;;                         Bare Uranyl
;;;                     ------------------------------
;;;                     SPC/Fw   TIP3P   TIP4P   TIP5P
;;; --------------------------------------------------
;;; qUU (e)              3.08    3.08    3.08    3.08
;;; qOU (e)              -0.54   -0.54   -0.54   -0.54
;;; σUU-OW (Å)           2.32    2.35    2.31    2.37
;;; σOU-OW (Å)           2.39    2.49    2.24    2.50
;;; εUU-OW (kJ/mol)      21.92   16.32   21.92   19.00
;;; εOU-OW (kJ/mol)      2.72    0.84    8.08    1.72
;;;
;;;             Solvated Uranyl (Recommended for Use)
;;;             --------------------------------------
;;;                     SPC/Fw   TIP3P   TIP4P   TIP5P
;;; --------------------------------------------------
;;; qUU (e)              2.50    2.50    2.50    2.50
;;; qOU (e)              -0.25   -0.25   -0.25   -0.25
;;; σUU-OW (Å)           3.25    3.25    2.78    3.26
;;; σOU-OW (Å)           3.00    3.00    3.07    2.92
;;; εUU-OW (kJ/mol)      0.27    0.27    1.39    0.26
;;; εOU-OW (kJ/mol)      1.08    1.08    0.85    1.47
;;;
;;; By applying LB mixing rule (which is used in our code), we could
;;; derive the LJ parameter of each site in UO2+ (using thermochemical
;;; kcal of 4.184 kJ).
;;;
;;;                       Solvated Uranyl
;;;             -----------------------------------
;;;                     SPC/Fw         TIP3P
;;; -----------------------------------------------
;;; σUU-OW (Å)           3.25          3.25
;;; σOU-OW (Å)           3.00          3.00
;;; σOW (Å)              3.1656        3.1506
;;; σUU (Å)              3.3344        3.3494
;;; σOU (Å)              2.8344        2.8494
;;; εUU-OW (kcal/mol)    0.0645315488  0.0645315488
;;; εOU-OW (kcal/mol)    0.258126195   0.258126195
;;; εOW (kcal/mol)       0.1553        0.1521
;;; εUU (kcal/mol)       0.0268146863  0.0273788349
;;; εOU (kcal/mol)       0.4290349811  0.4380613581
;;;
;;; [1] Force Field Development for Actinyl Ions via Quantum
;;;     Mechanical Calculations: An Approach to Account for Many Body
;;;     Solvation Effects, Rai, N.; Tiwari, S. P.; Maginn, E. J., The
;;;     Journal of Physical Chemistry B 2012, 116 (35), 10885-10897.
;;;     http://dx.doi.org/10.1021/jp3028275
;;;
("uranyl, SPC"
 (("U" (0.0 0.0 0.0) 3.3344 0.0268146863 2.5)
  ("O" (-1.76 0.0 0.0) 2.8344 0.4290349811 -0.25)
  ("O" (1.76 0.0 0.0) 2.8344 0.4290349811 -0.25)))

("uranyl, TIP3P"
 (("U" (0.0 0.0 0.0) 3.3494 0.0273788349 2.5)
  ("O" (-1.76 0.0 0.0) 2.8494 0.4380613581 -0.25)
  ("O" (1.76 0.0 0.0) 2.8494 0.4380613581 -0.25)))

;;;
;;; Another set of  force field parameters which is  also developed by
;;; E.  J.  Maginn  et. al.,  see [PM13]  for more  details.  Explicit
;;; parameters could  be found in  the Gromacs input file  provided by
;;; the authors:
;;;
;;;   http://www.rsc.org/suppdata/cp/c3/c3cp52444b/c3cp52444b_2.pdf
;;;
;;; The  file specifies  both,  the pair  parameters for  uranyl-water
;;; interactions and atomic parameters  for uranyl and water sites. It
;;; seems that  the pair interaction parameters are  consistent (up to
;;; three digits) to  those derived from the atomic  parameters by the
;;; LB rule. Thus, this database uses the atomic parameters.
;;;
;;; [PM13] Development and application of effective pairwise
;;;   potentials for UO2n+, NpO2n+, PuO2n+, and AmO2n+ (n = 1, 2) ions
;;;   with water, Vladimir Pomogaev, Surya Prakash Tiwari, Neeraj Rai,
;;;   George S. Goff, Wolfgang Runde, William F. Schneider and Edward
;;;   J. Maginn, Phys. Chem. Chem. Phys., 2013, Advance Article
;;;   http://dx.doi.org/10.1039/C3CP52444B
;;;
("uranyl, PM13" ;;; FIXME: need better name
 (("U" (0.0 0.0 0.0) (A 2.95) (kJ 0.530) 2.5)
  ("O" (-1.761 0.0 0.0) (A 3.83) (kJ 0.057) -0.25)
  ("O" (1.761 0.0 0.0) (A 3.83) (kJ 0.057) -0.25)))

;;;
;;; Two models modified by Kerisit and Liu based on GW model to
;;; strengthen uranyl-water interaction.
;;;
;;; [KL13] Structure, Kinetics, and Thermodynamics of the Aqueous
;;;   Uranyl(VI) Cation, Kerisit, S. and C. Liu, The Journal of Physical
;;;   Chemistry A 2013, 117(30): 6421-6432.
;;;   http://dx.doi.org/10.1021/jp404594p
;;;
("uranyl, KL1"
 (("U" (0.0 0.0 0.0) (r0 (A 1.60)) (kcal 0.12) 3.25)
  ("O" (-1.761 0.0 0.0) (r0 (A 1.75)) (kcal 0.2) -0.625)
  ("O" (1.761 0.0 0.0) (r0 (A 1.75)) (kcal 0.2) -0.625)))

("uranyl, KL2"
 (("U" (0.0 0.0 0.0) (r0 (A 1.58)) (kcal 0.30) 3.50)
  ("O" (-1.761 0.0 0.0) (r0 (A 1.75)) (kcal 0.2) -0.75)
  ("O" (1.761 0.0 0.0) (r0 (A 1.75)) (kcal 0.2) -0.75)))

;;;
;;; The  default   geometry  is  the  same  in   several  uranyl  aqua
;;; complexes. This one  is refered to at several  other entries. With
;;; SPC/E    water    FF    uranyl/water   interaction    energy    is
;;; -214.652929634494  kcal. FIXME:  Note that  water geometry  is not
;;; SPC/E and not even the same for all of them.
;;;
("UO2_5H2O"
  (("U"  (-0.4193  0.2123  0.1018))
   ("OU"  (1.1822  0.9792  0.1743))
   ("OU"  (-2.0068  -0.5497  -0.1304))
   ("OW"  (-1.6172  2.2191  -0.7126))
   ("HW"  (-2.5526  2.2260  -0.9995))
   ("HW"  (-1.2731  3.1280  -0.8244))
   ("OW"  (0.9046  -1.8307  -0.3086))
   ("HW"  (0.5819  -2.7447  -0.4413))
   ("HW"  (1.8777  -1.8501  -0.4076))
   ("OW"  (-0.1489  0.1448  -2.3461))
   ("HW"  (-0.8081  -0.1703  -2.9968))
   ("HW"  (0.6309  0.4518  -2.8510))
   ("OW"  (-0.6837  1.7241  2.0738))
   ("HW"  (-1.4562  1.8695  2.6561))
   ("HW"  (0.0204  2.3354  2.3698))
   ("OW"  (-0.6191  -1.1831  2.1599))
   ("HW"  (-1.3870  -1.7658  2.3269))
   ("HW"  (-0.0073  -1.2855  2.9163))))

("UO2_5H2O, SPC"                     ; RM12 uranyl with SH-SPC/E water
 (geometry "UO2_5H2O")
  (force-field (uranyl/spc water/spc)
   ("U" U) ;; SPC uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, PRSPC"                   ; RM12 uranyl with PR-SPC/E water
  (geometry "UO2_5H2O")
  (force-field (uranyl/spc water/pr-spc)
   ("U" U) ;; SPC uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, KL1-PRSPC"                ; KL1 uranyl with PR-SPC/E water
  (geometry "UO2_5H2O")
  (force-field (uranyl/kl1 water/pr-spc)
   ("U" U) ;; KL1 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, KL2-PRSPC"                ; KL2 uranyl with PR-SPC/E water
  (geometry "UO2_5H2O")
  (force-field (uranyl/kl2 water/pr-spc)
   ("U" U) ;; KL2 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, GW96-PRSPC"                ; GW uranyl with PR-SPC/E water
  (geometry "UO2_5H2O")
  (force-field (uranyl/gw96 water/pr-spc)
   ("U" U) ;; GW96 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, PM13-PRSPC"                ; PM13 uranyl with PR-SPC/E water
  (geometry "UO2_5H2O")
  (force-field (uranyl/pm13 water/pr-spc)
   ("U" U) ;; PM13 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, TIP3P"
  (geometry "UO2_5H2O")
  (force-field (uranyl/t3p water/t3p)
   ("U" U) ;; TIP3P uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_4H2O, D4H"
  (("U"  (-0.4146  0.2139  0.0486))
   ("OU"  (1.1864  0.9806  0.1210))
   ("OU"  (-2.0157  -0.5528  -0.0239))
   ("OW"  (-1.4189  2.3544  -0.4701))
   ("HW"  (-2.3937  2.5109  -0.6290))
   ("HW"  (-0.9232  3.2179  -0.5624))
   ("OW"  (0.5895  -1.9284  0.5652))
   ("HW"  (0.0929  -2.7910  0.6612))
   ("HW"  (1.5655  -2.0872  0.7146))
   ("OW"  (-0.1044  -0.2153  -2.3146))
   ("HW"  (-0.7631  -0.6773  -2.9084))
   ("HW"  (0.7018  0.0418  -2.8475))
   ("OW"  (-0.7255  0.6434  2.4112))
   ("HW"  (-1.5362  0.3956  2.9416))
   ("HW"  (-0.0620  1.0949  3.0078))))

("UO2_4H2O, D4H, KL2-PRSPC"
 (geometry "UO2_4H2O, D4H")
  (force-field (uranyl/kl2 water/pr-spc)
   ("U" U) ;; KL2 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_4H2O, D4H, KL1-PRSPC"
 (geometry "UO2_4H2O, D4H")
  (force-field (uranyl/kl1 water/pr-spc)
   ("U" U) ;; KL1 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_4H2O, D4H, PM13-PRSPC"
 (geometry "UO2_4H2O, D4H")
  (force-field (uranyl/pm13 water/pr-spc)
   ("U" U) ;; PM13 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_4H2O, D4H, PRSPC"
 (geometry "UO2_4H2O, D4H")
  (force-field (uranyl/spc water/pr-spc)
   ("U" U) ;; SPC uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_4H2O, D4H, GW96-PRSPC"
 (geometry "UO2_4H2O, D4H")
  (force-field (uranyl/gw96 water/pr-spc)
   ("U" U) ;; GW96 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, CS"
  (("U"  (-0.4146  0.2139  0.0486))
   ("OU"  (1.1864  0.9806  0.1210))
   ("OU"  (-2.0157  -0.5528  -0.0239))
   ("OW"  (-1.6850  2.1673  -0.8016))
   ("HW"  (-2.6429  2.1084  -1.0826))
   ("HW"  (-1.3570  3.1032  -0.9296))
   ("OW"  (0.8901  -1.8389  -0.4262))
   ("HW"  (0.5416  -2.7748  -0.4760))
   ("HW"  (1.8738  -1.8394  -0.6060))
   ("OW"  (-0.2050  0.1167  -2.4138))
   ("HW"  (-0.8829  -0.2596  -3.0452))
   ("HW"  (0.5771  0.4629  -2.9319))
   ("OW"  (-0.7306  1.7851  1.9526))
   ("HW"  (-1.5488  1.9426  2.5056))
   ("HW"  (0.0146  2.3612  2.2884))
   ("OW"  (-0.3670  -1.1429  2.1334))
   ("HW"  (-1.1495  -1.6887  2.4329))
   ("HW"  (0.3702  -1.2208  2.8045))))

("UO2_5H2O, CS, PRSPC"
 (geometry "UO2_5H2O, CS")
  (force-field (uranyl/spc water/pr-spc)
   ("U" U) ;; SPC uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, CS, KL2-PRSPC"
 (geometry "UO2_5H2O, CS")
  (force-field (uranyl/kl2 water/pr-spc)
   ("U" U) ;; KL2 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, CS, KL1-PRSPC"
 (geometry "UO2_5H2O, CS")
  (force-field (uranyl/kl1 water/pr-spc)
   ("U" U) ;; KL1 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, CS, GW96-PRSPC"
 (geometry "UO2_5H2O, CS")
  (force-field (uranyl/gw96 water/pr-spc)
   ("U" U) ;; GW96 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_5H2O, CS, PM13-PRSPC"
 (geometry "UO2_5H2O, CS")
  (force-field (uranyl/pm13 water/pr-spc)
   ("U" U) ;; PM13 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_6H2O, D6H"
  (("U"  (0.0000  0.0000  0.0000))
   ("OU"  (0.0000  0.0000  1.7700))
   ("OU"  (0.0000  0.0000  -1.7700))
   ("OW"  (-2.6162  0.0000  0.0000))
   ("HW"  (-3.1936  -0.0000  0.8165))
   ("HW"  (-3.1936  -0.0000  -0.8165))
   ("OW"  (1.3082  -2.2673  0.0000))
   ("HW"  (1.5970  -2.7673  0.8165))
   ("HW"  (1.5970  -2.7673  -0.8165))
   ("OW"  (1.3082  2.2673  0.0000))
   ("HW"  (1.5970  2.7673  0.8165))
   ("HW"  (1.5970  2.7673  -0.8165))
   ("OW"  (2.6162  0.0000  -0.0000))
   ("HW"  (3.1936  0.0000  -0.8165))
   ("HW"  (3.1936  0.0000  0.8165))
   ("OW"  (-1.3082  2.2673  -0.0000))
   ("HW"  (-1.5970  2.7673  -0.8165))
   ("HW"  (-1.5970  2.7673  0.8165))
   ("OW"  (-1.3082  -2.2673  -0.0000))
   ("HW"  (-1.5970  -2.7673  -0.8165))
   ("HW"  (-1.5970  -2.7673  0.8165))))

("UO2_6H2O, D6H, KL2-PRSPC"
 (geometry "UO2_6H2O, D6H")
  (force-field (uranyl/kl2 water/pr-spc)
   ("U" U) ;; KL2 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_6H2O, D6H, KL1-PRSPC"
 (geometry "UO2_6H2O, D6H")
  (force-field (uranyl/kl1 water/pr-spc)
   ("U" U) ;; KL1 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_6H2O, D6H, GW96-PRSPC"
 (geometry "UO2_6H2O, D6H")
  (force-field (uranyl/gw96 water/pr-spc)
   ("U" U) ;; GW96 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_6H2O, D6H, PM13-PRSPC"
 (geometry "UO2_6H2O, D6H")
  (force-field (uranyl/pm13 water/pr-spc)
   ("U" U) ;; PM13 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_6H2O, D6H, PRSPC"
 (geometry "UO2_6H2O, D6H")
  (force-field (uranyl/spc water/pr-spc)
   ("U" U) ;; SPC uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_6H2O, D3D"
  (("U"  (-0.1098  0.0247  0.2657))
   ("OU"  (1.2099  1.1978  0.3124))
   ("OU"  (-1.5524  -0.9973  0.1654))
   ("OW"  (-1.8010  1.5327  -0.8009))
   ("HW"  (-2.7284  1.2170  -1.0015))
   ("HW"  (-1.6949  2.4762  -1.1147))
   ("OW"  (0.5409  -2.1893  -0.8350))
   ("HW"  (-0.1429  -2.8938  -1.0250))
   ("HW"  (1.4233  -2.4711  -1.2116))
   ("OW"  (0.5353  0.3582  -2.1448))
   ("HW"  (0.1153  -0.0223  -2.9687))
   ("HW"  (1.2349  1.0239  -2.4044))
   ("OW"  (-0.8422  1.7093  1.9783))
   ("HW"  (-1.7262  1.7916  2.4385))
   ("HW"  (-0.3007  2.5330  2.1468))
   ("OW"  (-0.4315  -1.0882  2.5460))
   ("HW"  (-1.2186  -1.6839  2.7060))
   ("HW"  (0.0066  -0.8730  3.4189))
   ("OW"  (2.0224  -1.0724  1.2049))
   ("HW"  (2.8593  -0.5401  1.3322))
   ("HW"  (2.1367  -1.9767  1.6163))))

("UO2_6H2O, D3D, KL2-PRSPC"
 (geometry "UO2_6H2O, D3D")
  (force-field (uranyl/kl2 water/pr-spc)
   ("U" U) ;; KL2 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_6H2O, D3D, KL1-PRSPC"
 (geometry "UO2_6H2O, D3D")
  (force-field (uranyl/kl1 water/pr-spc)
   ("U" U) ;; KL1 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_6H2O, D3D, GW96-PRSPC"
 (geometry "UO2_6H2O, D3D")
  (force-field (uranyl/gw96 water/pr-spc)
   ("U" U) ;; GW96 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_6H2O, D3D, PM13-PRSPC"
 (geometry "UO2_6H2O, D3D")
  (force-field (uranyl/pm13 water/pr-spc)
   ("U" U) ;; PM13 uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))

("UO2_6H2O, D3D, PRSPC"
 (geometry "UO2_6H2O, D3D")
  (force-field (uranyl/spc water/pr-spc)
   ("U" U) ;; SPC uranyl
   ("OU" OU)
   ("HW" HW)
   ("OW" OW)))
;;;
;;; This  is a fake  solute which  tabulates some  half-way meaningful
;;; force field parameters for various sites:
;;;
("bgy3d"
 (("H" (0.0 0.0 0.0) 2.735 0.03971 0.2)
  ("Cl" (0 0.0 0.0) 3.353 0.51434 -0.2)
  ("Pd" (0.0 0.0 0.0) 2.5 0.1 0.0)))    ; fake numbers

;; Tabulating OPLS-AA force field from database in TINKER program
;; explicit parameters (sigma, epsilon and charge) are replaced by atom
;; indices which are used in oplsaa.prm
("methanol, oplsaa"
 (("CT" (-0.748 -0.015 0.024))
  ("HC" (-1.293 -0.202 -0.901))
  ("HC" (-1.263 0.754 0.6))
  ("HC" (-0.699 -0.934 0.609))
  ("OH" (0.558 0.42 -0.278))
  ("HO" (0.716 1.404 0.137)))
 (force-field (oplsaa)
  ("CT" 96) ;; "Alcohol CH3OH & RCH2OH"
  ("HC" 95) ;;"Methanol CH3-"
  ("OH" 93) ;;"Alcohol -OH"
  ("HO" 94))) ;;"Alcohol -OH", vdw parameters are
                 ;; 0.0 in original database

("methanol-staggered, oplsaa"
 (("CT" (-0.748 -0.015 0.024))
  ("HC" (-1.293 -0.202 -0.901))
  ("HC" (-1.263 0.754 0.6))
  ("HC" (-0.699 -0.934 0.609))
  ("OH" (0.558 0.42 -0.278))
  ("HO" (1.06949445 -0.34100606 -0.84781418)))
 (force-field (oplsaa)
  ("CT" 96)
  ("HC" 95)
  ("OH" 93)
  ("HO" 909)))

("butanoic acid, oplsaa"
 (("C" (1.422 -0.017 0.0))
  ("O" (1.422 1.353 0.0))
  ("OH" (2.643 -0.722 0.0))
  ("CT2" (0.1 -0.78 0.0))
  ("CT2" (-1.06 0.212 0.0))
  ("CT3" (-2.381 -0.551 0.0))
  ("HO" (3.21 -0.461 0.882))
  ("HC" (0.043 -1.407 0.89))
  ("HC" (0.043 -1.407 -0.89))
  ("HC" (-1.002 0.838 -0.89))
  ("HC" (-1.002 0.838 0.89))
  ("HC" (-2.439 -1.178 0.89))
  ("HC" (-2.439 -1.178 -0.89))
  ("HC" (-3.21 0.157 0.0)))
 (force-field (oplsaa)
  ("C" 206)  ;; "Carboxylic Acid -COOH"
  ("O" 207)  ;; "Carboxylic Acid C=O"
  ("OH" 208) ;; "Carboxylic Acid -OH"
  ("CT2" 78)  ;; "Alkane -CH2-"
  ("CT3" 77)  ;; "Alkane CH3-"
  ("HO" 209) ;; "Carboxylic Acid -COOH", original vdw parameters are
                 ;; 0.0
  ("HC" 82))) ;; "Alkane H-C"

("hexane, oplsaa"
 (("CT3" (1.709 -2.812 0.0))
  ("CT2" (1.684 -1.278 0.0))
  ("CT2" (0.245 -0.753 0.0))
  ("CT2" (0.241 0.779 0.0))
  ("CT2" (-1.198 1.304 0.0))
  ("CT3" (-1.206 2.834 0.0))
  ("HC" (2.236 -3.164 0.887))
  ("HC" (2.232 -3.164 -0.89))
  ("HC" (0.691 -3.204 0.003))
  ("HC" (2.202 -0.914 -0.888))
  ("HC" (2.201 -0.914 0.89))
  ("HC" (-0.273 -1.115 0.889))
  ("HC" (-0.272 -1.115 -0.89))
  ("HC" (0.757 1.142 -0.89))
  ("HC" (0.757 1.141 0.89))
  ("HC" (-1.716 0.944 0.89))
  ("HC" (-1.716 0.944 -0.89))
  ("HC" (-0.696 3.204 -0.89))
  ("HC" (-0.696 3.204 0.89))
  ("HC" (-2.236 3.19 0.0)))
 (force-field (oplsaa)
  ("CT3"  77) ;; "Alkane CH3-"
  ("CT2"  78) ;; "Alkane -CH2-"
  ("HC"  82))) ;; "Alkane H-C"

;; ********************
;; *                  *
;; *  NIST molecules  *
;; *                  *
;; ********************
;; Original molecule structure files are collected in ../nist-molecules
;; Currently the atomic symbols are converted to match the oplsaa force
;; field database by "artificial recognition" with the help of ViewerPro
;; There must be a better way to assign the atomic symbol automatically,
;; e.g., by identifying the functional group, calculating the bond order
;;
;; 1. Alcohols:
;;
("methanol, nist"
 (("HC"  (0.2453  0.8386  1.6056))
  ("CT"  (0.6776  0.9803  0.6074))
  ("HC"  (0.0000  1.5869  0.0000))
  ("HC"  (0.8133  0.0000  0.1338))
  ("OH"  (1.8631  1.7142  0.6464))
  ("HO"  (2.4856  1.2216  1.1660)))
 (force-field (oplsaa)
  ("CT" 96) ;; "Alcohol CH3OH & RCH2OH"
  ("HC" 95) ;; "Methanol CH3-"
  ("OH" 93) ;; "Alcohol -OH"
  ("HO" 94)))

("ethanol, nist"
 (("CT3"  (1.0303  0.8847  0.9763))
  ("CT"  (1.8847  1.9889  1.5717))
  ("OH"  (3.1883  1.4807  1.7425))
  ("HC"  (0.0000  1.2330  0.8324))
  ("HC"  (0.9949  0.0000  1.6255))
  ("HC"  (1.4097  0.5558  0.0000))
  ("HC"  (1.9050  2.8742  0.9059))
  ("HC"  (1.4753  2.3225  2.5456))
  ("HO"  (3.7056  2.1820  2.1139)))
 (force-field (oplsaa)
  ("CT" 96)    ;; "Alcohol CH3OH & RCH2OH"
  ("CT3" 77)   ;; "Alkane CH3-"
  ("OH" 93)    ;; "Alcohol -OH"
  ("HO" 94)    ;; "Alcohol -OH"
  ("HC" 82)))  ;; "Alkane H-C"

("propanol, nist"
 (("CT"  (0.7713  1.5705  1.3838))
  ("CT2"  (2.1696  1.0226  1.0958))
  ("CT3"  (3.2631  1.9141  1.6398))
  ("OH"  (0.3563  1.3950  2.7118))
  ("HC"  (0.7013  2.6399  1.1017))
  ("HC"  (0.0000  1.0213  0.8161))
  ("HC"  (2.2872  0.9154  0.0000))
  ("HC"  (2.2716  0.0000  1.5099))
  ("HC"  (3.2248  1.9921  2.7352))
  ("HC"  (3.1932  2.9340  1.2381))
  ("HC"  (4.2564  1.5270  1.3793))
  ("HO"  (1.0138  1.7939  3.2691)))
 (force-field (oplsaa)
  ("CT" 96)
  ("CT2" 78) ;; "Alkane -CH2-"
  ("CT3" 77)
  ("HC" 82)
  ("OH" 93)
  ("HO" 94)))

("butanol, nist"
 (("CT2"  (2.9651  2.0464  2.4042))
  ("CT2"  (2.3281  2.7103  1.1934))
  ("CT"  (0.8652  2.3158  0.9860))
  ("CT3"  (3.3109  0.5952  2.1559))
  ("OH"  (0.0000  2.7697  1.9922))
  ("HC"  (3.8836  2.6000  2.6826))
  ("HC"  (2.2886  2.1307  3.2779))
  ("HC"  (2.8987  2.4514  0.2796))
  ("HC"  (2.4074  3.8120  1.2923))
  ("HC"  (0.5089  2.6739  0.0000))
  ("HC"  (0.7244  1.2206  1.0177))
  ("HC"  (2.4210  0.0000  1.9086))
  ("HC"  (3.7712  0.1403  3.0423))
  ("HC"  (4.0193  0.4827  1.3240))
  ("HO"  (0.1229  3.7075  2.0732)))
 (force-field (oplsaa)
  ("CT" 96)
  ("CT2" 78)
  ("CT3" 77)
  ("HC" 82)
  ("OH" 93)
  ("HO" 94)))

("pentanol, nist"
 (("CT2"  (2.0264  2.2581  1.9056))
  ("CT2"  (3.2462  1.6488  1.2359))
  ("CT2"  (4.4464  1.6425  2.1693))
  ("CT"  (0.8116  2.2939  0.9786))
  ("CT3"  (5.6613  1.0380  1.5026))
  ("OH"  (0.3196  1.0290  0.6256))
  ("HC"  (2.2542  3.2923  2.2314))
  ("HC"  (1.7856  1.7003  2.8338))
  ("HC"  (3.0191  0.6161  0.9030))
  ("HC"  (3.4892  2.2109  0.3117))
  ("HC"  (4.6718  2.6757  2.5011))
  ("HC"  (4.2030  1.0807  3.0931))
  ("HC"  (1.0538  2.7449  0.0000))
  ("HC"  (0.0000  2.8939  1.4351))
  ("HC"  (5.4791  0.0000  1.1934))
  ("HC"  (5.9510  1.5990  0.6041))
  ("HC"  (6.5258  1.0323  2.1788))
  ("HO"  (0.1634  0.5491  1.4297)))
 (force-field (oplsaa)
  ("CT" 96)
  ("CT2" 78)
  ("CT3" 77)
  ("HC" 82)
  ("OH" 93)
  ("HO" 94)))

;;
;; 2. Carboxylate acids:
;;
("formic acid, nist"
 (("OH"  (1.9226  0.9205  1.8422))
  ("C"  (0.9364  0.3495  1.1289))
  ("O"  (1.2026  0.0000  0.0000))
  ("HO"  (2.7325  0.9756  1.3421))
  ("HCO"  (0.0000  0.2800  1.6921)))
 (force-field (oplsaa)
  ("C" 206)    ;; "Carboxylic Acid -COOH"
  ("O" 207)    ;; "Carboxylic Acid C=O"
  ("OH" 208)   ;; "Carboxylic Acid -OH"
  ("HO" 209)   ;; "Carboxylic Acid -COOH"
  ("HCO" 218)));; FIXME: "Aldehyde/Formamide HCO-"

("acetic acid, nist"
 (("CT3"  (0.7649  0.9627  1.0051))
  ("C"  (2.0422  1.7381  0.9144))
  ("OH"  (3.1225  1.0260  0.5122))
  ("O"  (2.2424  2.9166  1.1481))
  ("HC"  (0.0000  1.5240  1.5565))
  ("HC"  (0.3742  0.7552  0.0000))
  ("HC"  (0.9140  0.0000  1.5112))
  ("HO"  (3.8882  1.5907  0.4746)))
 (force-field (oplsaa)
  ("C" 206)
  ("CT3" 77)
  ("O" 207)
  ("OH" 208)
  ("HO" 209)
  ("HC" 82)))

;; Raw structure, coordinates of non-hydrogen atoms are from NIST, but
;; hydrogen atoms are generated by ViewerPro with C-H bond and O-H bond
;; length equal to 1.0 \AA. RISM/QM calculation of this strucure
;; diverged.
("propanoic acid, nist, raw"
 (("CT2"  (2.6488  2.4258  0.0002))
  ("C"  (1.3807  1.6210  0.0000))
  ("OH"  (1.4371  0.2295  0.0004))
  ("O"  (0.2857  2.1943  -0.0007))
  ("CT3"  (3.8613  1.4766  0.0006))
  ("HC"  (2.6776  3.0021  -0.8165))
  ("HC"  (2.6771  3.0023  0.8168))
  ("HO"  (0.5165  -0.1401  -0.0003))
  ("HC"  (4.7051  2.0131  0.0008))
  ("HC"  (3.8338  0.8996  -0.8157))
  ("HC"  (3.8334  0.8998  0.8171)))
 (force-field (oplsaa)
  ("C" 206)
  ("CT3" 77)
  ("CT2" 78)
  ("O" 207)
  ("OH" 208)
  ("HO" 209)
  ("HC" 82)))

;; Optimized in PG, angle and bond between non-hydrogen atoms are fixed.
;; C-H bond length is close to 1.10 \AA which is in consistent with
;; other carboxylate acid, but O-H bond length is 0.98 \AA which
;; slightly longer than 0.95 \AA in other molecules.
("propanoic acid, nist"
 (("OH"  (1.436311   0.222125  -0.018177))
  ("C"  (1.376332   1.613333   0.002239))
  ("O"  (0.279858   2.183798   0.007830))
  ("CT2"  (2.642361   2.421263   0.017252))
  ("CT3"  (3.857302   1.475244   0.005956))
  ("HO"  (0.513545  -0.119038  -0.025950))
  ("HC"  (2.633597   3.054770  -0.888930))
  ("HC"  (2.615882   3.086931   0.896967))
  ("HC"  (4.808148   2.038879  -0.032852))
  ("HC"  (3.819263   0.806527  -0.869903))
  ("HC"  (3.874655   0.840266   0.908270)))
 (force-field (oplsaa)
  ("C" 206)
  ("CT3" 77)
  ("CT2" 78)
  ("O" 207)
  ("OH" 208)
  ("HO" 209)
  ("HC" 82)))


("butanoic acid, nist"
 (("CT2"  (3.4057  1.8922  1.6847))
  ("C"  (4.6993  1.3329  1.1549))
  ("CT2"  (2.1503  1.2627  1.1037))
  ("OH"  (4.6276  0.2556  0.3391))
  ("O"  (5.8306  1.7309  1.3762))
  ("CT3"  (0.9097  1.9095  1.6789))
  ("HC"  (3.4147  1.7744  2.7880))
  ("HC"  (3.4050  2.9864  1.5011))
  ("HC"  (2.1505  1.3594  0.0000))
  ("HC"  (2.1380  0.1739  1.3082))
  ("HO"  (5.5056  0.0000  0.0723))
  ("HC"  (0.0000  1.4538  1.2674))
  ("HC"  (0.8647  1.8031  2.7710))
  ("HC"  (0.8690  2.9837  1.4539)))
 (force-field (oplsaa)
  ("C" 206)
  ("O" 207)
  ("OH" 208)
  ("HO" 209)
  ("CT2" 78)
  ("CT3" 77)
  ("HC" 82)))

("pentanoic acid, nist"
 (("CT3"  (0.9690  0.9499  1.0698))
  ("CT2"  (2.1126  1.8447  1.4909))
  ("CT2"  (3.4564  1.1801  1.2350))
  ("CT2"  (4.6009  2.1010  1.6329))
  ("C"  (5.9285  1.4318  1.3841))
  ("OH"  (6.4076  0.6912  2.4116))
  ("O"  (6.6289  1.4680  0.3886))
  ("HC"  (0.0000  1.4309  1.2543))
  ("HC"  (1.0152  0.7058  0.0000))
  ("HC"  (0.9735  0.0000  1.6208))
  ("HC"  (2.0574  2.8087  0.9469))
  ("HC"  (2.0170  2.1022  2.5646))
  ("HC"  (3.5191  0.2279  1.7986))
  ("HC"  (3.5430  0.9031  0.1652))
  ("HC"  (4.5558  3.0450  1.0542))
  ("HC"  (4.5113  2.3890  2.6998))
  ("HO"  (7.2417  0.3008  2.1695)))
 (force-field (oplsaa)
  ("C" 206)
  ("O" 207)
  ("OH" 208)
  ("HO" 209)
  ("CT2" 78)
  ("CT3" 77)
  ("HC" 82)))


;;
;; 3. Aldehydes:
;;
;; original molecule is translated to make C atom centeral
("formaldehyde, nist"
 (("O"  (-0.000 1.236 -0.001))
  ("C"  (0.000 0.000 0.000))
  ("HCO"  (-0.866 -0.500 0.001))
  ("HCO"  (0.866 -0.500 -0.000)))
 (force-field (oplsaa)
  ("C" 216)     ;; "Aldehyde/Acyl Halide C=O"
  ("O" 217)     ;; "Aldehyde/Acyl Halide C=O"
  ("HCO" 218))) ;; "Aldehyde/Formamide HCO-"

("acetaldehyde, nist"
 (("CT3"  (0.6107  0.8391  0.7069))
  ("C"  (2.0250  1.0943  1.1318))
  ("O"  (2.5872  2.1537  0.9742))
  ("HC"  (0.0000  0.5647  1.5771))
  ("HC"  (0.1354  1.7029  0.2243))
  ("HC"  (0.5749  0.0000  0.0000))
  ("HCO"  (2.5561  0.2555  1.6109)))
 (force-field (oplsaa)
  ("C" 216)
  ("O" 217)
  ("HCO" 218)
  ("CT3" 77)
  ("HC" 82)))

("propanal, nist"
 (("C"  (0.5931  1.2165  0.7591))
  ("O"  (0.0510  1.9862  0.0000))
  ("CT2"  (2.0767  0.9416  0.7333))
  ("CT3"  (2.6313  0.8177  2.1333))
  ("HCO"  (0.0000  0.6703  1.5114))
  ("HC"  (2.2352  0.0000  0.1701))
  ("HC"  (2.6247  1.7223  0.1693))
  ("HC"  (3.6982  0.5613  2.1144))
  ("HC"  (2.5307  1.7572  2.6922))
  ("HC"  (2.1168  0.0369  2.7095)))
 (force-field (oplsaa)
  ("C" 216)
  ("O" 217)
  ("HCO" 218)
  ("CT3" 77)
  ("CT2" 78)
  ("HC" 82)))

("butanal, nist"
 (("CT2"  (1.5108  2.8666  2.1213))
  ("CT3"  (0.9716  2.0212  3.2537))
  ("CT2"  (1.6991  2.0635  0.8454))
  ("C"  (2.8449  1.0876  0.9588))
  ("O"  (2.8392  0.0000  0.4293))
  ("HC"  (2.4701  3.3348  2.4211))
  ("HC"  (0.8190  3.7086  1.9227))
  ("HC"  (1.6532  1.1963  3.5016))
  ("HC"  (0.0000  1.5762  3.0008))
  ("HC"  (0.8329  2.6179  4.1642))
  ("HC"  (1.9254  2.7439  0.0000))
  ("HC"  (0.7556  1.5492  0.5728))
  ("HCO"  (3.7193  1.4109  1.5484)))
 (force-field (oplsaa)
  ("C" 216)
  ("O" 217)
  ("HCO" 218)
  ("CT3" 77)
  ("CT2" 78)
  ("HC" 82)))

("pentanal, nist"
 (("CT2"  (1.8430  1.1054  4.4714))
  ("CT2"  (1.9548  1.7583  3.1066))
  ("CT2"  (1.0505  1.0743  2.0933))
  ("C"  (2.7065  1.7213  5.5404))
  ("CT3"  (1.1556  1.7275  0.7336))
  ("O"  (3.4425  2.6662  5.3758))
  ("HC"  (0.7923  1.1357  4.8258))
  ("HC"  (2.1037  0.0294  4.4033))
  ("HC"  (1.6990  2.8343  3.1785))
  ("HC"  (3.0065  1.7280  2.7573))
  ("HC"  (1.3125  0.0000  2.0179))
  ("HC"  (0.0000  1.0998  2.4457))
  ("HCO"  (2.6391  1.2518  6.5367))
  ("HC"  (0.8551  2.7832  0.7673))
  ("HC"  (2.1828  1.6938  0.3468))
  ("HC"  (0.5115  1.2260  0.0000)))
 (force-field (oplsaa)
  ("C" 216)
  ("O" 217)
  ("HCO" 218)
  ("CT3" 77)
  ("CT2" 78)
  ("HC" 82)))

("ammonia, nist"
 (("NT"  (0.8726  0.5133  0.0000))
  ("H"  (0.0000  0.0000  0.0000))
  ("H"  (0.8726  1.4885  0.0000))
  ("H"  (1.7452  0.0000  0.0000)))
 (force-field (oplsaa)
  ("NT" 69)
  ;; ("H" 70)))
  ("H" 1003)))

("hydrogen sulfide, nist"
 (("SH"  (0.8727  0.5134  0.0000))
  ("HS"  (0.0000  0.0513  0.0000))
  ("HS"  (1.7454  0.0000  0.0000)))
 (force-field (oplsaa)
  ("SH" 140)
  ;; ("HS" 144)))
  ("HS" 1004)))

("hydrogen cyanide, nist"
  (("CZ"  (1.0256  0.0000  0.0000))
   ("HC"  (0.0000  0.0000  0.0000))
   ("NZ"  (2.0000  0.0000  0.0000)))
  (force-field (oplsaa)
   ("CZ" 692)
   ("NZ" 691)
   ("HC" 788)))

("methylamine, nist"
  (("HC"  (0.0000  1.0216  0.0000))
   ("CT"  (1.0216  1.0216  0.0000))
   ("HC"  (1.0216  1.9922  0.0000))
   ("HC"  (1.0216  0.0000  0.0000))
   ("NT"  (1.9922  1.0216  0.0000))
   ("H"  (2.5030  1.8389  0.0000))
   ("H"  (2.5030  0.1022  0.0000)))
  (force-field (oplsaa)
   ("NT" 730)
   ("CT" 733)
   ;; ("H" 739)
   ("H" 1005)
   ("HC" 82)))

;;
;; Solutes from [mobley09]. Except:
;;
;; 1. Geometries are from the charge mol2 files, except "fluoroethane"
;; and "111_trichloroethane, mobley12"
;;
;; 2. For "thioephene", "12_dichloroethane", "112_trichloroethane",
;; ff parameters and charges are assigned following the database in
;; http://virtualchemistry.org/gmld.php, also see the change in
;; force-fields.scm
;;
("methylamine"
 (("CTN"  (-0.7670  -0.6855  0.3708))
  ("NT"  (-1.8451  -0.8015  -0.5949))
  ("HC"  (-1.1443  -0.3070  1.3251))
  ("HC"  (0.0000  0.0000  0.0000))
  ("HC"  (-0.3040  -1.6619  0.5391))
  ("H"  (-2.2872  0.1083  -0.7200))
  ("H"  (-2.5666  -1.4187  -0.2245)))
 (force-field (oplsaa)
  ("NT" 730)
  ("CTN" 733)
  ("H" 739)
  ;; ("H" 1005)
  ("HC" 82)))

("ethylamine"
 (("CT3"  (-0.2605  0.7905  0.7119))
  ("CTN"  (0.0815  0.3848  2.1354))
  ("NT"  (1.5069  0.1164  2.2786))
  ("HC"  (0.2716  1.7026  0.4212))
  ("HC"  (0.0000  0.0000  -0.0000))
  ("HC"  (-1.3344  0.9839  0.6226))
  ("HC"  (-0.4834  -0.5110  2.4131))
  ("HC"  (-0.2063  1.1819  2.8287))
  ("H"  (2.0104  0.5793  1.5253))
  ("H"  (1.6713  -0.8810  2.1640)))
 (force-field (oplsaa)
  ("NT" 730)
  ("CT3" 77)
  ("CTN" 736)
  ("H" 739)
  ;; ("H" 1005)
  ("HC" 82)))

("n_propylamine"
 (("CT3"  (0.9952  -0.2871  0.3545))
  ("CT2"  (1.1802  0.0831  1.8180))
  ("CTN"  (0.1575  -0.6240  2.7038))
  ("NT"  (0.3420  -0.2375  4.0975))
  ("HC"  (1.1160  -1.3646  0.2041))
  ("HC"  (0.0000  -0.0000  0.0000))
  ("HC"  (1.7375  0.2276  -0.2635))
  ("HC"  (1.0845  1.1701  1.9267))
  ("HC"  (2.1976  -0.1851  2.1280))
  ("HC"  (0.2626  -1.7106  2.6114))
  ("HC"  (-0.8586  -0.3596  2.3917))
  ("H"  (0.9698  0.5617  4.1408))
  ("H"  (0.8051  -0.9948  4.5942)))
 (force-field (oplsaa)
  ("NT" 730)
  ("CT3" 77)
  ("CT2" 78)
  ("CTN" 736)
  ("H" 739)
  ;; ("H" 1005)
  ("HC" 82)))

("n_butylamine"
 (("CT3"  (1.0009  -0.0195  -0.3646))
  ("CT2"  (1.4148  1.2879  -1.0210))
  ("CT2"  (0.4956  1.6440  -2.1900))
  ("CTN"  (0.8962  2.9762  -2.8215))
  ("NT"  (0.0214  3.2872  -3.9457))
  ("HC"  (1.6717  -0.2548  0.4677))
  ("HC"  (-0.0188  0.0444  0.0283))
  ("HC"  (1.0429  -0.8478  -1.0792))
  ("HC"  (1.3928  2.0887  -0.2727))
  ("HC"  (2.4492  1.2014  -1.3736))
  ("HC"  (0.5384  0.8478  -2.9441))
  ("HC"  (-0.5416  1.6967  -1.8354))
  ("HC"  (0.8346  3.7806  -2.0804))
  ("HC"  (1.9309  2.9285  -3.1772))
  ("H"  (-0.2062  4.2784  -3.9259))
  ("H"  (-0.8572  2.7876  -3.8315)))
 (force-field (oplsaa)
  ("NT" 730)
  ("CT3" 77)
  ("CT2" 78)
  ("CTN" 736)
  ("H" 739)
  ;; ("H" 1005)
  ("HC" 82)))

("n_pentylamine"
 (("CT3"  (0.0375  -0.7160  0.8641))
  ("CT2"  (1.2741  -0.4303  1.7013))
  ("CT2"  (2.5561  -0.6357  0.8947))
  ("CT2"  (3.7991  -0.3799  1.7498))
  ("CTN"  (5.0775  -0.5515  0.9310))
  ("NT"  (6.2493  -0.3237  1.7683))
  ("HC"  (-0.8676  -0.5625  1.4602))
  ("HC"  (-0.0107  -0.0510  -0.0041))
  ("HC"  (0.0368  -1.7503  0.5056))
  ("HC"  (1.2754  -1.0894  2.5773))
  ("HC"  (1.2277  0.6007  2.0709))
  ("HC"  (2.5590  0.0403  0.0310))
  ("HC"  (2.5838  -1.6600  0.5034))
  ("HC"  (3.8069  -1.0736  2.6000))
  ("HC"  (3.7548  0.6360  2.1626))
  ("HC"  (5.0881  0.1533  0.0924))
  ("HC"  (5.1265  -1.5642  0.5167))
  ("H"  (5.9581  0.1262  2.6328))
  ("H"  (6.8715  0.3299  1.2988)))
 (force-field (oplsaa)
  ("NT" 730)
  ("CT3" 77)
  ("CT2" 78)
  ("CTN" 736)
  ("H" 739)
  ;; ("H" 1005)
  ("HC" 82)))

("methanethiol"
 (("CTS"  (0.2788  0.8226  0.6632))
  ("SH"  (-0.4279  2.3668  0.0544))
  ("HC"  (1.3683  0.9000  0.6965))
  ("HC"  (-0.0000  0.0000  0.0000))
  ("HC"  (-0.0956  0.6140  1.6684))
  ("HS"  (0.1845  2.3556  -1.1387)))
 (force-field (oplsaa)
  ("CTS" 156) ;; "Methanethiol CH3-SH"
  ("SH" 139)
  ("HS" 143)
  ("HC" 82)))

("ethanethiol"
 (("CT3"  (0.4412  -0.3496  0.9395))
  ("CTS"  (1.8518  -0.8681  0.7340))
  ("SH"  (2.9258  0.4446  0.0921))
  ("HC"  (0.4164  0.4722  1.6632))
  ("HC"  (0.0000  0.0000  0.0000))
  ("HC"  (-0.1973  -1.1505  1.3271))
  ("HC"  (1.8501  -1.7008  0.0240))
  ("HC"  (2.2645  -1.2294  1.6807))
  ("HS"  (3.5904  -0.3266  -0.7793)))
 (force-field (oplsaa)
  ("CT3" 77)
  ("CTS" 145) ;; "Thiol -CH2-SH"
  ("SH" 139)
  ("HS" 143)
  ("HC" 82)))

("n_butanethiol"
 (("CT3"  (0.8879  -0.0985  -0.6256))
  ("CT2"  (0.8885  -1.4198  -1.3796))
  ("CT2"  (-0.3384  -1.5430  -2.2871))
  ("CTS"  (-0.4452  -2.8933  -2.9954))
  ("SH"  (0.9538  -3.2248  -4.1064))
  ("HC"  (1.7691  -0.0301  0.0200))
  ("HC"  (-0.0027  -0.0060  0.0044))
  ("HC"  (0.9082  0.7482  -1.3192))
  ("HC"  (0.9009  -2.2445  -0.6576))
  ("HC"  (1.8112  -1.4872  -1.9659))
  ("HC"  (-0.3459  -0.7281  -3.0219))
  ("HC"  (-1.2397  -1.4132  -1.6745))
  ("HC"  (-0.5106  -3.7055  -2.2644))
  ("HC"  (-1.3592  -2.9163  -3.5978))
  ("HS"  (1.1284  -4.5003  -3.7355)))
 (force-field (oplsaa)
  ("CT3" 77)
  ("CT2" 78)
  ("CTS" 145)
  ("SH" 139)
  ("HS" 143)
  ("HC" 82)))

("n_propanethiol"
 (("CT3"  (0.1575  -0.4974  0.9622))
  ("CT2"  (1.6397  -0.5220  1.3061))
  ("CTS"  (1.8811  -1.2098  2.6451))
  ("SH"  (3.6533  -1.2345  3.0372))
  ("HC"  (-0.4127  0.0456  1.7227))
  ("HC"  (0.0000  0.0000  -0.0000))
  ("HC"  (-0.2452  -1.5128  0.8897))
  ("HC"  (2.0115  0.5095  1.3338))
  ("HC"  (2.1779  -1.0441  0.5058))
  ("HC"  (1.3554  -0.6852  3.4494))
  ("HC"  (1.5173  -2.2420  2.6204))
  ("HS"  (4.0131  -2.0270  2.0185)))
 (force-field (oplsaa)
  ("CT3" 77)
  ("CT2" 78)
  ("CTS" 145)
  ("SH" 139)
  ("HS" 143)
  ("HC" 82)))

("thiophenol"
 (("CA"  (1.8846  -1.0376  -0.1090))
  ("CA"  (0.8176  -0.4410  0.5626))
  ("CA"  (2.9348  -1.6038  0.6137))
  ("CA"  (0.8008  -0.4105  1.9570))
  ("CA"  (2.9180  -1.5734  2.0082))
  ("CAS"  (1.8509  -0.9767  2.6799))
  ("SH"  (1.8293  -0.9375  4.4617))
  ("HA"  (1.8979  -1.0618  -1.1946))
  ("HA"  (0.0005  0.0004  -0.0006))
  ("HA"  (3.7653  -2.0685  0.0904))
  ("HA"  (-0.0393  0.0592  2.4630))
  ("HA"  (3.7455  -2.0197  2.5544))
  ("HS"  (0.6757  -0.2675  4.6037)))
 (force-field (oplsaa)
  ("CA" 87)
  ("HA" 88)
  ("CAS" 673)
  ("SH" 672)
  ("HS" 143)))

("benzene"
 (("CA"  (1.8846  -1.0376  -0.1090))
  ("CA"  (0.8176  -0.4410  0.5626))
  ("CA"  (2.9348  -1.6038  0.6137))
  ("CA"  (0.8008  -0.4105  1.9570))
  ("CA"  (2.9180  -1.5734  2.0082))
  ("CA"  (1.8509  -0.9767  2.6799))
  ("HA"  (1.8980  -1.0609  -1.1947))
  ("HA"  (0.0001  0.0000  -0.0001))
  ("HA"  (3.7660  -2.0676  0.0907))
  ("HA"  (-0.0298  0.0544  2.4800))
  ("HA"  (3.7360  -2.0137  2.5709))
  ("HA"  (1.8378  -0.9527  3.7656)))
 (force-field (oplsaa)
  ("CA" 87)
  ("HA" 88)))

("pyridine"
 (("CA3"  (1.2858  -0.0601  1.7333))
  ("CA2"  (0.1842  -0.4388  0.9739))
  ("CA2"  (1.4967  -0.6399  2.9793))
  ("CA1"  (-0.6697  -1.3939  1.5014)) ;; the two C atoms next to N are most positive
  ("CA1"  (0.5838  -1.5858  3.4172))
  ("NC"  (-0.4976  -1.9792  2.7078))
  ("HA3"  (1.9788  0.6851  1.3547))
  ("HA2"  (0.0002  -0.0002  -0.0005))
  ("HA2"  (2.3487  -0.3604  3.5885))
  ("HA1"  (-1.5433  -1.7237  0.9484))
  ("HA1"  (0.7027  -2.0679  4.3821)))
 (force-field (oplsaa)
  ("NC" 458)
  ("CA1" 459)
  ("CA2" 460)
  ("CA3" 461)
  ("HA1" 462)
  ("HA2" 463)
  ("HA3" 464)))

("benzyl_alcohol"
 (("CA"  (1.8846  -1.0376  -0.1090))
  ("CA"  (0.8176  -0.4410  0.5626))
  ("CA"  (2.9348  -1.6038  0.6137))
  ("CA"  (0.8008  -0.4105  1.9570))
  ("CA"  (2.9180  -1.5734  2.0082))
  ("CAOH"  (1.8509  -0.9767  2.6799))
  ("CT"  (1.8328  -0.9437  4.1709))
  ("OH"  (2.4730  0.2424  4.6249))
  ("HA"  (1.8979  -1.0612  -1.1948))
  ("HA"  (0.0002  0.0000  -0.0002))
  ("HA"  (3.7655  -2.0682  0.0906))
  ("HA"  (-0.0347  0.0567  2.4718))
  ("HA"  (3.7404  -2.0171  2.5631))
  ("HC"  (0.8054  -0.9589  4.5597))
  ("HC"  (2.3572  -1.8117  4.5934))
  ("HO"  (2.4691  0.8644  3.8798)))
 (force-field (oplsaa)
  ("OH" 93)
  ("HO" 94)
  ("CT" 157) ;; "Benzyl Alcohol -CH2OH"
  ("CAOH" 160) ;; "Benzyl Alcohol/Nitrile" connected to -CH2OH
  ("CA" 87)  ;; Aromatic C"
  ("HA" 88)
  ("HC" 82)))

("benzaldehyde"
 (("CA"  (1.8846  -1.0376  -0.1090))
  ("CA"  (0.8176  -0.4410  0.5626))
  ("CA"  (2.9348  -1.6038  0.6137))
  ("CA"  (0.8008  -0.4105  1.9570))
  ("CA"  (2.9180  -1.5734  2.0082))
  ("CA"  (1.8509  -0.9767  2.6799))
  ("C"  (1.8333  -0.9446  4.1336))
  ("O"  (2.7308  -1.4281  4.8293))
  ("HA"  (1.8975  -1.0613  -1.1948))
  ("HA"  (0.0001  0.0000  -0.0000))
  ("HA"  (3.7651  -2.0690  0.0906))
  ("HA"  (-0.0358  0.0574  2.4697))
  ("HA"  (3.7417  -2.0179  2.5610))
  ("HC"  (0.9622  -0.4579  4.6242)))
 (force-field (oplsaa)
  ("C" 171)
  ("O" 217)
  ("CA" 87)
  ("HA" 88)
  ("HC" 218)))

;;
;; See http://virtualchemistry.org/molecules/110-02-1/OPLS/110-02-1.top
;; There parameters for "Furan" are used, except for "S", only the
;; atomic charge of "Furan O" is applied, LJ parameters are from normal
;; "S" in the database
;;
("thiophene"
 (("CS"  (-0.0619  0.6736  -2.1200))
  ("CS"  (-0.3790  -0.1513  -1.0023))
  ("CW"  (-0.6859  0.2469  -3.2715))
  ("CW"  (-1.2335  -1.1769  -1.3423))
  ("S"  (-1.6496  -1.1428  -3.0031))
  ("HAS"  (0.5911  1.5356  -2.0810))
  ("HAS"  (-0.0009  0.0012  0.0001))
  ("HAW"  (-0.6210  0.6847  -4.2576))
  ("HAW"  (-1.6313  -1.9457  -0.6951)))
 (force-field (oplsaa)
  ("S" 1008) ;; pseudo thiophene S
  ("CW" 505)
  ("CS" 506)
  ("HAW" 507)
  ("HAS" 508)))

("phenol"
 (("CA"  (1.8846  -1.0376  -0.1090))
  ("CA"  (0.8176  -0.4410  0.5626))
  ("CA"  (2.9348  -1.6038  0.6137))
  ("CA"  (0.8008  -0.4105  1.9570))
  ("CA"  (2.9180  -1.5734  2.0082))
  ("CAO"  (1.8509  -0.9767  2.6799))
  ("OH"  (1.8342  -0.9476  4.0403))
  ("HA"  (1.8990  -1.0598  -1.1948))
  ("HA"  (-0.0004  -0.0007  -0.0000))
  ("HA"  (3.7656  -2.0682  0.0905))
  ("HA"  (-0.0336  0.0549  2.4746))
  ("HA"  (3.7386  -2.0165  2.5657))
  ("HO"  (1.1790  -0.2881  4.3334)))
 (force-field (oplsaa)
  ("CAO" 105)
  ("OH" 106)
  ("HO" 107)
  ("CA" 87)
  ("HA" 88)))

("chlorobenzene"
 (("CA"  (1.8846  -1.0376  -0.1090))
  ("CA"  (0.8176  -0.4410  0.5626))
  ("CA"  (2.9348  -1.6038  0.6137))
  ("CA"  (0.8008  -0.4105  1.9570))
  ("CA"  (2.9180  -1.5734  2.0082))
  ("CACl"  (1.8509  -0.9767  2.6799))
  ("Cl"  (1.8302  -0.9386  4.4007))
  ("HA"  (1.8979  -1.0610  -1.1947))
  ("HA"  (0.0001  0.0000  -0.0001))
  ("HA"  (3.7658  -2.0678  0.0907))
  ("HA"  (-0.0373  0.0583  2.4666))
  ("HA"  (3.7436  -2.0183  2.5579)))
 (force-field (oplsaa)
  ("CACl" 202)
  ("Cl" 203)
  ("CA" 87)
  ("HA" 88)))

("2_chlorophenol"
 (("CA"  (1.8846  -1.0376  -0.1090))
  ("CA"  (0.8176  -0.4410  0.5626))
  ("CA"  (2.9348  -1.6038  0.6137))
  ("CA"  (0.8008  -0.4105  1.9570))
  ("CAO"  (2.9180  -1.5734  2.0082))
  ("CACl"  (1.8509  -0.9767  2.6799))
  ("OH"  (3.9570  -2.1349  2.6890))
  ("Cl"  (1.8053  -0.9221  4.4035))
  ("HA"  (1.8989  -1.0591  -1.1949))
  ("HA"  (0.0007  0.0019  0.0003))
  ("HA"  (3.7625  -2.0670  0.0834))
  ("HA"  (-0.0372  0.0588  2.4661))
  ("HO"  (3.8734  -1.9056  3.6319)))
 (force-field (oplsaa)
  ("CACl" 202)
  ("CAO" 105)
  ("OH" 106)
  ("HO" 107)
  ("Cl" 203)
  ("CA" 87)
  ("HA" 88)))

("3_chlorophenol"
 (("CA"  (0.8324  -0.3948  0.5751))
  ("CA"  (1.8991  -1.0135  -0.0768))
  ("CA"  (0.8349  -0.2829  1.9655))
  ("CA"  (2.9708  -1.4085  2.0522))
  ("CAO"  (2.9683  -1.5204  0.6618))
  ("CACl"  (1.9042  -0.7896  2.7040))
  ("OH"  (4.0086  -2.1248  0.0259))
  ("Cl"  (1.9058  -0.6518  4.4192))
  ("HA"  (0.0000  0.0002  0.0001))
  ("HA"  (1.8917  -1.0973  -1.1601))
  ("HA"  (-0.0029  0.2015  2.4604))
  ("HA"  (3.8079  -1.8070  2.6207))
  ("HO"  (4.6596  -2.4137  0.6913)))
 (force-field (oplsaa)
  ("CACl" 202)
  ("CAO" 105)
  ("OH" 106)
  ("HO" 107)
  ("Cl" 203)
  ("CA" 87)
  ("HA" 88)))

("4_chlorophenol"
 (("CA"  (0.8178  -0.4385  0.5644))
  ("CA"  (2.9371  -1.5972  0.6197))
  ("CA"  (0.7994  -0.4058  1.9588))
  ("CA"  (2.9185  -1.5646  2.0141))
  ("CAO"  (1.8867  -1.0341  -0.1051))
  ("CACl"  (1.8497  -0.9689  2.6836))
  ("OH"  (1.9059  -1.0642  -1.4656))
  ("Cl"  (1.8267  -0.9286  4.4038))
  ("HA"  (-0.0027  0.0029  0.0053))
  ("HA"  (3.7728  -2.0622  0.1040))
  ("HA"  (-0.0400  0.0617  2.4675))
  ("HA"  (3.7440  -2.0068  2.5662))
  ("HO"  (2.6776  -1.5828  -1.7582)))
 (force-field (oplsaa)
  ("CACl" 202)
  ("CAO" 105)
  ("OH" 106)
  ("HO" 107)
  ("Cl" 203)
  ("CA" 87)
  ("HA" 88)))

("chloroethane"
 (("CT3"  (0.1543  -0.7247  0.8062))
  ("CT"  (1.3672  -0.3680  1.6395))
  ("Cl"  (2.8321  -0.3570  0.6342))
  ("HC"  (-0.7442  -0.7295  1.4311))
  ("HC"  (0.0000  -0.0000  0.0000))
  ("HC"  (0.2588  -1.7188  0.3591))
  ("HCCl"  (1.5144  -1.0958  2.4424))
  ("HCCl"  (1.2550  0.6252  2.0832)))
 (force-field (oplsaa)
  ("CT3" 77)
  ("HC" 82)
  ("Cl" 800)
  ("CT" 801)
  ("HCCl" 802))) ;; "Alkyl Chloride H-C-Cl"

("bromoethane"
 (("CT3"  (0.0125  -0.7247  0.8210))
  ("CT"  (1.1399  -0.4427  1.7865))
  ("Br"  (2.8461  -0.5455  0.8594))
  ("HC"  (-0.9510  -0.6676  1.3370))
  ("HC"  (0.0000  -0.0000  0.0000))
  ("HC"  (0.1031  -1.7247  0.3836))
  ("HCBr"  (1.1722  -1.1779  2.5944))
  ("HCBr"  (1.0673  0.5629  2.2082)))
 (force-field (oplsaa)
  ("CT3" 77)
  ("HC" 82)
  ("Br" 805)
  ("CT" 806)
  ("HCBr" 807)))

("iodoethane"
 (("CT3"  (0.1447  -0.6955  0.8338))
  ("CT"  (1.4101  -0.3747  1.6104))
  ("I"  (3.1275  -0.5146  0.3556))
  ("HC"  (-0.7286  -0.6194  1.4899))
  ("HC"  (0.0000  -0.0000  0.0000))
  ("HC"  (0.1715  -1.7122  0.4272))
  ("HCI"  (1.5088  -1.0687  2.4520))
  ("HCI"  (1.3374  0.6351  2.0281)))
 (force-field (oplsaa)
  ("CT3" 77)
  ("HC" 82)
  ("I" 838)
  ("CT" 835)
  ("HCI" 839)))

;;
;; geometry from NIST webbook
;;
("fluoroethane, raw"
 (("CT3"  (-0.0002  9.9997  0.0001))
  ("CT"  (1.3026  9.2236  -0.0001))
  ("F"  (2.3785  10.1081  0.0001))
  ("HC"  (-0.8544  9.3128  -0.0001))
  ("HC"  (-0.0683  10.6366  0.8880))
  ("HC"  (-0.0683  10.6370  -0.8875))
  ("HCF"  (1.3873  8.5867  -0.8899))
  ("HCF"  (1.3873  8.5862  0.8894)))
 (force-field (oplsaa)
  ("CT3" 77)
  ("HC" 82)
  ("F" 786)
  ("CT" 787)
  ("HCF" 788)))

;;
;; z coordinates are shifted to make CT3 at the center
;;
("fluoroethane"
 (("CT3"  (0.0000  0.0000  0.0000))
  ("CT"  (1.3030  -0.7760  -0.0000))
  ("F"  (2.3790  0.1080  0.0000))
  ("HC"  (-0.8540  -0.6870  -0.0000))
  ("HC"  (-0.0680  0.6370  0.8880))
  ("HC"  (-0.0680  0.6370  -0.8880))
  ("HCF"  (1.3880  -1.4130  -0.8900))
  ("HCF"  (1.3880  -1.4140  0.8890)))
 (force-field (oplsaa)
  ("CT3" 77)
  ("HC" 82)
  ("F" 786)
  ("CT" 787)
  ("HCF" 788)))


;;
;; charge of Cl is the same as
;; http://virtualchemistry.org/molecules/75-34-3/OPLS/75-34-3.top
;;
("11_dichloroethane"
 (("CT3"  (-0.4712  0.6837  0.7145))
  ("CTCl"  (-0.0877  0.3252  2.1389))
  ("Cl"  (-0.8870  1.4225  3.2882))
  ("Cl"  (1.6769  0.4133  2.3441))
  ("HC"  (-0.1656  1.7038  0.4555))
  ("HC"  (-0.0000  0.0000  0.0000))
  ("HC"  (-1.5551  0.6111  0.5727))
  ("HCCl"  (-0.4016  -0.6948  2.3770)))
 (force-field (oplsaa)
  ("CT3" 77)
  ("HC" 82)
  ("Cl" 1009)
  ("CTCl" 803)
  ("HCCl" 802)))

("12_dichloroethane"
 (("CT"  (-0.3490  0.4174  -0.9492))
  ("CT"  (-1.3747  -0.5098  -1.5787))
  ("Cl"  (1.0832  0.6225  -1.9863))
  ("Cl"  (-1.9766  0.1239  -3.1295))
  ("HC"  (-0.7693  1.4085  -0.7518))
  ("HC"  (0.0000  0.0000  0.0000))
  ("HC"  (-2.2403  -0.6143  -0.9179))
  ("HC"  (-0.9602  -1.5067  -1.7579)))
 (force-field (oplsaa)
  ("Cl" 800)
  ("CT" 801)
  ("HC" 802)))

;;
;; No suitable charges for -CCl3, these are from charged mol2 file in [mobley09]
;;
("111_trichloroethane"
 (("CT3"  (1.0245  0.2811  -0.2703) 3.5 0.066 -0.1089)
  ("CTCl"  (1.2871  1.7565  0.0012) 3.5 0.066 0.2656)
  ("Cl"  (0.1596  2.7554  -0.9573) 3.4 0.3 -0.1242)
  ("Cl"  (1.0484  2.1022  1.7366) 3.4 0.3 -0.1242)
  ("Cl"  (2.9647  2.1642  -0.4534) 3.4 0.3 -0.1242)
  ("HC"  (1.7030  -0.3597  0.3051) 2.5 0.03 0.0720)
  ("HC"  (0.0000  -0.0000  0.0000) 2.5 0.03 0.0720)
  ("HC"  (1.1629  0.0376  -1.3305) 2.5 0.03 0.0720)))

;;
;; These are from [mobley12]http://dx.doi.org/10.1007/s10822-011-9528-8
;;
("111_trichloroethane, mobley12"
(("C1" ( 0.0000   -0.0000    1.7630) (A 3.39967) (kj 4.57730e-01) -0.491496)
 ("C2" (-0.0000   -0.0000    0.2540) (A 3.39967) (kj 4.57730e-01) -0.268343)
 ("Cl1"(-1.6470   -0.2770   -0.3620) (A 3.47094) (kj 1.10876e+00)  0.028497)
 ("Cl2"( 0.5840    1.5640   -0.3620) (A 3.47094) (kj 1.10876e+00)  0.028497)
 ("Cl3"( 1.0630   -1.2880   -0.3620) (A 3.47094) (kj 1.10876e+00)  0.028497)
 ("H1" ( 1.0140    0.1730    2.1180) (A 2.64953) (kj 6.56888e-02)  0.224783)
 ("H2" (-0.6570    0.7910    2.1190) (A 2.64953) (kj 6.56888e-02)  0.224783)
 ("H3" (-0.3570   -0.9650    2.1190) (A 2.64953) (kj 6.56888e-02)  0.224783)))


;;
;; charge of Cl on -CHCl2 is the same as
;; http://virtualchemistry.org/molecules/79-00-5/OPLS/79-00-5.top
;;
("112_trichloroethane"
 (("CT1"  (-0.3098  -0.9813  0.3742))
  ("CT2"  (0.3109  -2.0958  -0.4640))
  ("Cl1"  (-2.0927  -1.0020  0.3817))
  ("Cl2"  (-0.1403  -1.9651  -2.1822))
  ("Cl2"  (-0.1100  -3.7101  0.1605))
  ("HC"  (0.0182  -1.0557  1.4163))
  ("HC"  (-0.0000  0.0000  0.0000))
  ("HC"  (1.4010  -2.0167  -0.4190)))
 (force-field (oplsaa)
  ("CT1" 801)
  ("CT2" 803)
  ("HC" 802)
  ("Cl1" 800)
  ("Cl2" 1009)))


;;
;; Ions
;;
;; In [CFCRS08], Chuev et al. reported results of RISM-HNC/KH+RBC for
;; ions with OPLS force field. However, after checking the cited paper
;; [JUT04], LJ parameters in "oplsaa.prm" are different with those cited
;; by Chuev, parameters in another file "oplsaal.prm", seem to be more
;; close to those in [JUT04]. Note that the original OPLSAA paper
;; doesn't contain force field for alkali and halide ions and there is
;; no obvious clue indicating the source for those FF parameters of ions
;; in "oplsaa.prm"
;;
;; [CFCRS08] Hydration of ionic species studied by the reference
;;   interaction site model with a repulsive bridge correction, Gennady
;;   N.  Chuev, Maxim V. Fedorov, Sandro Chiodo, Nino Russo, Emilia
;;   Sicilia, Journal of Computational Chemistry, Volume 29, Issue 14,
;;   pages 2406–2415, 15 November 2008
;;   http://dx.doi.org/10.1002/jcc.20979
;;
;; [JUT04] Free Energies of Hydration from a Generalized Born Model and
;;   an All-Atom Force Field, W. L. Jorgensen, J. P. Ulmschneider and J.
;;   Tirado-Rives, The Journal of Physical Chemistry B, 2004, 104 (41),
;;   pp 16264-16270
;;   http://dx.doi.org/10.1021/jp0484579
;;
("Na+, oplsaa" (("Na+" (0.0 0.0 0.0))) (force-field (oplsaa) ("Na+" 346)))

("Cl-, oplsaa" (("Cl-" (0.0 0.0 0.0))) (force-field (oplsaa) ("Cl-" 341)))

;;
;; Monoatomic ions cited by [CFCRS08], parameters could be found in
;; [JUT04]. As far as I have observed, bivalent ions here have almost
;; the same parameters with those ./oplsaa-test table, only less
;; digits, but there are more significant difference in monovalent ions
;;
("Li+, JUT04" (("Li+" (0.0 0.0 0.0))) (force-field (jut04) ("Li+" Li)))
("Na+, JUT04" (("Na+" (0.0 0.0 0.0))) (force-field (jut04) ("Na+" Na)))
("K+, JUT04" (("K+" (0.0 0.0 0.0))) (force-field (jut04) ("K+" K)))
("Rb+, JUT04" (("Rb+" (0.0 0.0 0.0))) (force-field (jut04) ("Rb+" Rb)))
("Cs+, JUT04" (("Cs+" (0.0 0.0 0.0))) (force-field (jut04) ("Cs+" Cs)))

("Mg2+, JUT04" (("Mg2+" (0.0 0.0 0.0))) (force-field (jut04) ("Mg2+" Mg)))
("Ca2+, JUT04" (("Ca2+" (0.0 0.0 0.0))) (force-field (jut04) ("Ca2+" Ca)))
("Sr2+, JUT04" (("Sr2+" (0.0 0.0 0.0))) (force-field (jut04) ("Sr2+" Sr)))
("Ba2+, JUT04" (("Ba2+" (0.0 0.0 0.0))) (force-field (jut04) ("Ba2+" Ba)))

("F-, JUT04" (("F-" (0.0 0.0 0.0))) (force-field (jut04) ("F-" F)))
("Cl-, JUT04" (("Cl-" (0.0 0.0 0.0))) (force-field (jut04) ("Cl-" Cl)))
("Br-, JUT04" (("Br-" (0.0 0.0 0.0))) (force-field (jut04) ("Br-" Br)))
("I-, JUT04" (("I-" (0.0 0.0 0.0))) (force-field (jut04) ("I-" I)))

;;
;; CFCRS08 polyatomic ions:
;;
;; Several polyatomic ions which are studied in [CFCRS08] (see above).
;; Determined structures which could be found in NIST databse are for
;; OH-, CHOO-, CN-, NH4+.
;;
;; For CH3COO-, since there is only a 2D-mol file in NIST database in
;; which hydrogen atoms are missing, the structure here is optimized in
;; PG by removing H of -COOH and set the molecule with -1 charge,
;; minimal basis set (crenbl_ecp + ahlrichs coulomb fitting) is used.
;;
;; For H3O+, in the 2D-mol file ("hydronium cation-2d") of H3O+ from
;; NIST database all the atoms are in a plane and in the 3D-sdf
;; ("hydronium cation-3d") file Ow-H distance is larger than 3.0 A, but
;; according to Ion, H3O+ has a trigonal pyrameid geometry and the
;; strucutre he used is also included here ("hydronium cation"). FIXME:
;; using the minimal basis set I still get a nearly planar structure
;; after geometry optimization.
;;
;; It should be noticed that in [CFCRS08], atomic charges are derived
;; from QM calculation while LJ parameters are basically from OPLS force
;; field. Since there is no entry for H3O+ in OPLSAA, I think perhaps we
;; have to derive atomic charge for polyatomic ions by ourselves.
;;
("hydroxyl anion, CFCRS08"
  (("OW"  (-0.0203  0.0000  0.0000))
   ("HW"  (0.9623  0.0000  0.0000)))
  (force-field (cfcrs08)
   ("HW" 3)
   ("OW" 4)))

("formate, CFCRS08"
 (("O2"  (1.9226  0.9205  1.8422))
  ("C"  (0.9364  0.3495  1.1289))
  ("O2"  (1.2026  0.0000  0.0000))
  ("HC" (0.0000  0.2800  1.6921)))
 (force-field (cfcrs08)
   ("O2" 10)
   ("C"  9)
   ("HC" 11)))

("cyanide ion, CFCRS08"
  (("C"  (1.0256  0.0000  0.0000))
   ("N3"  (2.0000  0.0000  0.0000)))
  (force-field (cfcrs08)
    ("C" 7)
    ("N3" 8)))

("ammonium cation, CFCRS08"
  (("N3"  (0.0000  0.0000  -0.0001))
   ("H"  (1.0847  -0.0002  0.0270))
   ("H"  (-0.3619  1.0226  0.0267))
   ("H"  (-0.3395  -0.4804  -0.9124))
   ("H"  (-0.3832  -0.5417  0.8585)))
  (force-field (cfcrs08)
    ("N3" 6)
    ("H" 5)))

;; Optimized in PG
("acetate, CFCRS08"
 (("C"  (2.1371  1.7060  0.8952))
  ("O2"  (3.1306  0.9770  0.4943))
  ("O2"  (2.1315  2.9685  1.1820))
  ("CT"  (0.7645  0.9594  1.0071))
  ("HC"  (0.0241  1.5696  1.5617))
  ("HC"  (0.3692  0.7573  -0.0088))
  ("HC"  (0.9032  -0.0151  1.5161)))
 (force-field (cfcrs08)
  ("C" 12)
  ("O2" 13)
  ("CT" 14)
  ("HC" 15)))

;; Provided by Ion
("hydronium cation, CFCRS08"
  (("OW" (0.0958  0.1509  0.0011))
   ("HW" (-0.0258  -0.7738  0.3312))
   ("HW" (-0.3383  0.8303  0.5747))
   ("HW" (1.0366  0.3670  -0.2158)))
  (force-field (cfcrs08)
    ("OW" 2)
    ("HW" 1)))

;; From 2D mol file
("hydronium cation-2d, CFCRS08"
  (("HW"  (0.8664  1.4994  0.0000))
   ("OW"  (0.8664  0.4998  0.0000))
   ("HW"  (0.0000  0.0000  0.0000))
   ("HW"  (1.7327  0.0000  0.0000)))
  (force-field (cfcrs08)
    ("OW" 2)
    ("HW" 1)))

;; From 3D sdf file
("hydronium cation-3d, CFCRS08"
  (("HW"  (1.0931  1.8747  0.0677))
   ("OW"  (0.5910  1.0469  -0.0260))
   ("HW"  (-0.1327  1.1487  0.6157))
   ("HW"  (2.8317  0.6876  -1.4673)))
  (force-field (cfcrs08)
    ("OW" 2)
    ("HW" 1)))

;;;
;;; Force field parameters for these polyatomic ions could also be found
;;; in ./oplsaa-test table. Structures are the same with their "CFRCS08"
;;; counterparts listed above
;;;
("hydroxyl anion, oplsaa"
  (("OH"  (-0.0203  0.0000  0.0000))
   ("HO"  (0.9623  0.0000  0.0000)))
  (force-field (oplsaa)
   ("OH" 373)
   ;; ("HO" 374)))
   ("HO" 1007)))  ;; again, use finite-size hydrogen

("formate, oplsaa"
 (("O"  (1.9226  0.9205  1.8422))
  ("C"  (0.9364  0.3495  1.1289))
  ("O"  (1.2026  0.0000  0.0000))
  ("HCO"  (0.0000  0.2800  1.6921)))
 (force-field (oplsaa)
  ("C" 210)    ;; "Carboxylic Acid -COOH"
  ("O" 211)    ;; "Carboxylic Acid C=O"
  ("HCO" 218)));; FIXME: "Aldehyde/Formamide HCO-" is the only H in the
               ;; table which has +1 charge

("acetate, oplsaa"
 (("C"  (2.1371  1.7060  0.8952))
  ("O"  (3.1306  0.9770  0.4943))
  ("O"  (2.1315  2.9685  1.1820))
  ("CT"  (0.7645  0.9594  1.0071))
  ("HC"  (0.0241  1.5696  1.5617))
  ("HC"  (0.3692  0.7573  -0.0088))
  ("HC"  (0.9032  -0.0151  1.5161)))
 (force-field (oplsaa)
  ("C" 210)
  ("CT" 212)
  ("O" 211)
  ("HC" 82)))

;;;
;;; See  test-qm/uranyl/sym/*.gx for  the original  geometries  and PG
;;; inputs used to optimize them. Converted by:
;;;
;;;   for f in ?w,aq.xyz; do runbgy.scm read-xyz $f; done
;;;
;;; Self-energy, E(self):
;;;
;;;   runbgy.scm self-energy "uranyl, %dw, pcm"
;;;
;;; Complex solvation energy, E(RISM):
;;;
;;;   runbgy.scm energy --norm-tol 1e-14 --dielectric 78.4 --rho
;;;     0.0333295 --beta 1.6889 --L 160 --N 4096 --solvent "water,
;;;     PR-SPC/E" --solute "uranyl, %dw, pcm"
;;;
;;; n    E(self), kcal      E(RISM), kcal
;;; ---  -----------------  -----------------
;;; 0       0.0             -376.027350192114
;;; 4    -232.074402281419  -174.594711875714
;;; 5    -271.654767148776  -159.157414678584
;;; 6    -282.556908726866  -143.250721366546
;;;
("uranyl, 0w, pcm"
 (("U" (0.0 0.0 0.0))
  ("OU" (0.0 0.0 1.738053606507))
  ("OU" (0.0 0.0 -1.738053606507)))
 (force-field (uranyl/kl2 water/pr-spc)
  ("U" U)
  ("OU" OU)
  ("HW" HW)
  ("OW" OW)))

("uranyl, 4w, pcm"                      ; D4H
 (("U" (0.0 -0.0 0.0))
  ("OU" (-0.0 0.0 1.787391536814))
  ("OU" (-0.0 -0.0 -1.787391536814))
  ("OW" (2.372134493745 0.0 0.0))
  ("OW" (0.0 2.372134493745 -0.0))
  ("OW" (0.0 -2.372134493745 -0.0))
  ("OW" (-2.372134493745 -0.0 -0.0))
  ("HW" (2.949484697473 0.0 0.81649648835))
  ("HW" (0.0 2.949484697473 0.81649648835))
  ("HW" (0.0 -2.949484697473 0.81649648835))
  ("HW" (-2.949484697473 -0.0 0.81649648835))
  ("HW" (2.949484697473 0.0 -0.81649648835))
  ("HW" (-2.949484697473 0.0 -0.81649648835))
  ("HW" (0.0 2.949484697473 -0.81649648835))
  ("HW" (0.0 -2.949484697473 -0.81649648835)))
 (force-field (uranyl/kl2 water/pr-spc)
  ("U" U)
  ("OU" OU)
  ("HW" HW)
  ("OW" OW)))

("uranyl, 5w, pcm"                      ; D5H
 (("U" (0.0 0.0 0.0))
  ("OU" (0.0 0.0 1.79070730381))
  ("OU" (-0.0 0.0 -1.79070730381))
  ("OW" (2.434873537581 -0.0 -0.0))
  ("OW" (0.752417302266 2.315702344272 -0.0))
  ("OW" (0.752417302266 -2.315702344272 -0.0))
  ("OW" (-1.969854071057 1.431182756588 0.0))
  ("OW" (-1.969854071057 -1.431182756588 -0.0))
  ("HW" (3.012223741309 -0.0 0.81649648835))
  ("HW" (0.930828326924 2.864795017711 0.81649648835))
  ("HW" (0.930828326924 -2.864795017711 0.81649648835))
  ("HW" (-2.436940197579 1.770540691747 0.81649648835))
  ("HW" (-2.436940197579 -1.770540691747 0.81649648835))
  ("HW" (3.012223741309 -0.0 -0.81649648835))
  ("HW" (-2.436940197579 1.770540691747 -0.81649648835))
  ("HW" (0.930828326924 -2.864795017711 -0.81649648835))
  ("HW" (0.930828326924 2.864795017711 -0.81649648835))
  ("HW" (-2.436940197579 -1.770540691747 -0.81649648835)))
 (force-field (uranyl/kl2 water/pr-spc)
  ("U" U)
  ("OU" OU)
  ("HW" HW)
  ("OW" OW)))

("uranyl, 6w, pcm"                      ; D3D
 (("U" (0.0 0.0 0.0))
  ("OU" (0.0 0.0 1.799325495835))
  ("OU" (0.0 -0.0 -1.799325495835))
  ("OW" (0.0 2.382682819916 0.710764543024))
  ("OW" (-2.063463851208 -1.191341409958 0.710764543024))
  ("OW" (2.063463851208 -1.191341409958 0.710764543024))
  ("OW" (0.0 -2.382682819916 -0.710764543024))
  ("OW" (-2.063463851208 1.191341409958 -0.710764543024))
  ("OW" (2.063463851208 1.191341409958 -0.710764543024))
  ("HW" (0.0 2.711799660033 1.655053629471))
  ("HW" (-2.348487395563 -1.355899830017 1.655053629471))
  ("HW" (2.348487395563 -1.355899830017 1.655053629471))
  ("HW" (0.0 -2.711799660033 -1.655053629471))
  ("HW" (-2.348487395563 1.355899830017 -1.655053629471))
  ("HW" (2.348487395563 1.355899830017 -1.655053629471))
  ("HW" (0.0 3.163261495114 0.08570718161))
  ("HW" (-2.739464813582 -1.581630747558 0.08570718161))
  ("HW" (2.739464813582 -1.581630747557 0.08570718161))
  ("HW" (0.0 -3.163261495114 -0.08570718161))
  ("HW" (-2.739464813582 1.581630747558 -0.08570718161))
  ("HW" (2.739464813582 1.581630747558 -0.08570718161)))
 (force-field (uranyl/kl2 water/pr-spc)
  ("U" U)
  ("OU" OU)
  ("HW" HW)
  ("OW" OW)))

;; Acetonitrile
;; Geometry from NIST database
;; FF from http://dx.doi.org/10.1002/jcc.20721
;; mass = 41.053 u
;; density (298 K) = 776.4 kg/m^3 => 0.01138942364829424567 / A^3
("acetonitrile"
  (("CT"  (0.4219  0.8958  0.5193) (r0 (A 1.908)) (kcal 0.1094) -0.5503)
   ("YC"  (1.8612  0.9422  0.4897) (r0 (A 1.99)) (kcal 0.1341) 0.4917)
   ("YN"  (3.0198  0.9796  0.4659) (r0 (A 1.69)) (kcal 0.1331) -0.5126)
   ("HC"  (0.0000  1.7811  0.0257) (r0 (A 1.487)) (kcal 0.0157) 0.1904)
   ("HC"  (0.0568  0.0000  0.0000) (r0 (A 1.487)) (kcal 0.0157) 0.1904)
   ("HC"  (0.0607  0.8694  1.5558) (r0 (A 1.487)) (kcal 0.0157) 0.1904)))


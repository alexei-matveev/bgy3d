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
("K+" (("K+" (0.0 0.0 0.0) (A 3.14265) (kcal 0.087) 1.0)))
("Cl-" (("Cl-" (0.0 0.0 0.0) (A 4.04468) (kcal 0.15) -1.0)))

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
;; more than -27.6 kJ/mol of the original SPC model [1].
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
;; [1] "The Missing Term in Effective Pair Potentials", H.  J.  C.
;;     Berendsen, J.  R.  Grigera, and T.  P.  Straatsma, J.  Phys.
;;     Chem 1987, 91, 6269-6271. http://dx.doi.org/10.1021/j100308a038
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
("OW" (("OW" (0.0 0.0 0.0) 3.16 0.1549 0.0)))

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
;; By comparing atom by atom, we found that FFs for atom sites in methanol,
;; butanoic acid and hexane which are recorded in the original code are the
;; same with those in OPLS All-Atom (OPLS-AA) force field, see the supporting
;; information of [1]. We added place description of the matching entry in the
;; tables of supporting information of [1] for each type of atom in solutes,
;; hope that one can easily build new solutes (alcohols, carboxylic acids and
;; alkanes) by re-using these entries
;;
;; [1] "Development and Testing of the OPLS All-Atom Force Field on
;;     Conformational Energetics and Properties of Organic Liquids", W.
;;     L. Jorgensen, D. S. Maxwell and J. Tirado-Rives, Journal of the
;;     American Chemical Society, 1996, 118 (45), 11225-11236
;;     http://dx.doi.org/10.1021/ja9621760
("methanol"
 (("C" (-0.748 -0.015 0.024) 3.5 0.066 0.145)  ;; 18th entry of table 1
  ("HC1" (-1.293 -0.202 -0.901) 2.5 0.03 0.04) ;; 17th entry of table 1
  ("HC2" (-1.263 0.754 0.6) 2.5 0.03 0.04)
  ("HC3" (-0.699 -0.934 0.609) 2.5 0.03 0.04)
  ("O" (0.558 0.42 -0.278) 3.12 0.17 -0.683)   ;; 15th entry of table 1
  ("OH" (0.716 1.404 0.137) 0.4 0.04 0.418)))  ;; 16th entry of table 1
;;                                             FIXME: note that LJ parameters
;;                                             for "H" are actually "0.0" in
;;                                             the table, here small values are
;;                                             assigned possibly for numerical
;;                                             reason, similar to what has been
;;                                             done for "H" in "modified TIP3P"

;; H1 sigma and epsilon adopted:
("butanoic acid"
 (("C1" (1.422 -0.017 0.0) 3.75 0.105 0.52) ;; 19th entry of table 5
  ("O1" (1.422 1.353 0.0) 2.96 0.21 -0.44)  ;; 20th entry of table 5
  ("O2" (2.643 -0.722 0.0) 3.0 0.17 -0.53)  ;; 21th entry of table 5
  ("C2" (0.1 -0.78 0.0) 3.5 0.066 -0.12)    ;; 3rd entry of table 1
  ("C3" (-1.06 0.212 0.0) 3.5 0.066 -0.12)
  ("C4" (-2.381 -0.551 0.0) 3.5 0.066 -0.18);; 2nd entry of table 1
  ("OH" (3.21 -0.461 0.882) 3.4 0.046 0.45) ;; FIXME: expected to match
;;                                          the last entry of table 5 but
;;                                          shouldn't they be "0.4 0.046" even
;;                                          for numerical reason?
  ("H2" (0.043 -1.407 0.89) 2.5 0.03 0.06)  ;; 6th entry of table 1
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

;;;                    2+
;;; Force fields for UO  ion with the SPC/Fw, TIP3P, TIP4p and TIP5P
;;;                    2
;;; water models. See [1] for more details. For the geometry,  since
;;; it only mentions that U-O bond length is 1.76 Å, we choose the
;;; orientation along x-axis in consistent with other solutes
;;;
;;; [1] Force Field Development for Actinyl Ions via Quantum Mechanical
;;;     Calculations: An Approach to Account for Many Body Solvation
;;;     Effects, Rai, N.; Tiwari, S. P.; Maginn, E. J., The Journal of
;;;     Physical Chemistry B 2012, 116 (35), 10885-10897.
;;;     http://dx.doi.org/10.1021/jp3028275
("uranyl, SPC"
 (("U" (0.0 0.0 0.0) 3.25 0.27 2.5)
  ("O" (-1.76 0.0 0.0) 3.00 1.08 -0.25)
  ("O" (1.76 0.0 0.0) 3.00 1.08 -0.25)))

;; Non-bond FF parameters for TIP3P uranyl is the same with SPC one according
;; to [1], the difference is that they used flexible model for SPC/Fw and
;; treated TIP3P water as a rigid molecule.
("uranyl, TIP3P"
 (("U" (0.0 0.0 0.0) 3.25 0.27 2.5)
  ("O" (-1.76 0.0 0.0) 3.00 1.08 -0.25)
  ("O" (1.76 0.0 0.0) 3.00 1.08 -0.25)))

("uranyl, TIP4P"
 (("U" (0.0 0.0 0.0) 2.78 1.39 2.5)
  ("O" (-1.76 0.0 0.0) 3.07 0.85 -0.25)
  ("O" (1.76 0.0 0.0) 3.07 0.85 -0.25)))

("uranyl, TIP5P"
 (("U" (0.0 0.0 0.0) 3.26 0.26 2.5)
  ("O" (-1.76 0.0 0.0) 2.92 1.47 -0.25)
  ("O" (1.76 0.0 0.0) 2.92 1.47 -0.25)))

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
 (("CT" 96 96) ;; "Alcohol CH3OH & RCH2OH"
  ("HC" 95 95) ;;"Methanol CH3-"
  ("OH" 93 93) ;;"Alcohol -OH"
  ("HO" 94 94))) ;;"Alcohol -OH", vdw parameters are
                 ;; 0.0 in original database
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
 (("C" 206 206)  ;; "Carboxylic Acid -COOH"
  ("O" 208 208)  ;; FIXME: the descriptions in datafile
  ("OH" 207 207) ;; of these two may need to be swapped
  ("CT2" 78 78)  ;; "Alkane -CH2-"
  ("CT3" 77 77)  ;; "Alkane CH3-"
  ("HO" 209 209) ;; "Carboxylic Acid -COOH", original vdw parameters are
                 ;; 0.0
  ("HC" 82 82))) ;; "Alkane H-C"

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
 (("CT3" 77 77) ;; "Alkane CH3-"
  ("CT2" 78 78) ;; "Alkane -CH2-"
  ("HC" 82 82))) ;; "Alkane H-C"

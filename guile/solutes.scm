;;;
;;; Simplest  species,  single  neutral   site  LJ  particle.   To  be
;;; comparable  to  water  one  might   set  β  =  0.261610  and  ρ  =
;;; 1.054796. These  numbers have been  derived from LJ  parameters of
;;; oxygen  center  of the  two-site  water  model  (see below),  room
;;; temperature  and water density  (see comments  on the  TIP3P water
;;; model).
;;;
("LJ" (("LJ" (0.0 0.0 0.0) 1.0 1.0 0.0)))

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
;; Parametrisierung TIP3P :
;; sigma_H = 0.400   epsilon_H = 0.046  q_H = 0.417
;; sigma_O = 3.1506  epsilon_O = 0.1521 q_O =-0.834
;; r_OH = 0.9572
;; r_HH = 1.5139
;; theta_HOH = 104.52
;; mass H2O = 18.0154 u
;; density = 1 kg/l = 0.6022142 u/A^3 => 0.033427745 / A^3
;; temperature : T= 298.15 K (25 C) => 0.5921, => beta =1.6889
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

("methanol"
 (("C" (-0.748 -0.015 0.024) 3.5 0.066 0.145)
  ("HC1" (-1.293 -0.202 -0.901) 2.5 0.03 0.04)
  ("HC2" (-1.263 0.754 0.6) 2.5 0.03 0.04)
  ("HC3" (-0.699 -0.934 0.609) 2.5 0.03 0.04)
  ("O" (0.558 0.42 -0.278) 3.12 0.17 -0.683)
  ("OH" (0.716 1.404 0.137) 0.4 0.04 0.418)))

;; H1 sigma and epsilon adopted:
("butanoic acid"
 (("C1" (1.422 -0.017 0.0) 3.75 0.105 0.52)
  ("O1" (1.422 1.353 0.0) 2.96 0.21 -0.44)
  ("O2" (2.643 -0.722 0.0) 3.0 0.17 -0.53)
  ("C2" (0.1 -0.78 0.0) 3.5 0.066 -0.12)
  ("C3" (-1.06 0.212 0.0) 3.5 0.066 -0.12)
  ("C4" (-2.381 -0.551 0.0) 3.5 0.066 -0.18)
  ("OH" (3.21 -0.461 0.882) 3.4 0.046 0.45)
  ("H2" (0.043 -1.407 0.89) 2.5 0.03 0.06)
  ("H3" (0.043 -1.407 -0.89) 2.5 0.03 0.06)
  ("H4" (-1.002 0.838 -0.89) 2.5 0.03 0.06)
  ("H5" (-1.002 0.838 0.89) 2.5 0.03 0.06)
  ("H6" (-2.439 -1.178 0.89) 2.5 0.03 0.06)
  ("H7" (-2.439 -1.178 -0.89) 2.5 0.03 0.06)
  ("H8" (-3.21 0.157 0.0) 2.5 0.03 0.06)))

("hexane"
 (("C" (1.709 -2.812 0.0) 3.5 0.066 -0.18)
  ("C" (1.684 -1.278 0.0) 3.5 0.066 -0.12)
  ("C" (0.245 -0.753 0.0) 3.5 0.066 -0.12)
  ("C" (0.241 0.779 0.0) 3.5 0.066 -0.12)
  ("C" (-1.198 1.304 0.0) 3.5 0.066 -0.12)
  ("C" (-1.206 2.834 0.0) 3.5 0.066 -0.18)
  ("H" (2.236 -3.164 0.887) 2.5 0.03 0.06)
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
;;; This  is a fake  solute which  tabulates some  half-way meaningful
;;; force field parameters for various sites:
;;;
("bgy3d"
 (("H" (0.0 0.0 0.0) 2.735 0.03971 0.2)
  ("Cl" (0 0.0 0.0) 3.353 0.51434 -0.2)
  ("Pd" (0.0 0.0 0.0) 2.5 0.1 0.0)))    ; fake numbers

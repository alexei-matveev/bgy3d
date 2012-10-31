;;
;; Butanoic Acid, ECP, RKS.
;;
((operations
  (operations_symm #t)
  (operations_integral #t)
  (operations_scf #t)
  (operations_properties #t))
 ;; (tasks (task "Gradients"))
 (main_options
  (integrals_on_file #f)                ; This is faster
  (relativistic "false")                ; This is an ECP calculation
  (spin_restricted #t))                 ; HCl is closed shell
 (geo
  ("C1" (2.68719131551167391874 -0.03212535327967542659 0.0) (z 6))
  ("O1" (2.68719131551167391874 2.55680017572946189315 0.0) (z 8))
  ("O2" (4.99454757165777367597 -1.36438265105445047069 0.0) (z 8))
  ("C2" (0.18897266635103192114 -1.47398679753804898496 0.0) (z 6))
  ("C3" (-2.00311026332093836418 0.40062205266418767283 0.0) (z 6))
  ("C4" (-4.49943918581807004256 -1.04123939159418588553 0.0) (z 6))
  ("OH" (6.06602258986812466889 -0.87116399187825715649 1.66673891721610154453) (z 1))
  ("H2" (0.08125824653094372609 -2.65884541555901913057 1.68185673052418409822) (z 1))
  ("H3" (.08125824653094372609 -2.65884541555901913057 -1.68185673052418409822) (z 1))
  ("H4" (-1.89350611683733984991 1.58359094402164749923 -1.68185673052418409822) (z 1))
  ("H5" (-1.89350611683733984991 1.58359094402164749923 1.68185673052418409822) (z 1))
  ("H6" (-4.60904333230166855683 -2.22609826201819051092 1.68185673052418409822) (z 1))
  ("H7" (-4.60904333230166855683 -2.22609826201819051092 -1.68185673052418409822) (z 1))
  ("H8" (-6.06602258986812466889 0.29668708617112011620 0.0) (z 1)))
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (xc_control (xc "pbe"))
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (properties
  (plot_orbitals #t))
 (orbital_plot
  (n_input_lines 0)
  (density_plot #t)
  (density_tot #t)))

;;;
;;; BGY3d input  follows. The parameter solute specifies  the table to
;;; take the  site-specific force field  parameters. The names  of the
;;; sites are those specified in PG input above:
;;;
((solute "butanoic acid"))

;;
;; Butanoic Acid, ECP, RKS.
;;
((operations
  (operations_symm #t)
  (operations_integral #t)
  (operations_scf #t)
  (operations_properties #t))
 (main_options
  (integrals_on_file #f)                ; This is faster
  (relativistic "false")                ; This is an ECP calculation
  (spin_restricted #t))                 ; HCl is closed shell
 (geo
  (units angstrom)
  ("C1" (1.422 -0.017 0.0) (z 6))
  ("O1" (1.422 1.353 0.0) (z 8))
  ("O2" (2.643 -0.722 0.0) (z 8))
  ("C2" (0.1 -0.78 0.0) (z 6))
  ("C3" (-1.06 0.212 0.0) (z 6))
  ("C4" (-2.381 -0.551 0.0) (z 6))
  ("OH" (3.21 -0.461 0.882) (z 1))
  ("H2" (0.043 -1.407 0.89) (z 1))
  ("H3" (0.043 -1.407 -0.89) (z 1))
  ("H4" (-1.002 0.838 -0.89) (z 1))
  ("H5" (-1.002 0.838 0.89) (z 1))
  ("H6" (-2.439 -1.178 0.89) (z 1))
  ("H7" (-2.439 -1.178 -0.89) (z 1))
  ("H8" (-3.21 0.157 0.0) (z 1)))
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (rep 14 (gridatom (nrad 30) (nang 131)))
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

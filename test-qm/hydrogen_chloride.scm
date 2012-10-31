;;
;; HCl, ECP, RKS.
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
  ("H" (0.0 0.0 0.0))
  ("Cl" (-2.3753855294576236 0.0 0.0)))
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (xc_control (xc "pbe"))
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "Cl" "crenbl_ecp" "ahlrichs_coulomb_fitting")
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
((solute "hydrogen chloride"))

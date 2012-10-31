;;
;; CS2, ECP, RKS.
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
  (spin_restricted #t))                 ; CS2 is closed shell
 (symmetry_group (point_group "C1"))
 (unique_atom_number (n_unique_atoms 3))
 (unique_atom (name "C") (z 6) (n_equal_atoms 1))
 (0.0 0.0 0.0)
 (unique_atom (name "S1") (z 16) (n_equal_atoms 1))
 (2.94797359507609796993 0.0 0.0)
 (unique_atom (name "S2") (z 16) (n_equal_atoms 1))
 (-2.94797359507609796993 0.0 0.0)
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (xc_control (xc "pbe"))
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "S" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "S" "crenbl_ecp" "ahlrichs_coulomb_fitting")
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
((solute "carbon disulfide"))

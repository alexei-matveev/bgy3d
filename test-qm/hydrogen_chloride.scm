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
 (symmetry_group (point_group "C1"))
 (unique_atom_number (n_unique_atoms 2))
 (unique_atom (name "H") (z 1) (n_equal_atoms 1))
 (0.0 0.0 0.0)
 (unique_atom (name "Cl") (z 17) (n_equal_atoms 1))
 (-2.3753855294576236 0.0 0.0 )
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (xc_control (xc "pbe"))
 (basis "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "Cl" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (properties
  (plot_orbitals #t))
 (orbital_plot
  (n_input_lines 0)
  (density_plot #t)
  (density_tot #t)))

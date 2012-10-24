;;
;; H2O, ECP, RKS.
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
  (spin_restricted #t))                 ; H2O is closed shell
 (symmetry_group (point_group "C1"))
 (unique_atom_number (n_unique_atoms 3))
 (unique_atom (name "O") (z 8) (n_equal_atoms 1))
 (-.55350093974217249704 0.0 0.0)
 (unique_atom (name "OH") (z 1) (n_equal_atoms 1))
 (.55350093974217249704 1.43052308427731164310 0.0)
 (unique_atom (name "OH") (z 1) (n_equal_atoms 1))
 (.55350093974217249704 -1.43052308427731164310 0.0)
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (gridatom (nrad 30) (nang 131))
 (xc_control (xc "pbe"))
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (basis "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (basis "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (properties
  (plot_orbitals #t))
 (orbital_plot
  (n_input_lines 0)
  (density_plot #t)
  (density_tot #t)))

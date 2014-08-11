;;
;; butanoic acid ECP, RKS.
;;
((operations
  (operations_symm #t)
  (operations_integral #t)
  (operations_scf #t)
  (operations_dipole #t)
  (operations_properties #t))
 (main_options
  (integrals_on_file #f)                ; This is faster
  (relativistic "false")                ; This is an ECP calculation
  (spin_restricted #t))                 ; butanoic acid is closed shell
 (geo
  (units angstrom)
  ("CT2"  (3.4057  1.8922  1.6847) (z 6))
  ("C"  (4.6993  1.3329  1.1549) (z 6))
  ("CT2"  (2.1503  1.2627  1.1037) (z 6))
  ("OH"  (4.6276  0.2556  0.3391) (z 8))
  ("O"  (5.8306  1.7309  1.3762) (z 8))
  ("CT3"  (0.9097  1.9095  1.6789) (z 6))
  ("HC"  (3.4147  1.7744  2.7880) (z 1))
  ("HC"  (3.4050  2.9864  1.5011) (z 1))
  ("HC"  (2.1505  1.3594  0.0000) (z 1))
  ("HC"  (2.1380  0.1739  1.3082) (z 1))
  ("HO"  (5.5056  0.0000  0.0723) (z 1))
  ("HC"  (0.0000  1.4538  1.2674) (z 1))
  ("HC"  (0.8647  1.8031  2.7710) (z 1))
  ("HC"  (0.8690  2.9837  1.4539) (z 1)))
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (rep 14 (gridatom (nrad 30) (nang 131)))
 (xc_control (xc "pbe"))
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
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
((solute "butanoic acid, nist")
 (solvent "water")
 (L 10.0)                               ; box size, A
 (N 128)                                 ; grid dimension
 ;; (rho 0.033427745)                   ; solvent number density, A^-3
 ;; (beta 1.6889)                       ; inverse temperature, kcal^-1
 (norm-tol 1.0e-12)
 (closure KH))

;;
;; propanoic acid ECP, RKS.
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
  (spin_restricted #t))                 ; propanoic acid is closed shell
 (geo
  (units angstrom)
  ("OH"  (1.4371  0.2295  0.0004) (z 8))
  ("C"  (1.3807  1.6210  0.0000) (z 6))
  ("O"  (0.2857  2.1943  -0.0007) (z 8))
  ("CT2"  (2.6488  2.4258  0.0002) (z 6))
  ("CT3"  (3.8613  1.4766  0.0006) (z 6))
  ("HC"  (2.6776  3.0021  -0.8165) (z 1))
  ("HC"  (2.6771  3.0023  0.8168) (z 1))
  ("HO"  (0.5165  -0.1401  -0.0003) (z 1))
  ("HC"  (4.7051  2.0131  0.0008) (z 1))
  ("HC"  (3.8338  0.8996  -0.8157) (z 1))
  ("HC"  (3.8334  0.8998  0.8171) (z 1)))
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (rep 11 (gridatom (nrad 30) (nang 131)))
 (xc_control (xc "pbe"))
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
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
((solute "propanoic acid, nist")
 (solvent "water")
 (L 10.0)                               ; box size, A
 (N 128)                                 ; grid dimension
 ;; (rho 0.033427745)                   ; solvent number density, A^-3
 ;; (beta 1.6889)                       ; inverse temperature, kcal^-1
 (norm-tol 1.0e-12)
 (closure KH))

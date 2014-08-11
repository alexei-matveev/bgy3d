;;
;; pentanoic acid ECP, RKS.
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
  (spin_restricted #t))                 ; pentanoic acid is closed shell
 (geo
  (units angstrom)
  ("CT3"  (0.9690  0.9499  1.0698) (z 6))
  ("CT2"  (2.1126  1.8447  1.4909) (z 6))
  ("CT2"  (3.4564  1.1801  1.2350) (z 6))
  ("CT2"  (4.6009  2.1010  1.6329) (z 6))
  ("C"  (5.9285  1.4318  1.3841) (z 6))
  ("OH"  (6.4076  0.6912  2.4116) (z 8))
  ("O"  (6.6289  1.4680  0.3886) (z 8))
  ("HC"  (0.0000  1.4309  1.2543) (z 1))
  ("HC"  (1.0152  0.7058  0.0000) (z 1))
  ("HC"  (0.9735  0.0000  1.6208) (z 1))
  ("HC"  (2.0574  2.8087  0.9469) (z 1))
  ("HC"  (2.0170  2.1022  2.5646) (z 1))
  ("HC"  (3.5191  0.2279  1.7986) (z 1))
  ("HC"  (3.5430  0.9031  0.1652) (z 1))
  ("HC"  (4.5558  3.0450  1.0542) (z 1))
  ("HC"  (4.5113  2.3890  2.6998) (z 1))
  ("HO"  (7.2417  0.3008  2.1695) (z 1)))
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (rep 17 (gridatom (nrad 30) (nang 131)))
 (xc_control (xc "pbe"))
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
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
((solute "pentanoic acid, nist")
 (solvent "water")
 (L 10.0)                               ; box size, A
 (N 128)                                 ; grid dimension
 ;; (rho 0.033427745)                   ; solvent number density, A^-3
 ;; (beta 1.6889)                       ; inverse temperature, kcal^-1
 (norm-tol 1.0e-12)
 (closure KH))

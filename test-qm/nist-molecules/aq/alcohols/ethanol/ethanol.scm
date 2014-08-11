;;
;; methanol, ECP, RKS.
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
  (spin_restricted #t))                 ; ethanol is closed shell
 (geo
  (units angstrom)
  ("CT3"  (1.0303  0.8847  0.9763) (z 6))
  ("CT"  (1.8847  1.9889  1.5717) (z 6))
  ("OH"  (3.1883  1.4807  1.7425) (z 8))
  ("HC"  (0.0000  1.2330  0.8324) (z 1))
  ("HC"  (0.9949  0.0000  1.6255) (z 1))
  ("HC"  (1.4097  0.5558  0.0000) (z 1))
  ("HC"  (1.9050  2.8742  0.9059) (z 1))
  ("HC"  (1.4753  2.3225  2.5456) (z 1))
  ("HO"  (3.7056  2.1820  2.1139) (z 1)))
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (rep 9 (gridatom (nrad 30) (nang 131)))
 (xc_control (xc "pbe"))
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
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
((solute "ethanol, nist")
 (solvent "water")
 (L 10.0)                               ; box size, A
 (N 128)                                 ; grid dimension
 ;; (rho 0.033427745)                   ; solvent number density, A^-3
 ;; (beta 1.6889)                       ; inverse temperature, kcal^-1
 (norm-tol 1.0e-12)
 (closure KH))

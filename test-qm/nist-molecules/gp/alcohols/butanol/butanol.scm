;;
;; butanol ECP, RKS.
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
  (spin_restricted #t))                 ; butanol is closed shell
 (geo
  (units angstrom)
  ("CT2"  (2.9651  2.0464  2.4042) (z 6))
  ("CT2"  (2.3281  2.7103  1.1934) (z 6))
  ("CT"  (0.8652  2.3158  0.9860) (z 6))
  ("CT3"  (3.3109  0.5952  2.1559) (z 6))
  ("OH"  (0.0000  2.7697  1.9922) (z 8))
  ("HC"  (3.8836  2.6000  2.6826) (z 1))
  ("HC"  (2.2886  2.1307  3.2779) (z 1))
  ("HC"  (2.8987  2.4514  0.2796) (z 1))
  ("HC"  (2.4074  3.8120  1.2923) (z 1))
  ("HC"  (0.5089  2.6739  0.0000) (z 1))
  ("HC"  (0.7244  1.2206  1.0177) (z 1))
  ("HC"  (2.4210  0.0000  1.9086) (z 1))
  ("HC"  (3.7712  0.1403  3.0423) (z 1))
  ("HC"  (4.0193  0.4827  1.3240) (z 1))
  ("HO"  (0.1229  3.7075  2.0732) (z 1)))
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (rep 15 (gridatom (nrad 30) (nang 131)))
 (xc_control (xc "pbe"))
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
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
((solute "butanol, nist")
 (solvent "water")
 (L 10.0)                               ; box size, A
 (N 128)                                 ; grid dimension
 ;; (rho 0.033427745)                   ; solvent number density, A^-3
 ;; (beta 1.6889)                       ; inverse temperature, kcal^-1
 (norm-tol 1.0e-12)
 (closure KH))

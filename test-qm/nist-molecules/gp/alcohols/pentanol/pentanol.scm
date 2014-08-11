;;
;; pentanol ECP, RKS.
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
  (spin_restricted #t))                 ; pentanol is closed shell
 (geo
  (units angstrom)
  ("CT2"  (2.0264  2.2581  1.9056) (z 6))
  ("CT2"  (3.2462  1.6488  1.2359) (z 6))
  ("CT2"  (4.4464  1.6425  2.1693) (z 6))
  ("CT"  (0.8116  2.2939  0.9786) (z 6))
  ("CT3"  (5.6613  1.0380  1.5026) (z 6))
  ("OH"  (0.3196  1.0290  0.6256) (z 8))
  ("HC"  (2.2542  3.2923  2.2314) (z 1))
  ("HC"  (1.7856  1.7003  2.8338) (z 1))
  ("HC"  (3.0191  0.6161  0.9030) (z 1))
  ("HC"  (3.4892  2.2109  0.3117) (z 1))
  ("HC"  (4.6718  2.6757  2.5011) (z 1))
  ("HC"  (4.2030  1.0807  3.0931) (z 1))
  ("HC"  (1.0538  2.7449  0.0000) (z 1))
  ("HC"  (0.0000  2.8939  1.4351) (z 1))
  ("HC"  (5.4791  0.0000  1.1934) (z 1))
  ("HC"  (5.9510  1.5990  0.6041) (z 1))
  ("HC"  (6.5258  1.0323  2.1788) (z 1))
  ("HO"  (0.1634  0.5491  1.4297) (z 1)))
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (rep 18 (gridatom (nrad 30) (nang 131)))
 (xc_control (xc "pbe"))
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
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
((solute "pentanol, nist")
 (solvent "water")
 (L 10.0)                               ; box size, A
 (N 128)                                 ; grid dimension
 ;; (rho 0.033427745)                   ; solvent number density, A^-3
 ;; (beta 1.6889)                       ; inverse temperature, kcal^-1
 (norm-tol 1.0e-12)
 (closure KH))

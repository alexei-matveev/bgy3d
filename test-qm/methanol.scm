;;
;; methanol, ECP, RKS.
;;
((operations
  (operations_symm #t)
  (operations_integral #t)
  (operations_scf #t)
  (operations_properties #t))
 (main_options
  (integrals_on_file #f)                ; This is faster
  (relativistic "false")                ; This is an ECP calculation
  (spin_restricted #t))                 ; methanol is closed shell
 (geo
  ("C" (-1.41351570457521774377 -0.02834590316661532908 0.04535344506658452653))
  ("HC1" (-2.44341685296224136725 -0.38172482931041976503 -1.70264391687469410046) (z 1))
  ("HC2" (-2.38672504662901070908 1.42485406584186387541 1.13383612666461316345) (z 1))
  ("HC3" (-1.32091908756427433543 -1.76500490384124782445 1.15084366856458236091) (z 1))
  ("O" (1.05446759779809024201 0.79368528866522921442 -0.52534407202127076573))
  ("OH" (1.35304444448643837506 2.65317653639519480249 0.25889258225508667232) (z 1)))
 (mixing (chmix 0.5) (start_after_cycle 5))
 (grid (sym_reduce #t) (weight_grads #t))
 (rep 6 (gridatom (nrad 30) (nang 131)))
 (xc_control (xc "pbe"))
 (ecp "nwchem" "C" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (bas "nwchem" "H" "crenbl_ecp" "ahlrichs_coulomb_fitting")
 (ecp "nwchem" "O" "crenbl_ecp" "ahlrichs_coulomb_fitting")
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
((solute "methanol"))

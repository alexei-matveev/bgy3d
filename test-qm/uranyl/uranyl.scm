;;
;; UO22+, RKS.
;;
((operations
  (operations-symm #t)
  (operations-integral #t)
  (operations-scf #t)
  (operations-geo-opt #t)
  (operations-dipole #t)
  (operations-properties #f)
  (operations-solvation-effect #f))
 (main-options
  (integrals-on-file #f)                ; This is faster
  (relativistic "true")                 ; This is an AE calculation
  (spin-restricted #t))                 ; UO22+ is closed shell
 (geo
  (units angstrom)
  ("U"  (0.0 0.0 0.0))
  ("OU"  (0.0 0.0  1.7179211))
  ("OU"  (0.0 0.0 -1.7179211)))
 ;;
 ;; GP min  at 1.71792111 A, E = -28110.790780703865 Hartree
 ;; PCM min at 1.73538978 A
 ;; GP+RISM at 1.7421507  A
 ;;
 (~ "U_24.19.16.11_10.7.7.4")
 (~ "O_9.5.1_5.4.1")
 (~ "O_9.5.1_5.4.1")
 (mixing (chmix 0.1) (start-after-cycle 5))
 (grid (sym-reduce #t) (weight-grads #t))
 (rep 3 (gridatom (nrad 50) (nang 291)))
 (xc-control (xc "bp"))
 (occupation (charge 2.0))
 (convergence-list
  (max-geo-iteration 0))
 (diis (diis-on #f))
 (solvation
  (smoothing "fixpva")
  (sol-start-cycle 10)
  (correction-type "None")))


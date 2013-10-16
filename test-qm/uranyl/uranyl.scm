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
 (point-group "D8H")
 (unique-atom-number 2)
 (unique-atom (name "U") (z 92) (n-equal-atoms 1))
 (0.0 0.0 0.0)
 (unique-atom (name "O") (z 8) (n-equal-atoms 2))
 (0.0 0.0 3.24640132737424)
 ;;
 ;; GP min at 1.71792111 A = 3.24640132737424 au
 ;; PCM min at 1.73538978 A = 3.27941233884931 au
 ;;
 (~ "U_24.19.16.11_10.7.7.4")
 (~ "O_9.5.1_5.4.1")
 (mixing (chmix 0.1) (start-after-cycle 5))
 (grid (sym-reduce #t) (weight-grads #t))
 (rep 2 (gridatom (nrad 50) (nang 291)))
 (xc-control (xc "pbe"))
 (occupation (charge 2.0))
 (convergence-list
  (max-geo-iteration 0))
 (diis (diis-on #f))
 (solvation
  (smoothing "fixpva")
  (sol-start-cycle 10)
  (correction-type "None")))


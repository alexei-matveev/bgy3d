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
 (point-group "D5H")                    ; auto count unique atoms
 (unique-atom (name "U") (z 92) (n-equal-atoms 1))
 (0.0 0.0 0.0)
 (unique-atom (name "OU") (z 8) (n-equal-atoms 2))
 (0.0 0.0 3.353788537052)
 (unique-atom (name "OW") (z 8) (n-equal-atoms 5))
 (4.699628346156 0.0 0.0)
 (unique-atom (name "HW") (z 1) (n-equal-atoms 10))
 (5.790662544029 0.0 1.542955359644)
 (~ "U_24.19.16.11_10.7.7.4")
 (~ "O_9.5.1_5.4.1")
 (~ "O_9.5.1_5.4.1")
 (~ "H_6.1_4.1")
 (mixing (chmix 0.1) (start-after-cycle 5))
 (grid (sym-reduce #t) (weight-grads #t))
 (gridatom (nrad 50) (nang 291))       ; the same for all unique atoms
 (xc-control (xc "bp"))
 (occupation (charge 2.0))
 (convergence-list
  (max-geo-iteration 20))
 (diis (diis-on #f))
 (solvation
  (smoothing "fixpva")
  (sol-start-cycle 10)
  (correction-type "None")))


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
 (point-group "D3D")                    ; auto count unique atoms
 (unique-atom (name "U") (z 92) (n-equal-atoms 1))
 (0.0 0.0 0.0)
 (unique-atom (name "OU") (z 8) (n-equal-atoms 2))
 (0.0 0.0 3.372889924808)
 (unique-atom (name "OW") (z 8) (n-equal-atoms 6))
 (0.0 4.563295023703 1.432395610352)
 (unique-atom (name "HW") (z 1) (n-equal-atoms 6))
 (0.0 4.988110260442 3.273753499904)
 (unique-atom (name "HW") (z 1) (n-equal-atoms 6))
 (0.0 6.157738811850 0.418090000971)
 (~ "U_24.19.16.11_10.7.7.4")
 (~ "O_9.5.1_5.4.1")
 (~ "O_9.5.1_5.4.1")
 (~ "H_6.1_4.1")
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


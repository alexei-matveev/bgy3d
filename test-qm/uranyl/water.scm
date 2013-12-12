;;
;; Water, RKS.
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
  (spin-restricted #t))                 ; water is closed shell
 (geo
  (units angstrom)
  ("OW"  (-1.6172  2.2191  -0.7126))
  ("HW"  (-2.5526  2.2260  -0.9995))
  ("HW"  (-1.2731  3.1280  -0.8244)))
 (~ "O_9.5.1_5.4.1")
 (~ "H_6.1_4.1")
 (~ "H_6.1_4.1")
 (mixing (chmix 0.1) (start-after-cycle 5))
 (grid (sym-reduce #t) (weight-grads #t))
 (rep 3 (gridatom (nrad 50) (nang 291)))
 (xc-control (xc "bp"))
 (occupation (charge 0.0))
 (convergence-list
  (max-geo-iteration 0))
 (diis (diis-on #f))
 (solvation
  (smoothing "fixpva")
  (sol-start-cycle 10)
  (correction-type "None")))


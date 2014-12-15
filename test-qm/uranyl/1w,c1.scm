;;
;; Water in SPC/E geometry, RKS.
;;
((operations
  (operations-symm #t)
  (operations-integral #t)
  (operations-scf #t)
  (operations-geo-opt #f)
  (operations-dipole #t)
  (operations-properties #f)
  (operations-solvation-effect #f))
 (main-options
  (integrals-on-file #f)                ; This is faster
  (relativistic "ADKH")                 ; This is an AE calculation
  (spin-restricted #t))                 ; Closed shell
 (geo                                   ; This is D3D geometry
  (units angstrom)
  ("OW" ( 0.000000000000  2.382682819916  0.710764543024))
  ("HW" ( 0.000000000000  2.711799660033  1.655053629471))
  ("HW" ( 0.000000000000  3.163261495114  0.085707181610)))
 ;;
 ;; This feature is new in GPL version (Dec 2014):
 ;;
 (~ ("U_24.19.16.11_10.7.7.4" "U")
    ("O_9.5.1_5.4.1" "O" "OW")
    ("H_6.1_4.1" "HW"))
 (mixing (chmix 0.1) (start-after-cycle 5))
 (grid (sym-reduce #t) (weight-grads #t))
 (gridatom (nrad 50) (nang 291))       ; the same for all unique atoms
 (xc-control (xc "bp"))
 (occupation (charge 0.0))
 (convergence-list
  (max-geo-iteration 0))
 (diis
  (mmax 8)
  (loop-start 15)
  (start-steps 5)
  (cfix 0.095)
  (threshold 15.0)
  (diis-on #t))
 (solvation
  (smoothing "fixpva")
  (sol-start-cycle 10)
  (correction-type "None")))


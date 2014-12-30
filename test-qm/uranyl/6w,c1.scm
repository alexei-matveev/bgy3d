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
 ;;
 ;; ADKH versus DKH energies for this geometry are
 ;;
 ;;     -28570.264843374138 (ADKH)
 ;;     -28570.265512536022 (DKH)
 ;;          0.42 kcal of 0.7 mH diff
 ;;
 (main-options
  (integrals-on-file #f)                ; This is faster
  (relativistic "ADKH")                 ; This is an AE calculation
  (spin-restricted #t))                 ; UO22+ is closed shell
 (geo                                   ; This is D3D geometry
  (units angstrom)
  ("U"  ( 0.000000000000  0.000000000000  0.000000000000))
  ("O"  ( 0.000000000000  0.000000000000  1.799325495835))
  ("O"  ( 0.000000000000 -0.000000000000 -1.799325495835))

  ("OW" ( 0.000000000000  2.382682819916  0.710764543024))
  ("HW" ( 0.000000000000  2.711799660033  1.655053629471))
  ("HW" ( 0.000000000000  3.163261495114  0.085707181610))

  ("OW" (-2.063463851208 -1.191341409958  0.710764543024))
  ("HW" (-2.348487395563 -1.355899830017  1.655053629471))
  ("HW" (-2.739464813582 -1.581630747558  0.085707181610))

  ("OW" ( 2.063463851208 -1.191341409958  0.710764543024))
  ("HW" ( 2.348487395563 -1.355899830017  1.655053629471))
  ("HW" ( 2.739464813582 -1.581630747557  0.085707181610))

  ("OW" ( 0.000000000000 -2.382682819916 -0.710764543024))
  ("HW" ( 0.000000000000 -2.711799660033 -1.655053629471))
  ("HW" ( 0.000000000000 -3.163261495114 -0.085707181610))

  ("OW" (-2.063463851208  1.191341409958 -0.710764543024))
  ("HW" (-2.348487395563  1.355899830017 -1.655053629471))
  ("HW" (-2.739464813582  1.581630747558 -0.085707181610))

  ("OW" ( 2.063463851208  1.191341409958 -0.710764543024))
  ("HW" ( 2.348487395563  1.355899830017 -1.655053629471))
  ("HW" ( 2.739464813582  1.581630747558 -0.085707181610)))
 ;;
 ;; This feature is new in GPL version (Dec 2014):
 ;;
 (~ ("U_24.19.16.11_10.7.7.4" "U")
    ("O_9.5.1_5.4.1" "O" "OW")
    ("H_6.1_4.1" "HW"))
 (mixing (chmix 0.1) (start-after-cycle 5))
 (grid (sym-reduce #t) (weight-grads #t))
 (gridatom (nrad 30) (nang 131))       ; the same for all unique atoms
 (xc-control (xc "bp"))
 (occupation (charge 2.0))
 (convergence-list
  (max-geo-iteration 0))
 (diis
  (mmax 8)
  (loop-start 15)
  (start-steps 5)
  (cfix 0.095)
  (threshold 15.0)
  (diis-on #f))
 (solvation
  (smoothing "fixpva")
  (sol-start-cycle 10)
  (correction-type "None")))


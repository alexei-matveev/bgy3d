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
  ("U"  (-0.4193  0.2123  0.1018))
  ("OU"  (1.1822  0.9792  0.1743))
  ("OU"  (-2.0068  -0.5497  -0.1304))
  ("OW"  (-1.6172  2.2191  -0.7126))
  ("HW"  (-2.5526  2.2260  -0.9995))
  ("HW"  (-1.2731  3.1280  -0.8244))
  ("OW"  (0.9046  -1.8307  -0.3086))
  ("HW"  (0.5819  -2.7447  -0.4413))
  ("HW"  (1.8777  -1.8501  -0.4076))
  ("OW"  (-0.1489  0.1448  -2.3461))
  ("HW"  (-0.8081  -0.1703  -2.9968))
  ("HW"  (0.6309  0.4518  -2.8510))
  ("OW"  (-0.6837  1.7241  2.0738))
  ("HW"  (-1.4562  1.8695  2.6561))
  ("HW"  (0.0204  2.3354  2.3698))
  ("OW"  (-0.6191  -1.1831  2.1599))
  ("HW"  (-1.3870  -1.7658  2.3269))
  ("HW"  (-0.0073  -1.2855  2.9163)))
 (~ "U_24.19.16.11_10.7.7.4")
 (~ "O_9.5.1_5.4.1")
 (~ "O_9.5.1_5.4.1")
 (rep 5
      (~ "O_9.5.1_5.4.1")
      (~ "H_6.1_4.1")
      (~ "H_6.1_4.1"))
 (mixing (chmix 0.1) (start-after-cycle 5))
 (grid (sym-reduce #t) (weight-grads #t))
 (rep 18 (gridatom (nrad 50) (nang 291)))
 (xc-control (xc "bp"))
 (occupation (charge 2.0))
 (convergence-list
  (max-geo-iteration 0))
 (diis (diis-on #f))
 (solvation
  (smoothing "fixpva")
  (sol-start-cycle 10)
  (correction-type "None")))


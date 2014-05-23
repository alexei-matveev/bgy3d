;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;;  ../../bgy3d -L ../../ -s ./fit.scm
;;;
(use-modules
 (guile bgy3d)              ; rism-solute
 (guile molecule)           ; find-molecule
 (guile utils)              ; ddd
 (srfi srfi-1)
 (ice-9 match)
 (ice-9 pretty-print))

(define *solute* (find-molecule "UO2_5H2O, KL2-PRSPC"))
(define *species* (molecule-species *solute* 1.0))

(pretty-print (molecule-self-energy *solute* *species*))

(define *settings*
  '((L . 20.0)
    (N . 512)
    (rho . 0.0333295)
    (beta . 1.6889)
    (dielectric . 78.4)
    (norm-tol . 1e-14)))

(define BETA (assoc-ref *settings* 'beta))

;;;
;;; These were obtained in QM calculations: CAR value is U-OW distance
;;; of  a moving  water,  CDR  is the  self-energy  of the  penta-aqua
;;; complex.  The subtracted  value is  a sum  of uranyl  in  the SAME
;;; geometry as in  the aqua complex and 5 water  energies in the SAME
;;; geometry as in aqua complex.
;;;
(define *pes*
  (map (lambda (xy)
         (let ((x (car xy))
               (y (cdr xy)))
           (cons x (- y -17879815.1717215155)))) ;  -17879819.0704410
       '((1.995446975115876320e+00 . -1.788005483369186148e+07)
         (2.125540169872842444e+00 . -1.788007364587581903e+07)
         (2.264114792366613838e+00 . -1.788008380361054093e+07)
         (2.411723789402665297e+00 . -1.788008795517193526e+07)
         (2.568956157161547349e+00 . -1.788008829860143363e+07)
         (2.736439291438425681e+00 . -1.788008642712172493e+07)
         (2.914841491106517246e+00 . -1.788008347216470167e+07)
         (3.104874624793863003e+00 . -1.788008013535733894e+07)
         (3.307296971414129594e+00 . -1.788007684228642657e+07)
         (3.522916245885863606e+00 . -1.788007384180522710e+07)
         (3.752592822113550053e+00 . -1.788007126449393481e+07)
         (3.997243166090973521e+00 . -1.788006918173339590e+07))))

(pretty-print *pes*)

;;; Names are meaningless here, only coordinates:
(define *m0* (with-input-from-file "5w,scale1.xyz" read-xyz))
(define *m1* (with-input-from-file "5w,scale2.xyz" read-xyz))

(define (rc m)
  (let ((sites (molecule-sites m)))
    (site-distance (list-ref sites 0)    ; U
                   (list-ref sites 3)))) ; OW

(pretty-print (rc *m0*))
(pretty-print (rc *m1*))

;;;
;;; FIXME: linearity assumed:
;;;
(define (make-mol d)
  (let ((d0 (rc *m0*))
        (d1 (rc *m1*)))
    (let ((a (/ (- d d0)
                (- d1 d0))))
      ;; Solute with an interpolated geometry:
      (move-molecule *solute*
                     (molecule-positions (interpolate a *m0* *m1*))))))

(pretty-print (rc (make-mol 2.0)))
(pretty-print (rc (make-mol 3.0)))

(define (pes d)
  (molecule-self-energy (make-mol d) *species*))

(pretty-print (pes 2.46))

(define (make-ff p)
  ;; (lambda  (r)
  ;;   (let* ((x (/ 1 r))
  ;;          (p6 (expt x 6))
  ;;          (p12 (* p6 p6))
  ;;          (e (* 4 p6 (- p6 1)))
  ;;          (e' (* -4 x (- (* 12 p12) (* 6 p6)))))
  ;;     (list e e')))
  (lambda (r)
    (match p
     ((c6 c3)
      (let* ((x (/ 1 r))
             (p2 (* x x))
             (p3 (* p2 x))
             (p4 (* p2 p2))
             (p6 (* p4 p2))
             (p9 (* p6 p3))
             (p12 (* p6 p6))
             (e (+ (* c6 p6)
                   (* c3 p3)))
             (e' (* -1 x (+ (* 6 c6 p6)
                            (* 3 c3 p3)))))
        (list e e')))))
  )

;; (let* ((f (make-ff '(47.8502957719869 -72.529255806734)))
;;        (f0 (lambda (x) (car (f x))))
;;        (f1 (lambda (x) (cadr (f x)))))
;;   (pretty-print (f 1.3))
;;   (pretty-print (ddd f0 1.3)))
;; (exit 0)


(define (fvec p)
  (with-fluids ((*ff* (make-ff p)))
   (pretty-print p)
   (map (lambda (xy)
          (let ((x (car xy))
                (y (cdr xy)))
            ;;
            ;;
            ;; If  you use weighted  least squares with  these weigths
            ;; the  fit can  only be accurate  at the location  of the
            ;; minimum:
            ;;
            ;;  2   2   -βE(r)
            ;; w = r * e
            ;;
            (let ((w (* x (exp (* -1/2 BETA (- y -269.228160399944)))))
                  (y1 (pes x)))
              (* 1.0 (- y1 y)))))
        *pes*)))


;;;
;;; Returns a pair of d(U-O) and corresponding self-energy:
;;;
(define (property mol)
  (let* ((sites (molecule-sites mol))
         (distance (rc mol)))
    (cons distance
          (molecule-self-energy mol *species*))))

;;; 6-2: 39.154432273093 -68.5820693441408 -- ok 5
;;; 6-3: 47.8502957719869 -72.529255806734 -- ok 3
;;; 9: 5.13454712162678 -- not ok
;;; 9-3: 10.0618648533996 -19.6235591106314 -- ok 2
(let ((p* (least-squares fvec '(47.8502957719869 -72.529255806734)
                         )))
  (pretty-print (fvec p*))
  (pretty-print p*)
  (with-fluids ((*ff* (make-ff p*)))
   (let ((scan (scan-property property
                              (make-mol 1.0)
                              (make-mol 8.0))))
     (with-output-to-file "tst4.txt"
       (lambda ()
         (for-each (lambda (xy)
                     (display (car xy))
                     (display " ")
                     (display (cdr xy))
                     (newline))
                   scan))))))

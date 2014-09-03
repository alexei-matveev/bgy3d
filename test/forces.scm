;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;;  ../bgy3d -L ../ -s ./forces.scm
;;;
(use-modules
 (guile bgy3d)                       ; hnc3d-run-solute, %null-pointer
 (guile molecule)                    ; find-molecule
 (guile utils)                       ; ldot, ddd
 (srfi srfi-1)
 (ice-9 match)
 (ice-9 pretty-print)
 (ice-9 format))

(define *settings*
  '((L . 10.0)
    (N . 96)
    (rho . 0.0333295)
    (beta . 1.6889)
    (norm-tol . 1.0e-14)
    (derivatives . #t)
    (response . #t)
    (closure . HNC)
    (verbosity . 0)))

(define *solvent* (find-molecule "OW"))

;;; With d/2 =  2 it should give mu =  14.0074968008389.  OW/OW by the
;;; same code  gives 6.86940765594207.  United  atom UA/OW (OW  with 4
;;; times the epsilon) gives 5.82278007673305,
(define (make-solute d/2)
  `("OW2" (("OW" (0.0 0.0 ,(+ d/2)) 3.16 0.1549 0.0)
           ("OW" (0.0 0.0 ,(- d/2)) 3.16 0.1549 0.0))))

;;;
;;; Call sequence for 3D RISM:
;;;
(define (do-3d solute)
  (let ((alist (hnc3d-run-solute solute
                                 *solvent*
                                 *settings*
                                 %null-pointer)))
    (map vec-destroy (assoc-ref alist 'GUV))
    (bgy3d-pot-destroy (assoc-ref alist 'POTENTIAL))
    (bgy3d-restart-destroy (assoc-ref alist 'RESTART))
    ;; Return free energy and gradients computed in two different
    ;; ways:
    (let ((e (assoc-ref alist 'free-energy))
          (g (assoc-ref alist 'free-energy-gradient))
          (r (assoc-ref alist 'free-energy-response)))
      (list e g r))))

(define (pes-3d x)
  (do-3d (make-solute x)))

;; (pretty-print (pes-3d 2.0))
;; (exit 1)

;;; Excess  chemical potential for  the solvent  in solvent  and united
;;; atom limit:
(if #f
    (let ((ow (do-3d *solvent*))
          (ua (do-3d (find-molecule "OW2-UA")))
          (e0 (pes-3d 0.0))
          (e8 (pes-3d (* 0.5 (env-ref *settings* 'L))))) ; "infinite" separation
      (pretty-print/serial (list 'OW: ow))
      (pretty-print/serial (list 'UA: ua))
      (pretty-print/serial (list 'E0: e0))
      (pretty-print/serial (list 'E8: e8))
      (exit 0)))

(if #t
    (let* ((grid (mesh/lin 0.0 4.0 40))
           (results (map (lambda (x) (list x (pes-3d x)))
                         grid)))
      ;; Only one process should do IO:
      (begin/serial
       (pretty-print results)
       (with-output-to-file "forces.out"
         (lambda ()
           (for-each
            (match-lambda
             ((d/2 (e ((_ _ ga) (_ _ gb)) ((_ _ ra) (_ _ rb)))) ; ->
              (format #t "~A ~A ~A ~A ~A ~A\n" d/2 e ga gb ra rb)))
            results))))
      (exit 0)))

;;;
;;; SPC/E water geometry  as recorded is recovered by  (zmatrix 1 (* 2
;;; (atan (sqrt 2)))) or (zmatrix 1.0 1.91063323624902).
;;;
(define (zmatrix oh hoh)
  (let ((x (* oh (cos (* 1/2 hoh)) 1/2))
        (y (* oh (sin (* 1/2 hoh)))))
    (list (list (- x) 0.0 0.0)
          (list (+ x) (+ y) 0.0)
          (list (+ x) (- y) 0.0))))

;;; Stretching:
(define (stretching x)
  (zmatrix (+ 1 x) (* 2 (atan (sqrt 2)))))

;;; Bending:
(define (bending x)
  (zmatrix 1 (+ x (* 2 (atan (sqrt 2))))))

;;; Translation:
(define (translate vs dv)
  (map (lambda (v) (map + v dv)) vs))

(define (make-translation v)
  (let ((ref (zmatrix 1 (* 2 (atan (sqrt 2))))))
    (lambda (s)
      (translate ref (map (lambda (x) (* s x)) v)))))

;;; Translations:
(define trans-x (make-translation '(1 0 0)))
(define trans-y (make-translation '(0 1 0)))
(define trans-z (make-translation '(0 0 1)))

(define mode trans-z)

;;;
;;; Re-define solvent, and solute:
;;;
(set! *settings* (env-set 'dielectric 78.4 *settings*))
(set! *solvent* (find-molecule "water, cSPC/E"))
(set! make-solute
      (let ((solute (find-molecule "water, cSPC/E"))) ; let over lambda
        (lambda (x)
          (move-molecule solute (mode x)))))

;;;
;;; dot (g, df/dx) = d/dx dot (g, f)
;;;
(define (mode-gradient g x)
  (and g                                ; return #f is g is #f
       (let ((f (lambda (x)
                  (ldot (mode x) g))))
         (ddd f x))))                   ; not at x = 0.0!

(let* ((grid (linspace -0.125 0.125 15))
       (results (map (lambda (x) (list x (pes-3d x)))
                     grid)))
  ;; Only one process should do IO:
  (begin/serial
   (pretty-print results)
   (with-output-to-file "forces-water.out"
     (lambda ()
       (for-each
        (match-lambda
         ((x (e grad resp)) ; ->
          (format #t "~A ~A ~A\n" x e (mode-gradient grad x))))
        results)))))


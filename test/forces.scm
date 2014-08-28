;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;;  ../bgy3d -L ../ -s ./forces.scm
;;;
(use-modules
 (guile bgy3d)                       ; hnc3d-run-solute, %null-pointer
 (guile molecule)                    ; find-molecule
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
    ;; Return free energy:
    (let ((e (assoc-ref alist 'free-energy))
          (g (assoc-ref alist 'free-energy-gradient))
          (r (assoc-ref alist 'free-energy-response)))
      (list e g r))))

(define (pes-3d d/2)
  (do-3d (make-solute d/2)))

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
      (pretty-print/serial (list 'E8: e8))))
;; (exit 1)

(define *grid* (mesh/lin 0.0 4.0 40))

(define *results* (map (lambda (x)
                         (list x (pes-3d x)))
                       *grid*))

;;; Only one process should do IO:
(begin/serial
 (pretty-print *results*)
 (with-output-to-file "forces.out"
   (lambda ()
     (for-each
      (lambda (row)
        (match row
          ((d/2 (e ((_ _ ga) (_ _ gb)) ((_ _ ra) (_ _ rb))))
                (format #t "~A ~A ~A ~A ~A ~A\n" d/2 e ga gb ra rb))))
      *results*))))


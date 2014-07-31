;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;; FIXME: The default solver does not converge, see snes-solver in
;;; settings.
;;;
;;;  ../bgy3d -L ../ -s ./run-ions.scm
;;;
(use-modules
 (guile bgy3d)          ; rism-solute, hnc3d-run-solute, %null-pointer
 (guile molecule)       ; find-molecule
 (ice-9 pretty-print))

(define *settings*
  '((L . 160.0)
    (N . 4096)
    (rho . 0.0333295)
    (beta . 1.6889)
    (dielectric . 78.4)
    (norm-tol . 1.0e-14)
    (max-iter . 1500)
    (damp-start . 1.0)
    (lambda . 0.02)
    (bond-length-thresh . 2.0)
    (derivatives . #f)
    (snes-solver . "trial")))

(define *solvent*
  (find-molecule "water, cSPC/E"))

(define *names*
  '("Li+, JC08, SPC/E"
    "Na+, JC08, SPC/E"
    "K+, JC08, SPC/E"
    "Rb+, JC08, SPC/E"
    "Cs+, JC08, SPC/E"
    "F-, JC08, SPC/E"
    "Cl-, JC08, SPC/E"
    "Br-, JC08, SPC/E"
    "I-, JC08, SPC/E"))

(define *ions*
  (map find-molecule *names*))

;;; Solvent susceptibility  for both closures.  Two  pure solvent runs
;;; here.  Previous version  used  the same  (HNC) susceptibility  for
;;; doing KH and  HNC ions in water calculations.  The KH numbers were
;;; thus slightly wrong.
(define *chi*
  (map (lambda (closure)
         (let ((alist (rism-solvent *solvent*
                                    (env-set 'closure closure *settings*))))
           (cons closure
                 (assoc-ref alist 'susceptibility))))
       (list 'HNC 'KH)))

;;;
;;; Call sequence for 1D RISM:
;;;
(define (run-1d solute closure)
  (let* ((settings (env-set 'closure closure *settings*))
         (chi (assoc-ref *chi* closure))
         (results (rism-solute solute *solvent* settings chi))
         (solute-properties (assoc-ref results 'solute))
         (free-energy (assoc-ref solute-properties 'XXX)))
    free-energy))

;;;
;;; Call sequence for 3D RISM:
;;;
(define (run-3d solute closure)
  ;; 3D cannot handle large dimensions. FIXME: literals here:
  (let ((settings (env-set 'N 96 (env-set 'L 10.0 (env-set 'closure closure *settings*)))))
    (let ((alist (hnc3d-run-solute solute
                                   *solvent*
                                   settings
                                   %null-pointer)))
      (map vec-destroy (assoc-ref alist 'GUV))
      (bgy3d-pot-destroy (assoc-ref alist 'POTENTIAL))
      (bgy3d-restart-destroy (assoc-ref alist 'RESTART))
      ;; Return free energy:
      (assoc-ref alist 'XXX))))

(define *hnc-1d* (map (lambda (m) (run-1d m 'HNC)) *ions*))
(define *kh-1d* (map (lambda (m) (run-1d m 'KH)) *ions*))
(define *hnc-3d* (map (lambda (m) (run-3d m 'HNC)) *ions*))
(define *kh-3d* (map (lambda (m) (run-3d m 'KH)) *ions*))

;;; FIXME: tty is polluted with all  kinds of prints, so we also write
;;; selected data into a separate file. See ./out/run-ions-out for the
;;; reference output.
(let* ((hnc-1d (map cons *names* *hnc-1d*))
       (kh-1d (map cons *names* *kh-1d*))
       (hnc-3d (map cons *names* *hnc-3d*))
       (kh-3d (map cons *names* *kh-3d*))
       (out `((KH-1D ,kh-1d)
              (KH-3D ,kh-3d)
              (HNC-1D ,hnc-1d)
              (HNC-3D ,hnc-3d)
              )))
  (pretty-print/serial out)
  ;;
  ;; One should not be writing to the same file on multiple workers:
  (begin/serial
   (with-output-to-file "run-ions-out"
     (lambda () (pretty-print out)))))


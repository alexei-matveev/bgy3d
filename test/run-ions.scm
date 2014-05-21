;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;; FIXME: The default solver does not converge, see snes-solver in
;;; settings.
;;;
;;;  ../bgy3d -L ../ -s ./run-ions.scm
;;;
(use-modules
 (guile bgy3d)              ; rism-solute
 (guile molecule)           ; find-molecule
 (ice-9 pretty-print))

(define *settings*
  '((L . 160.0)
    (N . 4096)
    (rho . 0.033295)
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

;;; Pure solvent run here:
(define *solvent-properties*
  (rism-solvent *solvent* *settings*))

;;; Solvent susceptibility:
(define *chi*
  (assoc-ref *solvent-properties* 'susceptibility))


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

(define (run-one solute closure)
  (let* ((settings (acons 'closure closure *settings*))
         (results (rism-solute solute *solvent* settings *chi*))
         (solute-properties (assoc-ref results 'solute))
         (free-energy (assoc-ref solute-properties 'XXX)))
    free-energy))

(define *hnc* (map (lambda (m) (run-one m 'HNC)) *ions*))
(define *kh* (map (lambda (m) (run-one m 'KH)) *ions*))

;;; FIXME: tty is polluted with all  kinds of prints, so we also write
;;; selected data into a separate file. See ./out/run-ions-out for the
;;; reference output.
(let* ((hnc (map cons *names* *hnc*))
       (kh (map cons *names* *kh*))
       (out `((KH ,kh) (HNC ,hnc))))
  (pretty-print out)
  (with-output-to-file "run-ions-out"
    (lambda () (pretty-print out))))


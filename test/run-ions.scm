;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;; FIXME: The default  solver does not converge, but  there is no way
;;; to choose except  from command line by means  of PETSC environment
;;; (beware of typos in command line, there is no checking):
;;;
;;;  ../bgy3d -L ../ -s ./run-ions.scm --snes-solver trial
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
    (closure . HNC)
    (derivatives . #f)))

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

(define (run-one solute)
  ;; (pretty-print solute)
  ;; (pretty-print *solvent*)
  ;; (pretty-print *settings*)
  (let* ((results (rism-solute solute *solvent* *settings* *chi*))
         (solute-properties (assoc-ref results 'solute))
         (free-energy (assoc-ref solute-properties 'XXX)))
    free-energy))

(define *energies* (map run-one *ions*))

;;; FIXME: tty is polluted with all  kinds of prints, so we also write
;;; selected data into a separate file
(let ((out (map cons *names* *energies*)))
  (pretty-print out)
  (with-output-to-file "run-ions-out"
    (lambda () (pretty-print out))))


;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;; FIXME: The default solver does not converge, see snes-solver in
;;; settings.
;;;
;;;  ../bgy3d -L ../ -s ./run-ions.scm
;;;
(use-modules
 (srfi srfi-1)
 (guile bgy3d)          ; rism-solute, hnc3d-run-solute
 (guile molecule)       ; find-molecule
 (guile utils)          ; numbers->strings
 (guile amber)          ; from mobley09
 (ice-9 pretty-print))

(when #f
      (let* ((results (with-input-from-file "species" read))
             (entries (map car results))
             (results (map cdr results))
             (energies (map (lambda (d) (assoc-ref d 'free-energy))
                            results))
             (kb-integrals (map (lambda (d) (assoc-ref d 'excess-coordination))
                                results)))
        (for-each
         (lambda (name energy kb)
           (format #t "~S,~A,~A\n" name energy kb))
         entries
         energies
         kb-integrals))
      (exit 0))

(define *settings*
  '((L . 10.0)
    (N . 96)
    (rho . 0.0333295)
    (beta . 1.6889)
    (dielectric . 78.4)
    (closure . KH)
    (norm-tol . 1.0e-14)
    (max-iter . 1500)
    (damp-start . 1.0)
    (lambda . 0.02)
    (bond-length-thresh . 2.0)
    (derivatives . #f)
    (snes-solver . "newton")))

(define *solvent*
  (find-molecule "water, PR-SPC/E"))

;;; Solvent susceptibility:
(define *chi*
  (let ((alist (rism-solvent *solvent* *settings*)))
    (with-output-to-file "solvent"
      (lambda () (pretty-print alist)))
    (assoc-ref alist 'susceptibility)))

;;;
;;; Call sequence for 1D RISM:
;;;
(define (run-1d entry)
  (let* ((alist (rism-solute (make-solute entry)
                             *solvent*
                             *settings*
                             *chi*)))
    (cons entry alist)))

;;;
;;; Call   sequence   for   3D   RISM.    3D   cannot   handle   large
;;; dimensions. Check the settings!
;;;
(define (run-3d entry)
  (let ((alist (hnc3d-run-solute (make-solute entry)
                                 *solvent*
                                 *settings*)))
    (destroy alist)
    (cons entry (list (assoc 'free-energy alist)))))


(let ((results (map run-3d entries)))
  (with-output-to-file "species"
    (lambda () (pretty-print results))))
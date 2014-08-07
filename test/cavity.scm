;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;; The default solver may not converge, see snes-solver in settings.
;;;
;;;  ../bgy3d -L ../ -s ./cavity.scm
;;;
(use-modules
 (guile bgy3d)          ; rism-solute, hnc3d-run-solute, %null-pointer
 (guile molecule)       ; find-molecule
 (guile utils)          ; numbers->strings
 (ice-9 match)
 (ice-9 pretty-print)
 (ice-9 format))

;;;
;;; Make a  hard-sphere potential if any  of the sites  is named "HS",
;;; otherwise return #f telling to use the default recepie:
;;;
(define (make-ff a b)
  (define (hs a b)
    (and (equal? (site-name a) "HS")
         (let ((d (site-sigma a)))
           (lambda (rab)
             (if (> rab d)
                 0.0
                 +inf.0)))))
  (or (hs a b)
      (hs b a)))

(define *settings*
  `((L . 320.0)
    (N . 16384)
    (rho . 0.0333295)
    (beta . 1.6889)
    ;; (dielectric . 78.4)
    (norm-tol . 1.0e-14)
    (max-iter . 1500)
    (damp-start . 1.0)
    (lambda . 0.02)
    (bond-length-thresh . 2.0)
    (derivatives . #f)
    (closure . HNC)
    (snes-solver . "trial")
    (force-field-short . ,make-ff)
    (verbosity . 0)))

(define *solvent*
  (find-molecule "OW"))

;;; Precompute solvent susceptibility. Pure solvent run here:
(define *properties* (rism-solvent *solvent* *settings*))
(define chi (assoc-ref *properties* 'susceptibility))
(define kappa (assoc-ref *properties* 'compressibility))
(define kappa-error (assoc-ref *properties* 'compressibility-error))
(define beta (assoc-ref *settings* 'beta))
(define rho (assoc-ref *settings* 'rho))

;;;
;;;                    -1        3              -3
;;; 1 - β/ρκ, β in kcal  , κ in A / kcal, ρ in A
;;;
(format #t "beta = ~A kcal^-1\n" beta)
(format #t "rho = ~A A^-3\n" rho)
(format #t "kappa = ~A A^3/kcal, error = ~A A^3/kcal (1)\n" kappa kappa-error)

(define excess-coordination (assoc-ref *properties* 'excess-coordination))
(define excess-coordination-error (assoc-ref *properties* 'excess-coordination-error))
(define kappa (/ (* beta (+ 1 excess-coordination)) rho))
(define kappa-error (/ (* beta excess-coordination-error) rho))
(format #t "kappa = ~A A^3/kcal, error = ~A A^3/kcal (2)\n" kappa kappa-error)

(format #t "rho * kappa = ~A kcal^-1\n" (* rho kappa))
(format #t "beta / (rho * kappa) = ~A\n" (/ beta (* rho kappa)))
(define nc0 (- 1 (/ beta (* rho kappa))))
(format #t "rho * c(0) = 1 - beta / (rho * kappa) = ~A\n" nc0)

;;;
;;; A = -2.99 kcal (HNC) or -4.19 kcal (KH):
;;;
(define A (/ nc0 beta 2))
(format #t "A = kT * rho * c(0) / 2 = ~A kcal\n" A)

;;; With sigma = 3.16 it should give mu = 6.44 kcal (KH):
(define (make-solute sigma)
  `("HS" (("HS" (0.0 0.0 0.0) ,sigma .59210136775415951210 0.0)))) ; 0.1549

;;;
;;; Call sequence for 1D RISM:
;;;
;;; (set! *settings* (env-set 'verbosity 2 *settings*))
(define (run-1d sigma)
  (let* ((solute (make-solute sigma))
         (alist (rism-solute solute *solvent* *settings* chi))
         (free-energy (assoc-ref alist 'free-energy))
         (excess-coordination (assoc-ref alist 'excess-coordination)))
    (pretty-print alist)
    (list sigma
          free-energy
          excess-coordination
          (+ free-energy (* A (- excess-coordination))))))

;; (run-1d 3.16)
;; (exit 1)


(define *grid*
  (list
        0.125 0.25 0.5 1. 2. 3.0 3.16 4. 5. 6. 7. 8. 9. 10.
        20. 40.
        60. 80.
        ))


(define *results*
  (map run-1d *grid*))

(pretty-print *results*)
(for-each
 (lambda (row)
   (match row
          ((s mu n e)
           (format #t "~A ~A ~A ~A\n" s mu n e))))
 *results*)


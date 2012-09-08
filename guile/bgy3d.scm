#!/users/alexei/darcs/bgy3d/bgy3d -s
!#

(set! %load-path (cons "/users/alexei/darcs/bgy3d/guile" %load-path))

(use-modules (srfi srfi-1)              ; list manipulation
             (ice-9 pretty-print))

;;
;; FIXME: at the moment this function only emulates the minimum of the
;; functionality of the original  executable. The new functionality is
;; in flux:
;;
(define (old-main argv)
 (cond
  ;;
  ((member "--BGY2Site" argv)
   (bgy3d-run-solvent '()))             ; Use defaults and Petsc env
  ;;
  ((member "--BGYM2Site" argv)
   (let ((h-cl                       ; Assuming first entry is for HCl
          (car (slurp (find-file "solutes.scm")))))
     (bgy3d-run-solute h-cl '()))) ; Use defaults and Petsc env
  ;;
  (else
   (new-main argv))))

;;
;; Find a file in the search patch, or die:
;;
(define (find-file file)
  (or (search-path %load-path file)
      (error "Not found:" file)))

;;
;; Slurps the whole file into a list:
;;
(define (slurp file)
  (with-input-from-file file
    (lambda ()
      (let loop ((acc (list)))
        (let ((sexp (read)))
          (if (eof-object? sexp)
              (reverse acc)
              (loop (cons sexp acc))))))))

;; (write (slurp "guile/solutes.scm"))
;; (newline)

;;
;; Solute parameters are currently represented by a name and a list of
;; sites:
;;
;; ("water"
;;  (("O" (-0.2929 0.0 0.0) 3.1506 0.1521 -0.834)
;;   ("OH" (0.2929 0.757 0.0) 0.4 0.046 0.417)
;;   ("OH" (0.2929 -0.757 0.0) 0.4 0.046 0.417)))
;;
(define (site-name site) (first site))
(define (site-position site) (second site))
(define (site-x site) (first (site-position site)))
(define (site-y site) (second (site-position site)))
(define (site-z site) (third (site-position site)))

(define (print-xyz solute)
  (let ((sites (second solute)))
    (format #t "~a\n" (length sites))
    (format #t "# ~a\n" (first solute))
    (for-each (lambda (site)
                (format #t "~a ~a ~a ~a\n"
                        (site-name site)
                        (site-x site)
                        (site-y site)
                        (site-z site)))
              sites)))

;; (for-each print-xyz (slurp (find-file "solutes.scm")))
;; (exit 0)

;;;
;;; Ignores command line argumens. Petsc environment respects them:
;;;
(define (new-main argv)
  (let ((settings       ; Settings are handled as an association list:
         '((N . 32)
           (rho . 0.018)
           (beta . 1.1989)
           (norm-tol . 1.0e-2)
           (max-iter . 200)
           (L . 10.0)
           (zpad . 10.0)
           (damp-start . 1.0)
           (lambda . 0.02)))
        (solutes                        ; all entries in the file
         (slurp (find-file "solutes.scm"))))
    ;;
    ;; At the moment both functions, bgy3d-run-solvent and
    ;; bgy3d-run-solute, echo settings as is:
    ;;
    (pretty-print solutes)
    (display "Processing pure solvent. Settings:\n")
    (pretty-print settings)
    (bgy3d-run-solvent settings)
    (for-each (lambda (solute)          ; process each solute ...
                (display "Processing solute description:\n")
                (pretty-print solute)
                (print-xyz solute)
                (bgy3d-run-solute solute settings))
              solutes)))                ; ... from this list

;;;
;;; Trying to  emulate behaviour  of old executable  unless we  find a
;;; better interface:
;;;
(old-main (command-line))

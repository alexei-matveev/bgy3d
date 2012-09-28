;;;
;;; Scheme  interface  to  BGY3d  code.   Not to  pollute  the  global
;;; namespace  we put bgy3d-*  functions into  this module.
;;;
(define-module (guile bgy3d)
  #:export (new-main
            old-main
            bgy3d-run-solvent
            bgy3d-run-solute
            bgy3d-run))

(use-modules (srfi srfi-1)              ; list manipulation
             (ice-9 pretty-print))

;;;
;;; This name has to be defined on guile startup, see the C sources of
;;; bgy3d_guile_init() in bgy3d-guile.c:
;;;
(define guile-bgy3d-module-init
  (@@ (guile-user) guile-bgy3d-module-init))

;;;
;;; The list of the procedures defined by the next call includes:
;;;
;;;   bgy3d-run-solvent
;;;   bgy3d-run-solute
;;;
;;; and posissibly more, depending on the compilation options.
;;;
(guile-bgy3d-module-init)

;;;
;;; Settings  are  handled  as  an  association list,  these  are  the
;;; settings used in regression tests:
;;;
(define bgy3d-settings
  '((N . 32)                            ; grid dimension
    (rho . 0.018)                       ; solvent density
    (beta . 1.1989)                     ; inverse temperature
    (norm-tol . 1.0e-2)                 ; convergence threshold
    (max-iter . 200)                    ; max number of iterations
    (L . 10.0)                          ; [-L, L] gives the box size
    (zpad . 10.0)                       ; ??
    (damp-start . 1.0)                  ; scaling factor?
    (lambda . 0.02)))                   ; not the scheme lambda
;;
;; FIXME: at the moment this function only emulates the minimum of the
;; functionality of the original  executable. The new functionality is
;; in flux.
;;
;; Also note  that find-file assumes the %load-path  contains the repo
;; directory in order to find "guile/solutes.scm".
;;
(define (old-main argv)
 (cond
  ;;
  ((member "--BGY2Site" argv)
   (bgy3d-run-solvent '()))             ; Use defaults and Petsc env
  ;;
  ((member "--BGYM2Site" argv)
   (let ((h-cl   ; find entry with "hydrogen cloride" in car position:
          (assoc "hydrogen chloride"
                 (slurp (find-file "guile/solutes.scm"))))
         (g1-files
          (list "g0.bin" "g1.bin")))
     (map bgy3d-vec-save
          g1-files
          (bgy3d-run-solute h-cl '())))) ; Use defaults and Petsc env
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
(define (make-solute name sites)
  (list name sites))

(define (solute-name solute) (first solute))
(define (solute-sites solute) (second solute))

(define (make-site name position sigma epsilon charge)
  (list name position sigma epsilon charge))

(define (site-name site) (first site))
(define (site-position site) (second site))
(define (site-sigma site) (third site))
(define (site-epsilon site) (fourth site))
(define (site-charge site) (fifth site))

(define (site-x site) (first (site-position site)))
(define (site-y site) (second (site-position site)))
(define (site-z site) (third (site-position site)))

(define (print-xyz solute)
  (let ((sites (solute-sites solute)))
    (format #t "~a\n" (length sites))
    (format #t "# ~a\n" (solute-name solute))
    (for-each (lambda (site)
                (format #t "~a ~a ~a ~a\n"
                        (site-name site)
                        (site-x site)
                        (site-y site)
                        (site-z site)))
              sites)))

;; (for-each print-xyz (slurp (find-file "solutes.scm")))
;; (exit 0)

(define (update-sites solute)
  (let* ((solutes
          (slurp (find-file "guile/solutes.scm")))
         (table                     ; fake solute with site-parameters
          (second (assoc "bgy3d" solutes)))
         (update-one (lambda (site)
                       ;; (pretty-print site)
                       (let* ((name (site-name site))
                              (table-site (assoc name table)))
                         (make-site name ; original, same as in table
                                    (site-position site)    ; original
                                    (site-sigma table-site) ; from table
                                    (site-epsilon table-site) ; from table
                                    (site-charge table-site)))))) ; from table
    ;; (pretty-print solutes)
    ;; (pretty-print table)
    (make-solute (solute-name solute)
                 (map update-one
                      (solute-sites solute)))))

;;
;; Will   not  work   with  distributed   vectors.   The   problem  is
;; bgy3d-vec-ref (VecGetValues) does not either.
;;
(define (vec-fold kons knil vec)
  "A (left) fold of a (Petsc) vector.  E.g. (vec-fold + 0.0 vec)
computes the sum of all vector elements."
  (let ((len (bgy3d-vec-length vec))
        (vec (lambda (i) (bgy3d-vec-ref vec i))))
    (let loop ((knil knil)
               (i 0))
      (if (< i len)
          (loop (kons (vec i) knil)
                (+ 1 i))
          knil))))

;;
;; BGY3d  code operates  in angstroms,  QM codes  use atomic  units by
;; convention:
;;
(define (bohr->angstrom x) (* x 0.52917706))
(define (angstrom->bohr x) (/ x 0.52917706))

;;
;; Brain-dead implementation of cubic root:
;;
(define (cubic-root n)
  (let loop ((x 0))
    (if (< n (* x x x))
        (- x 1)
        (loop (+ 1 x)))))
;;
;; This writes a GAMESS-UK punch file to the current output port:
;;
(define (write-punch-file solute vecs settings)
  (define (header alist)
    (for-each (lambda (pair)
                (format #t "~a=~a " (car pair) (cdr pair)))
              alist)
    (format #t "\n"))
  ;;
  ;; File title:
  ;;
  (header '((block . title) (records . 1)))
  (format #t "This file was generated by BGY3d\n")
  ;;
  ;; Description of the solute:
  ;;
  (header '((block . fragment) (records . 0)))
  (header '((block . title) (records . 1)))
  (format #t "~a\n" (solute-name solute))
  (header `((block . coordinates)
            (records . (unquote (length (solute-sites solute))))))
  (for-each (lambda (site)
              (format #t "~a ~a ~a ~a\n"
                      (site-name site)
                      (angstrom->bohr (site-x site)) ; punch file is in AU
                      (angstrom->bohr (site-y site))
                      (angstrom->bohr (site-z site))))
            (solute-sites solute))
  ;;
  ;; Description of the bonds:
  ;;
  (header '((block . connectivity) (records . 0)))
  ;; (format #t "1 2\n")                   ; FIXME: real ones
  ;;
  ;; Grid data:
  ;;
  (let loop ((vecs vecs)
             (vec-id 0))
    (if (not (null? vecs))
        (let* ((vec (first vecs))
               (len (bgy3d-vec-length vec))
               (n (cubic-root len))
               (L (or (assq-ref settings 'L) 10.0))
               (L (angstrom->bohr L))   ; punch file is in AU
               (S (* L (/ (- n 1) n)))) ; FIXME: off-by-one?
          (header '((block . data) (records . 0)))
          (header '((block . grid_title) (records . 1)))
          (format #t "Solvent site hole density ~a\n" vec-id)
          (header '((block . grid_axes) (records . 3)))
          (format #t "~a   0.000000 ~a 0 au xaxis\n" n (* 2 S))
          (format #t "~a   0.000000 ~a 0 au yaxis\n" n (* 2 S))
          (format #t "~a   0.000000 ~a 0 au zaxis\n" n (* 2 S))
          (header '((block . grid_mapping) (records . 3)))
          (format #t "~a ~a ~a  ~a ~a ~a\n" (- L) (- L) (- L) (+ S) (- S) (- S))
          (format #t "~a ~a ~a  ~a ~a ~a\n" (- L) (- L) (- L) (- S) (+ S) (- S))
          (format #t "~a ~a ~a  ~a ~a ~a\n" (- L) (- L) (- L) (- S) (- S) (+ S))
          ;;
          ;; Actual numeric data:
          ;;
          (header `((block . grid_data)
                    (records . (unquote len))
                    (elements . 1)))
          (vec-fold (lambda (x seed) (format #t "~a\n" (- x 1.0)))
                    #f
                    vec)
          (loop (cdr vecs)
                (+ 1 vec-id))))))

;;;
;;; This hook is called from PG:
;;;
(define (bgy3d-run . args)
  "To be called from QM code."
  (let ((settings bgy3d-settings)
        (solute args))                  ; FIXME!
    (if (equal? (solute-name solute) "bgy3d") ; only when called by PG
        (set! solute (update-sites solute))) ; update force-field params
    (pretty-print solute)
    (force-output)
    ;;
    ;; At the moment  the function bgy3d-run-solvent echos settings as
    ;; is, the output is written to disk instead:
    ;;
    (bgy3d-run-solvent settings)        ; writes g??.bin files to disk
    ;;
    ;; The function bgy3d-run-solute allocates and returns two Petsc
    ;; Vecs in a list, it is the callers responsibility to destroy
    ;; them:
    ;;
    (let ((g1 (bgy3d-run-solute solute settings))) ; reads g??.bin
      ;;
      ;; Use g1 vectors to produce a *.pun file for visualization:
      ;;
      (with-output-to-file "plot.pun"
        (lambda () (write-punch-file solute g1 settings)))
      ;;
      ;; Dont forget to destroy them after use:
      ;;
      (map bgy3d-vec-destroy g1))))

;;;
;;; Ignores command line argumens. Petsc environment respects them:
;;;
(define (new-main argv)
  ;; (pretty-print argv)
  (let* ((solutes (slurp (find-file "guile/solutes.scm")))
         (run (lambda (name)
                (let ((solute (assoc name solutes)))
                  (if solute
                      (apply bgy3d-run solute)
                      (error "no such solute: " name))))))
    (map run (cdr argv))))


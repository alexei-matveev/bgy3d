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
             (ice-9 pretty-print)
             (ice-9 getopt-long))

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
;;;   bgy3d-vec-destroy
;;;   bgy3d-vec-save
;;;   bgy3d-vec-load
;;;   bgy3d-vec-length
;;;   bgy3d-vec-ref
;;;   bgy3d-rank
;;;   bgy3d-size
;;;
;;; and posissibly more, depending on the compilation options.
;;;
(guile-bgy3d-module-init)

;;;
;;; Specifications of command line flags for use with getopt-long:
;;;
(define option-spec
  '((solute (value #t))                 ; a string
    (N (value #t))                      ; grid dimension
    (rho (value #t))                    ; solvent density
    (beta (value #t))                   ; inverse temperature
    (norm-tol (value #t))               ; convergence threshold
    (max-iter (value #t))               ; max number of iterations
    (L (value #t))                      ; [-L, L] gives the box size
    (zpad (value #t))                   ; ??
    (damp-start (value #t))             ; scaling factor?
    (lambda (value #t))                 ; mixing parameter
    (BGY2Site (value #f))               ; pure solvent run
    (BGYM2Site (value #f))))            ; solute + solvent run

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

;;;
;;; Find a solute in a database or die:
;;;
(define (find-solute name)
  (or (assoc name (slurp (find-file "guile/solutes.scm")))
      (error "No such solute: " name)))

;;;
;;; Find  a file in  the search  patch, or  die.  Note  that find-file
;;; assumes  the %load-path contains  the repo  directory in  order to
;;; find "guile/solutes.scm".
;;;
(define (find-file file)
  (or (search-path %load-path file)
      (error "Not found:" file)))

;;;
;;; Slurps the whole file into a list:
;;;
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

;;;
;;; Solute parameters are  currently represented by a name  and a list
;;; of sites:
;;;
;;; ("water"
;;;  (("O" (-0.2929 0.0 0.0) 3.1506 0.1521 -0.834)
;;;   ("OH" (0.2929 0.757 0.0) 0.4 0.046 0.417)
;;;   ("OH" (0.2929 -0.757 0.0) 0.4 0.046 0.417)))
;;;
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
  (let ((name (solute-name solute))
        (sites (solute-sites solute)))
    (format #t "~a\n" (length sites))
    (format #t "# ~a\n" name)
    (for-each (lambda (site)
                (format #t "~a ~a ~a ~a\n"
                        (site-name site)
                        (site-x site)
                        (site-y site)
                        (site-z site)))
              sites)))

;; (for-each print-xyz (slurp (find-file "guile/solutes.scm")))
;; (exit 0)

;;;
;;; Update site force field parameters by replacing them with the data
;;; from a table. The name of the site serves as a key to lookup table
;;; entries. The table is stored in guile/solutes.scm as a fake solute
;;; with all possible sites listed once.
;;;
(define (update-sites solute)
  (let* ((table                     ; fake solute with site-parameters
          (solute-sites (find-solute "bgy3d")))
         (update-one
          (lambda (site)
            ;; (pretty-print site)
            (let* ((name (site-name site))
                   (table-site (assoc name table)))
              (make-site name                          ; original
                         (site-position site)          ; original
                         (site-sigma table-site)       ; from table
                         (site-epsilon table-site)     ; from table
                         (site-charge site))))))       ; original
    ;; (pretty-print solutes)
    ;; (pretty-print table)
    (make-solute (solute-name solute)
                 (map update-one
                      (solute-sites solute)))))

;;;
;;; Will  not   work  with   distributed  vectors.   The   problem  is
;;; bgy3d-vec-ref (VecGetValues) does not either.
;;;
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

;;;
;;; BGY3d code  operates in  angstroms, QM codes  use atomic  units by
;;; convention:
;;;
(define (bohr->angstrom x) (* x 0.52917706))
(define (angstrom->bohr x) (/ x 0.52917706))

;;;
;;; Brain-dead implementation of cubic root:
;;;
(define (cubic-root n)
  (let loop ((x 0))
    (if (< n (* x x x))
        (- x 1)
        (loop (+ 1 x)))))
;;;
;;; This writes a GAMESS-UK punch file to the current output port. See
;;; interfaces/filepunch.py in CCP1 GUI repo:
;;;
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
               (S (* L (/ (- n 2) n)))) ; off-by-one on a 2L interval
          (header '((block . data) (records . 0)))
          (header '((block . grid_title) (records . 1)))
          (format #t "Solvent site hole density ~a\n" vec-id)
          ;;
          ;; Only  the number  of points "n"  is being  interpreted by
          ;; CCP1 GUI:
          ;;
          (header '((block . grid_axes) (records . 3)))
          (format #t "~a   0.000000 ~a 0 au xaxis\n" n (+ L S))
          (format #t "~a   0.000000 ~a 0 au yaxis\n" n (+ L S))
          (format #t "~a   0.000000 ~a 0 au zaxis\n" n (+ L S))
          ;;
          ;; The grid mapping block consists of six triples of floats:
          ;;
          ;;   O    A
          ;;   .    B
          ;;   .    C
          ;;
          ;; with  positions marked by dot apparently  ignored by CCP1
          ;; GUI.  All four points,  O, A, B,  and C, correspond  to a
          ;; grid point  of a regular grid at  four corners. Beware of
          ;; off-by-one errors when  converting periodic grids to this
          ;; representation:  if, say, the  left corner is at  -L then
          ;; its  periodic  image  to  the  right  is  at +L  but  the
          ;; rightmost point of the unique  grid portion is at S = L -
          ;; h, with h = 2L/n:
          ;;
          (header '((block . grid_mapping) (records . 3)))
          (format #t "~a ~a ~a  ~a ~a ~a\n" (- L) (- L) (- L) (+ S) (- L) (- L))
          (format #t "~a ~a ~a  ~a ~a ~a\n" (- L) (- L) (- L) (- L) (+ S) (- L))
          (format #t "~a ~a ~a  ~a ~a ~a\n" (- L) (- L) (- L) (- L) (- L) (+ S))
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
(define (bgy3d-run name sites . rest)
  "To be called from QM code."
  (let ((settings bgy3d-settings)
        (solute (list name sites))    ; incomplete solute descr
        (funptr (if (null? rest)      ; integer function pointer
                    0                 ; only present if called from PG
                    (first rest))))
    (if (equal? name "bgy3d")                ; only when called by PG
        (set! solute (update-sites solute))) ; update force-field params
    ;;
    ;; Extend settings by an  entry with the funciton pointer that can
    ;; be used to compute additional solute charge density:
    ;;
    (set! settings (acons 'qm-density   ; key
                          funptr        ; value
                          settings))    ; alist
    ;; Print on master only:
    (if (zero? (bgy3d-rank))
        (begin
          (pretty-print solute)
          (pretty-print settings)))
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
      (let ((path (if (zero? (bgy3d-rank))
                      "plot.pun"
                      "/dev/null")))    ; discard output of slaves
        (with-output-to-file path
          (lambda () (write-punch-file solute g1 settings))))
      ;;
      ;; Dont forget to destroy them after use:
      ;;
      (map bgy3d-vec-destroy g1))))

;;;
;;; FIXME: at  the moment this  function only emulates the  minimum of
;;; the   functionality   of   the   original  executable.   The   new
;;; functionality is in flux.
;;;
(define (old-main argv)
  (let ((opts (getopt-long argv option-spec)))
    ;; (pretty-print opts)
    (cond
     ;;
     ;; Pure solvent:
     ;;
     ((option-ref opts 'BGY2Site #f)
      (bgy3d-run-solvent '()))          ; Use defaults and Petsc env
     ;;
     ;; Solute  with solvent.  Note  that at variance with  the legacy
     ;; code the function find-solute uses on-disk database of solutes
     ;; in ./solutes.scm and not  the compiled in set of (currently) 6
     ;; preset solutes from bgy3d-solutes.c:
     ;;
     ((option-ref opts 'BGYM2Site #f)
      (let ((name (option-ref opts 'solute "hydrogen chloride")))
        (let ((g1 (bgy3d-run-solute (find-solute name)
                                    '()))) ; Use defaults and Petsc env
          (map bgy3d-vec-save (list "g0.bin" "g1.bin") g1)
          (map bgy3d-vec-destroy g1)))) ; dont forget to destroy them
     ;;
     ;; Fall through to the new variant:
     ;;
     (else
      (new-main argv)))))

;;;
;;; Interpretes each argument as the name of the solute and runs first
;;; a pure solvent then a solute calculations:
;;;
(define (new-main argv)
  ;; (pretty-print argv)
  (let* ((run (lambda (name)
                (apply bgy3d-run (find-solute name)))))
    (map run (cdr argv))))



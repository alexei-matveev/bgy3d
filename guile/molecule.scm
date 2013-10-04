;;;
;;; Representation of molecular solvents and solutes and site-specific
;;; force field parameters.
;;;
(define-module (guile molecule)
  #:use-module (srfi srfi-1)            ; list manipulation
  #:use-module (srfi srfi-2)            ; and-let*
  #:use-module (ice-9 match)            ; match-lambda
  #:use-module (guile tinker)           ; tinker-table
  #:use-module (guile utils)            ; memoize
  #:export
  (make-molecule
   move-molecule
   molecule-name
   molecule-sites
   molecule-positions
   molecule-charge
   molecule-dipole
   find-molecule
   print-molecule/xyz
   from-ab
   from-re
   from-cc
   make-site
   move-site
   site-name
   site-position
   site-sigma
   site-epsilon
   site-charge
   site-x
   site-y
   site-z))

;;;
;;; ITC Calorie here. See also bgy3d.h:
;;;
(define (kj->kcal num)
  (/ num 4.1868))
(define (kcal->kj num)
  (* num 4.1868))

;;;
;;; Database  contains  entries in  non-native  units  such as  kJ/mol
;;; represented by an s-expression  like (kJ 0.3640).  Convert them to
;;; working units kcal & angstrom:
;;;
(define (eval-units ast)
  (if (not (pair? ast))
      ast
      (match ast
        (('A num)
         (eval-units num))              ; for completeness, default
        (('r0 num)
         (* (eval-units num)
            (expt 2 (/ 5 6))))          ; (r0 x) -> x * 2^(5/6)
        (('^2 num)                      ; (^2 x) -> x^2
         (let ((num (eval-units num)))
           (* num num)))
        (('kcal num)
         (eval-units num))              ; for completeness, default
        (('kJ num)
         (kj->kcal (eval-units num)))   ; (kJ 0.3640) -> 0.08694
        (ast
         (map eval-units ast)))))


;;;
;;; LJ-type VdW interaction may take several forms.
;;;
;;; σε form: 4ε [-(σ/r)^6 + (σ/r)^12]
;;; Rε form: ε [-2(R/r)^6 + (R/r)^12]
;;; AB form: -(A/r)^6 + (B/r)^12
;;; 6-12 form: -C6 / r^6 + C12 / r^12
;;;
(define (from-ab params)
  (let ((A (first params))
        (B (second params)))
    (list (/ (* B B) A)                 ; σ = B^2/A
          (/ (expt (/ A B) 12) 4))))    ; ε = (A/B)^12 / 4

(define (from-re params)
  (let ((r (first params))
        (e (second params)))
    (list (* r (expt 2 (/ 5 6)))        ; σ = r * 2^(5/6)
          e)))                          ; ε = ε

(define (from-cc params)
  (let ((c6 (first params))
        (c12 (second params)))
    (list (expt (/ c12 c6) (/ 1 6))     ; σ = (c12/c6)^(1/6)
          (/ (* c6 c6) (* 4 c12)))))    ; ε = c6^2/4c12


;;;
;;; Find a solute/solvent in a database or die:
;;;
(define (find-molecule name)
  (and-let* ((solutes (slurp/cached (find-file "guile/solutes.scm")))
             (molecule (or (assoc name solutes)
                           (error "Not in the list:"
                                  name
                                  (map first solutes)))))
    (tabulate-ff (eval-units molecule))))

;; (set! find-molecule (memoize find-molecule))

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

;;;
;;; Premature optimization  and potential memory hog  here. Cache file
;;; contents in memory:
;;;
(define slurp/cached (memoize slurp))

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
(define (make-molecule name sites)
  (list name sites))

(define (molecule-name solute) (first solute))
(define (molecule-sites solute) (second solute))
(define (molecule-positions solute)
  (map site-position (molecule-sites solute)))

(define (move-molecule solute positions)
  (let ((name (molecule-name solute))
        (sites (molecule-sites solute)))
    (make-molecule name (map move-site sites positions))))

;;;
;;; Return  the   force  field   parameter  table  in   each  molecule
;;; description.  FIXME:  A missing table  is interpreted as  an empty
;;; table:
;;;
(define (molecule-ff-tab solute)
  (match solute
    ((name sites table) table) ; return the third entry, if there is such
    (_ '())))                  ; otherwise return an empty list

(define (make-site name position sigma epsilon charge)
  (list name position sigma epsilon charge))

(define (site-name site) (first site))
(define (site-position site) (second site))
(define (site-sigma site) (third site))
(define (site-epsilon site) (fourth site))
(define (site-charge site) (fifth site))

(define (move-site site new-position)
  (match site
    ((name old-position sigma epsilon charge)
     (make-site name new-position sigma epsilon charge))))

(define (site-x site) (first (site-position site)))
(define (site-y site) (second (site-position site)))
(define (site-z site) (third (site-position site)))

;;;
;;; This extracts the  charge field from every site  and sums them up.
;;; The  molecule needs  to be  in  the "canonical"  format, with  the
;;; force-field parameters expanded to numbers:
;;;
(define (molecule-charge mol)
  (apply + (map site-charge (molecule-sites mol))))

;;;
;;; Dipole moment. Same constraints as for molecule-charge:
;;;
(define (molecule-dipole mol)
  (define (dot a b)
    (apply + (map * a b)))
  (let* ((sites (molecule-sites mol))
         (q (map site-charge sites))
         (x (map site-x sites))
         (y (map site-y sites))
         (z (map site-z sites)))
    (list (dot q x)
          (dot q y)
          (dot q z))))

(define (print-molecule/xyz solute)
  (let ((name   (molecule-name solute))
        (sites  (molecule-sites solute)))
    (format #t "~a\n" (length sites))
    (format #t "# ~a\n" name)
    (for-each (lambda (site)
                (format #t "~a ~a ~a ~a\n"
                        (site-name site)
                        (site-x site)
                        (site-y site)
                        (site-z site)))
              sites)))

;; (for-each print-molecule/xyz (slurp (find-file "guile/solutes.scm")))
;; (exit 0)

;;;
;;; Example to load the ff database as a list, which could be obtained
;;; by tinker-test module.  This is  for test, ff parameters appear in
;;; each line
;;;
(define *tinker-ff-parameter-file*
  "guile/oplsaa-test")

(define (load-ff-file file)
  (let* ((contents (slurp/cached (find-file file)))
         (force-field (assoc "oplsaa" contents))) ; FIXME: literal here
    ;;
    ;; Force-field name in CAR  position is irrelevant for the rest of
    ;; the code, return only the list of entries:
    ;;
    (cdr force-field)))

;;;
;;; Get particular column from a data row as stored in the force-field
;;; data file.   Beware of the  fixed order of "sigma",  "epsilon" and
;;; "charge" parameters:
;;;
;;;   (82 HC "Alkane H-C" 1.008 2.5 0.03 0.06)
;;;
;;; or the longer version
;;;
;;;   (157 CT "Benzyl Alcohol -CH2OH" 12.011 3.5 0.066 0.2 6 13 4)
;;;
;;; with each field having the following meaning:
;;;
;;;   (type symbol descr weight sigma epsilon charge atnum class ligand)
;;;
;;; See e.g. ./oplsaa-test or ./tinker.scm for details.
;;;
(define (ff-symbol row) (second row))
(define (ff-sigma row) (fifth row))
(define (ff-epsilon row) (sixth row))
(define (ff-charge row) (seventh row))

;;;
;;; This will be called from outside, augment each site in solute with
;;; force field information (if missing):
;;;
(define (tabulate-ff solute)
  (let ((ff-tab (molecule-ff-tab solute))
        ;; (whole-file (tinker-table 'oplsaa))
        (whole-file (load-ff-file *tinker-ff-parameter-file*)))
    ;;
    ;; This  will derive the  atom type 77  from a name "CT3"  using a
    ;; list of entries (("CT3" 77) ...)
    ;;
    (define (ff-type name)
      (second (or (assoc name ff-tab)
                  (error "Not in the table:" name))))
    ;;
    ;; This  will return  a row from  the parameter file  matching the
    ;; given atom type:
    ;;
    (define (ff-row type)
      (or (assoc type whole-file)
          (error "Not in the file:" type)))

    ;;
    ;; Returns  a list of  non-bonding force field  parameters: (sigma
    ;; epsilon charge)
    ;;
    (define (non-bonding name)
      (let ((row (ff-row (ff-type name))))
        ;; (if (not (equal? name (symbol->string (ff-symbol row))))
        ;;     (pk "WARNING: do not match" name row))
        (list (ff-sigma row)
              (ff-epsilon row)
              (ff-charge row))))
    ;;
    ;; Build new sites and use them to construct a new molecule:
    ;;
    (let ((new-sites
           (map (lambda (site)
                  (match site
                    ((name position ('c6/c12 . cc) charge) ; -> derive σε from C6 & C12
                     (let* ((ff (from-cc cc)) ; cc is a list of two numbers
                            (sigma (first ff))
                            (epsilon (second ff)))
                       (make-site name position sigma epsilon charge)))
                    ((name position sigma epsilon charge) ; canonical form
                     ;; -> keep a valid site as is
                     site)
                    ((name position)   ; ff is missing
                     ;; -> append force field params
                     (apply make-site name position (non-bonding name)))))
                (molecule-sites solute))))
      (make-molecule (molecule-name solute) new-sites))))

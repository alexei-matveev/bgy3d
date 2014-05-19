;;;
;;; Copyright (c) 2013, 2014 Alexei Matveev
;;; Copyright (c) 2013 Bo Li
;;;
;;; Representation of molecular solvents and solutes and site-specific
;;; force field parameters.
;;;
(define-module (guile molecule)
  #:use-module (srfi srfi-1)            ; list manipulation
  #:use-module (srfi srfi-2)            ; and-let*
  #:use-module (srfi srfi-11)           ; let-values
  #:use-module (ice-9 match)            ; match-lambda
  #:use-module (ice-9 rdelim)           ; read-line
  #:use-module (ice-9 pretty-print)     ; pretty-print
  #:use-module (guile tinker)           ; tinker-table
  #:use-module (guile utils)            ; memoize
  #:use-module (guile atoms)            ; covalent-radius, canonical-name
  #:use-module (guile compat)           ; unless
  #:export
  (make-molecule
   move-molecule
   molecule-name
   molecule-sites
   molecule-positions
   molecule-charge
   molecule-dipole
   molecule-species
   molecule-self-energy
   scan-self-energy
   find-molecule
   print-molecule/xyz
   read-xyz
   kj->kcal
   kcal->kj
   kelvin->kcal
   kcal->kelvin
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
   update-param
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

(define (kelvin->kcal T)
  (let ((kboltzman (/ 8.3144621 (kcal->kj 1000))))
    (* kboltzman T)))
(define (kcal->kelvin E)
  (let ((kboltzman (/ 8.3144621 (kcal->kj 1000))))
    (/ E kboltzman)))

;;;
;;; Database  contains  entries in  non-native  units  such as  kJ/mol
;;; represented by an s-expression  like (kJ 0.3640).  Convert them to
;;; working  units kcal  &  angstrom. FIXME:  mutually recursive  with
;;; find-entry!
;;;
(define (expand-forms ast)
  (if (not (pair? ast))
      ast
      (match ast
        (('A num)
         (expand-forms num))            ; for completeness, default
        (('r0 num)
         (* (expand-forms num)
            (expt 2 5/6)))              ; (r0 x) -> x * 2^(5/6)
        (('^2 num)                      ; (^2 x) -> x^2
         (let ((num (expand-forms num)))
           (* num num)))
        (('kcal num)
         (expand-forms num))              ; for completeness, default
        (('kJ num)
         (kj->kcal (expand-forms num))) ; (kJ 0.3640) -> 0.08694
        (('geometry name)               ; Do not copy FF params
         (let ((sites (molecule-sites (find-entry name)))) ; Recursion here!
           (map (lambda (s)
                  (make-site (site-name s)
                             (site-position s)))
                sites)))
        (ast
         (map expand-forms ast)))))


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
    (list (* r (expt 2 5/6))            ; σ = r * 2^(5/6)
          e)))                          ; ε = ε

(define (from-cc params)
  (let ((c6 (first params))
        (c12 (second params)))
    (list (expt (/ c12 c6) 1/6)         ; σ = (c12/c6)^(1/6)
          (/ (* c6 c6) (* 4 c12)))))    ; ε = c6^2/4c12


;;;
;;; Find  a  solute/solvent in  a  database  or  die. FIXME:  mutually
;;; recursive with expand-forms!
;;;
(define (find-entry name)
  (and-let* ((solutes (slurp/cached (find-file "guile/solutes.scm")))
             (molecule (or (assoc name solutes)
                           (error "Not in the list:"
                                  name
                                  (map first solutes)))))
            (expand-forms molecule)))

(define (find-molecule name)
  (tabulate-ff (find-entry name)))

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

(define molecule-name first)
(define molecule-sites second)
;;; There  is also  an *optional*  third  field in  the molecule,  see
;;; molecule-force-field.

(define (molecule-positions solute)
  (map site-position (molecule-sites solute)))

(define (move-molecule solute positions)
  (let ((name (molecule-name solute))
        (sites (molecule-sites solute)))
    (make-molecule name (map move-site sites positions))))

;;;
;;; Return  the   force  field   parameter  table  in   each  molecule
;;; description.  FIXME:  A missing table  is interpreted as  an empty
;;; table that refers to an empty force field:
;;;
(define (molecule-force-field solute)
  (or (assoc 'force-field (cdr solute))
      '(force-field ())))               ; fake empty force field table

;;;
;;; Two  legal ways to  call make-site:  with and  without force-field
;;; parameters. FIXME:  This is road to  hell, a site  created by this
;;; constructor when passed to force-field accessors may fail ...
;;;
(define make-site
  (case-lambda
   ((name position sigma epsilon charge)
    (list name position sigma epsilon charge))
   ((name position)
    (list name position))))

(define site-name first)
(define site-position second)
(define site-sigma third)
(define site-epsilon fourth)
(define site-charge fifth)

;;;
;;; update sigma, epsilon and charge with new values
;;;

(define (find-site sites site-name)
  ;; (remove-if-not (lambda (x) (equal? (car x) site-name)) sites))
  (let ((site (assoc site-name sites)))
    (if (not site)
      (error "Not in sites" site-name))))

(define (update-param molecule target-site param value)
  (let ((name (molecule-name molecule))
        (sites (molecule-sites molecule)))

    (define (update-site site)
      (match param
             ("sigma"
              (if (equal? target-site (site-name site))
                (update-sigma site value)
                (update-sigma site (site-sigma site))))
             ("epsilon"
              (if (equal? target-site (site-name site))
                (update-epsilon site value)
                (update-epsilon site (site-epsilon site))))
             ("charge"
              (if (equal? target-site (site-name site))
                (update-charge site value)
                (update-charge site (site-charge site))))))

    (find-site sites target-site)
    (make-molecule name (map update-site sites))))

(define (update-sigma site new-sigma)
  (match site
    ((name position old-sigma epsilon charge)
     (make-site name position new-sigma epsilon charge))))

(define (update-epsilon site new-epsilon)
  (match site
    ((name position sigma old-epsilon charge)
     (make-site name position sigma new-epsilon charge))))

(define (update-charge site new-charge)
  (match site
    ((name position sigma epsilon old-charge)
     (make-site name position sigma epsilon new-charge))))

(define (move-site site new-position)
  (match site
    ((name old-position sigma epsilon charge)
     (make-site name new-position sigma epsilon charge))))

(define (site-x site) (first (site-position site)))
(define (site-y site) (second (site-position site)))
(define (site-z site) (third (site-position site)))

(define (dot a b)
  (apply + (map * a b)))

(define (site-distance a b)
  (let ((xa (site-position a))
        (xb (site-position b)))
    (let ((d (map - xa xb)))
      (sqrt (dot d d)))))

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
  (let* ((sites (molecule-sites mol))
         (q (map site-charge sites))
         (x (map site-x sites))
         (y (map site-y sites))
         (z (map site-z sites)))
    (list (dot q x)
          (dot q y)
          (dot q z))))

(define (molecule-species solute scale)
  ;;
  ;; Estimate the length of a typical bond between two such sites:
  ;;
  (define (bond-length a b)
    (+ (covalent-radius (site-name a))
       (covalent-radius (site-name b))))
  ;;
  ;; Comparator for two sites. Returns true if sites are closer than
  ;; scaled estmate for a typical bond between two such atoms. It
  ;; should not be used to compare the keys in a dictionary.  Because
  ;; this is NOT an equivalence relation, not even transitive.
  ;;
  (define (close? a b)
    (let ((distance (site-distance a b))
          (estimate (bond-length a b)))
      (< distance (* scale estimate))))
  ;;
  ;; A  site is connected to  a species if it  is close to  any of the
  ;; sites of that species:
  ;;
  (define (connected? a species)
    (any (lambda (b) (close? a b))
         species))
  ;;
  ;; Group sites into species recursively. Each new site is connected
  ;; to zero or more species and thus makes another, eventually
  ;; larger, species.
  ;;
  (define (classify sites)
    (let loop ((sites sites)
               (species '()))
      (if (null? sites)
          (reverse species) ; reverse is optional for a set
          (let ((a (car sites)))
            ;; Partition the species into ones connected to site "a"
            ;; and other species.
            (let-values (((conn other) (partition (lambda (s)
                                                    (connected? a s))
                                                  species)))
              ;; A new species will contain the site "a" and zero or
              ;; more sites of the connected species:
              (loop (cdr sites)
                    (cons (cons a (concatenate conn))
                          other)))))))
  ;;
  ;; Assign a site a numeric ID of the species or #f if not found:
  ;;
  (define (lookup a species)
    (let loop ((n 0)
               (species species))
      (if (null? species)
          #f                            ; not found
          (if (memq a (car species))
              n
              (loop (+ 1 n) (cdr species))))))
  ;;
  ;; Return  a list of numeric  IDs. This list is  passed to C/Fortran
  ;; code occasionally:
  ;;
  (let* ((sites (molecule-sites solute))
         (species (classify sites)))
    (map (lambda (a)
           (lookup a species))
         sites)))



(define (print-molecule/xyz solute)
  (let ((name   (molecule-name solute))
        (sites  (molecule-sites solute)))
    (format #t "~a\n" (length sites))
    (format #t "# ~a\n" name)
    (for-each (lambda (site)
                (format #t "~a ~a ~a ~a\n"
                        (canonical-name (site-name site))
                        (site-x site)
                        (site-y site)
                        (site-z site)))
              sites)))

;; (for-each print-molecule/xyz (slurp (find-file "guile/solutes.scm")))
;; (exit 0)

(define (energy a b)
  ;;
  (define (lj r)
    (let ((sr6 (expt (/ 1 r) 6)))
      (* 4 sr6 (- sr6 1))))
  ;;
  (define (combine-sigmas sa sb)
    (/ (+ sa sb) 2))
  ;;
  (define (combine-epsilons ea eb)
    (sqrt (* ea eb)))
  ;;
  (define EPSILON0INV 331.84164)       ; FIXME: literal from units.f90
  ;;
  (let ((sa (site-sigma a))
        (ea (site-epsilon a))
        (qa (site-charge a))
        (sb (site-sigma b))
        (eb (site-epsilon b))
        (qb (site-charge b)))
    (let ((sab (combine-sigmas sa sb))
          (eab (combine-epsilons ea eb))
          (qab (* qa qb))
          (rab (site-distance a b)))
      (let ((e-coul (* EPSILON0INV (/ qab rab)))
            (e-lj (if (zero? sab)
                      0.0
                      (* eab (lj (/ rab sab))))))
        (+ e-coul e-lj)))))


(define (molecule-self-energy m species)
  (let ((sites (molecule-sites m)))
    (let lp1 ((sites sites)
              (species species)
              (e 0.0))
      (if (null? sites)
          e
          (let ((a (car sites))
                (a-species (car species))
                (sites (cdr sites))
                (species (cdr species)))
            (lp1 sites
                 species
                 (let lp2 ((sites sites)
                           (species species)
                           (e e))
                   (if (null? sites)
                       e
                       (let ((b (car sites))
                             (b-species (car species))
                             (sites (cdr sites))
                             (species (cdr species)))
                         (lp2 sites
                              species
                              (+ e (if (equal? a-species b-species)
                                       0.0
                                       (energy a b)))))))))))))


(define (interpolate x a b)
  (let ((wb x)
        (wa (- 1 x)))
    (let go ((a a)
             (b b))
      (cond
       ((and (null? a) (null? b))
        '())
       ((and (pair? a) (pair? b))
        (cons (go (car a)
                  (car b))
              (go (cdr a)
                  (cdr b))))
       (else
        (+ (* wa a) (* wb b)))))))


(define (scan-self-energy mol x0 x1)
  (let ((species (molecule-species mol 1.0))
        (n 20))
    (map (lambda (i)
           (let ((mol' (move-molecule mol
                                      (interpolate (/ i (- n 1))
                                                   x0
                                                   x1))))
             (molecule-self-energy mol' species)))
         (iota n))))


;;;
;;; This reads  the current  input port and  returns a  structure that
;;; resembles a molecule:
;;;
;;;  ("Comment line goes here"
;;;    (("AtomName1" (x y z)) ...))
;;;
(define (read-xyz)
  ;; Here let*  is used to force evaluation  order of expressions that
  ;; have side-effects:
  (let* ((natoms (read))                 ; number of atoms
         (rest-of-line (read-line))      ; discard the rest of the line
         (comment-line (read-line)))     ; any text
    (let loop ((n natoms)
               (acc '()))
      (if (zero? n)
          (make-molecule comment-line (reverse acc))
          (let* ((name (symbol->string (read))) ; yes, let*
                 (x (read))
                 (y (read))
                 (z (read)))
            (loop (- n 1)
                  (cons (make-site name (list x y z)) acc)))))))

;; (let ((m (with-input-from-file "a1.xyz" read-xyz)))
;;   (pretty-print m)
;;   (sorted-distances m))
;; (exit 0)

;;;
;;; Example to load the ff database as a list, which could be obtained
;;; by tinker-test module.  This is  for test, ff parameters appear in
;;; each line
;;;
(define *tinker-ff-parameter-file*
  "guile/force-fields.scm")

(define (find-force-field force-field)
  (let ((contents (slurp/cached (find-file *tinker-ff-parameter-file*))))
    ;;
    ;; Force-field name in CAR  position is irrelevant for the rest of
    ;; the code, return only the list of entries:
    ;;
    (assoc-ref contents force-field)))

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
;;; See e.g. ./force-fields.scm or ./tinker.scm for details.
;;;
(define ff-symbol second)
(define ff-sigma fifth)
(define ff-epsilon sixth)
(define ff-charge seventh)

;;;
;;; This will be called from outside, augment each site in solute with
;;; force field information (if missing):
;;;
(define (tabulate-ff solute)
  ;;
  ;; The structure of FF form:
  ;;
  ;;  (force-field (ff-sym1 ff-sym2 ...)
  ;;    (site-name1 site-sym1) ...)
  ;;
  ;; The car position is ignored. The CADR position is a list of force
  ;; field symbols  to combine. CDDR is the  translation table of site
  ;; names  to symbols/numbers used in  FF tables. FIXME:  There is no
  ;; checking if a site symbols is repeated in multiple rows.
  ;;
  (let* ((force-field-form (molecule-force-field solute))
         (force-field-symbols (cadr force-field-form))
         (translation-table (cddr force-field-form))
         (rows (concatenate (map find-force-field force-field-symbols))))
    ;;
    ;; This  will derive the  atom type 77  from a name "CT3"  using a
    ;; list of entries (("CT3" 77) ...) in the FF table.
    ;;
    (define (ff-type name)
      (second (or (assoc name translation-table)
                  (error "Not in the table:" name))))
    ;;
    ;; This  will return  a row from  the parameter file  matching the
    ;; given atom type:
    ;;
    (define (ff-row type)
      (or (assoc type rows)
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
    ;; This  function adds missing  foce-field parameters to  the site
    ;; description:
    ;;
    (define (new-site site)
      (match site
        ;; 1) Canonical form -> keep a valid site as is
        ((name position sigma epsilon charge)
         site)
        ;; 2) FF is missing -> append force field params
        ((name position)
         (apply make-site name position (non-bonding name)))
        ;; 3) A form  instead of literal sigma and epsilon
        ;; -> derive σε from C6 & C12
        ((name position ('c6/c12 . cc) charge)
         (let* ((ff (from-cc cc))       ; cc is a list of two numbers
                (sigma (first ff))
                (epsilon (second ff)))
           (make-site name position sigma epsilon charge)))))
    ;;
    ;; FIXME: This check is so far down because of the "define"s above
    ;; that have to be put at the beginning of the scope.
    ;;
    (unless rows
      (error "Did not find force-field:" force-field-symbols))
    ;;
    ;; Build new sites and use them to construct a new molecule:
    ;;
    (let ((new-sites (map new-site (molecule-sites solute))))
      (make-molecule (molecule-name solute) new-sites))))

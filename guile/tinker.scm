;;;
;;; Tools to handle TINKER input such as parameter files.
;;;
(define-module (guile tinker)
  #:use-module (srfi srfi-1)            ; list manipulation
  #:use-module (ice-9 rdelim)           ; read-line
  #:use-module (ice-9 match)            ; match-lambda
  #:use-module (ice-9 pretty-print)     ; pretty-print
  #:export
  (tinker-slurp
   tinker-table))

;;;
;;; Reads a file and returns a list of lines:
;;;
(define (lines file)
  (with-input-from-file file
    (lambda ()
      (let loop ((acc (list)))
        (let ((line (read-line)))
          (if (eof-object? line)
              (reverse acc)
              (loop (cons line acc))))))))

;;;
;;; Destructively  replaces character  '#' with  ';'. At  the  time of
;;; writing semicolons were used only in #-commented sections.
;;;
;;; FIXME: unfortunately there is at  least one use of the #-character
;;; in  atomic class  symbol ---  C#.   See atom  types 897  & 898  in
;;; oplsaa.prm and smoothaa.prm. It happens  that this C# atom type is
;;; the only  place where the  char # is  not preceeded by a  space or
;;; another # (cf.  grep "[^ #]#" *.prm).
;;;
(define (rewrite! str)
  (let ((n (string-length str)))
    (let loop ((i 0)                    ; pointer into the string
               (f #t)) ; should we interprete eventual # as a comment sign?
      (if (< i n)
          (let ((c (string-ref str i))
                (i++ (+ 1 i)))
            (if (and f (equal? c #\#))
                (begin
                  (string-set! str i #\;)
                  (loop i++ #t))
                (loop i++ (char-whitespace? c))))))))


;;;
;;; (Incomplete) list of TINKER keywords used in parameter files:
;;;
(define *keywords*
  '(
    angle
    atom
    biotype
    bond
    charge
    chg-14-scale
    dielectric
    electric
    epsilonrule
    forcefield
    imptors
    imptorunit
    opbend
    pitors
    polarization
    polarize
    radiusrule
    radiussize
    radiustype
    strbnd
    torsion
    torsionunit
    vdw
    vdw-14-scale
    vdwindex
    vdwtype
    ;; added in the second iteration:
    angang
    dipole
    dipole4
    dipole5
    electneg
    hbond
    improper
    mmffangle
    mmffarom
    mmffbci
    mmffbond
    mmffbonder
    mmffcovrad
    mmffdefstbn
    mmffequiv
    mmffopbend
    mmffpbci
    mmffprop
    mmffstrbnd
    mmfftorsion
    mmffvdw
    multipole
    piatom
    pibond
    pibond5
    strtors
    tortors
    ureybrad
    ))


;;;
;;; Regular  expression  for  valid   lines  starting  with  a  TINKER
;;; keyword. The real  problem is that the parameter  files include an
;;; uncommented     section     with     a     human-readable     text
;;; (Literature).  Everything  that  does  not match  this  regexp  is
;;; considered such a text.
;;;
(define *rx*
  (make-regexp
   (string-append
    "^("
    (string-join (map symbol->string *keywords*)
                 "|")
    ")")))


;;;
;;; "Read" all fields turning them into Scheme objects.
;;;
(define (lex line)
  ;; (pk line)
  (with-input-from-string line
    (lambda ()
      (let loop ((acc (list)))
        (let ((x (read)))
          (if (eof-object? x)
              (reverse acc)
              (loop (cons x acc))))))))


;;;
;;; A hack to either tokenize a valid line or leave the prosa as (text
;;; "Prosa"):
;;;
(define (tinker-slurp file)
  (let ((text (lines file)))
    (for-each rewrite! text)
    (map (lambda (line)
           (if (regexp-exec *rx* line)  ; if starts with a valid kw ...
               (lex line)               ; then lex the line ..
               (list 'text line)))      ; ... otherwise dont do that.
         text)))


;;;
;;; Extract entries  from the parameter  file contents as  returned by
;;; tinker-slurp having  specific symbol  in car position.   Strip the
;;; redundant CAR position in result.
;;;
(define (filter-contents symbol contents)
  (map cdr ; strip redundant symbol in car position common for all entries
       (filter (lambda (row)
                 (equal? symbol (first row)))
               contents)))


;;;
;;; Make a table by combining  atom-, vdw-, and charge entries for all
;;; atom types:
;;;
(define (make-table contents)
  (let ((atm-tab (filter-contents 'atom contents))
        (vdw-tab (filter-contents 'vdw contents))
        (chg-tab (filter-contents 'charge contents)))
    (map (match-lambda*
          ;;
          ;; Case 1. Type and class:
          ;;
          ;; Argument list is a 3-list having the following structure:
          ;;
          (((type class symbol descr atnum weight ligand)
            (type sigma epsilon)
            (type charge))
           ;; ->
           (list
            type                        ; 1
            symbol                      ; 2
            descr                       ; 3
            weight                      ; 4
            sigma                       ; 5
            epsilon                     ; 6
            charge                      ; 7
            atnum                       ; 8
            class                       ; 9
            ligand                      ; 10
            ))
          ;;
          ;; Case 2. Type but no class:
          ;;
          (((type symbol descr atnum weight ligand)
            (type sigma epsilon)
            (type charge))
           ;; ->
           (list
            type                        ; 1
            symbol                      ; 2
            descr                       ; 3
            weight                      ; 4
            sigma                       ; 5
            epsilon                     ; 6
            charge                      ; 7
            atnum                       ; 8
            type                        ; 9 FIXME: type == class?
            ligand                      ; 10
            )))
         atm-tab
         vdw-tab
         chg-tab)))

;;;
;;; In general, atom class is not exactly in one-to-one correspondence
;;; with atom  symbol. Some TINKER  force fields not  even distinguish
;;; between atom type and atom  class.  However, for some force fields
;;; there is such a correspondence, namely for:
;;;
;;;   charmm19 charmm22cmap charmm22 hoch iwater oplsaa smoothaa water
;;;
;;; For AMBER force  fields there is just one  conflict see "Formyl H"
;;; entry that quotes 0 as  the atom class whereas all other hydrogens
;;; have class 29.
;;;
(define (verify-class/symbol-equivalence atom-table)
  (let ((classes (map second atom-table))
        (symbols (map third atom-table))
        (class->symbol/table (make-hash-table (length atom-table)))
        (symbol->class/table (make-hash-table (length atom-table))))
    (for-each
        (lambda (class symbol)
          (hash-set! class->symbol/table class symbol)
          (hash-set! symbol->class/table symbol class))
      classes
      symbols)
    (for-each
        (lambda (c s)
          (let ((s' (hash-ref class->symbol/table c))
                (c' (hash-ref symbol->class/table s)))
            (unless (and (equal? c c')
                         (equal? s s'))
                    (pretty-print (list "ERROR: Do not match" c c' s s')))))
      classes
      symbols)))

;;;
;;; These are distributed with TINKER:
;;;
(define *force-fields*
  '(amber94
    amber96
    amber98
    amber99
    amber99sb
    amoeba04
    amoeba09
    amoebabio09
    amoebapro04
    amoebapro13
    charmm19
    charmm22cmap
    charmm22
    dang
    hoch
    iwater
    ;;
    ;; These three have charge entries only for a few of atom
    ;; types. This breaks the simplified logic of make-table:
    ;;
    ;; mm2
    ;; mm3
    ;; mm3pro
    mmff
    oplsaal
    oplsaa
    oplsua
    smoothaa
    smoothua
    tiny
    water))


(define *tinker-parameter-dir*
  "/home/alexei/devel/tinker/params/")

(define (tinker-parameter-file force-field)
  (string-append *tinker-parameter-dir*
                 (symbol->string force-field)
                 ".prm"))


;;;
;;; Reads and parses a file to return a list of rows:
;;;
(define (tinker-table force-field)
  (make-table (tinker-slurp (tinker-parameter-file force-field))))


;;;
;;; Let-over-lambda  here.   Given a  function  (f  x)  The result  of
;;; (memoize f)  is a  function (f' x)  == (f  x) that will  cache all
;;; results.  FIXME: make it work for arbitrary number of arguments.
;;;
(define (memoize f)
  (let ((*cache* '()))                  ; empty cache
    (lambda (x)
      (cdr (or (assoc x *cache*)        ; first check the cache
               (let* ((y (f x))         ; otherwise invoke f
                      (p (cons x y)))   ; make dictionary pair
                 (set! *cache* (cons p *cache*)) ; cache new pair
                 p))))))                         ; and return it too

;;;
;;; The  new  tinker-table should  cache  the  results (hopefully  the
;;; force-field is not going to change between two calls):
;;;
(set! tinker-table (memoize tinker-table))


;;;
;;; Tools to handle TINKER input such as parameter files.
;;;
(define-module (guile tinker)
  #:use-module (srfi srfi-1)            ; list manipulation
  #:use-module (ice-9 rdelim)           ; read-line
  #:use-module (ice-9 pretty-print)     ; pretty-print
  #:export
  (tinker-slurp
   tinker-test))

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
;;; in atomic  class symbols C#.   See e.g.  atom  types 897 &  898 in
;;; oplsaa.prm.
;;;
(define (rewrite! str)
  (string-map! (lambda (c)
                 (case c
                   ((#\#) #\;)          ; # -> ;
                   (else c)))           ; leave the rest alone
               str))

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
;;; These are distributed with TINKER:
;;;
(define *files*
  '("amber94.prm"
    "amber96.prm"
    "amber98.prm"
    "amber99.prm"
    "amber99sb.prm"
    "amoeba04.prm"
    "amoeba09.prm"
    "amoebabio09.prm"
    "amoebapro04.prm"
    "amoebapro13.prm"
    "charmm19.prm"
    "charmm22cmap.prm"
    "charmm22.prm"
    "dang.prm"
    "hoch.prm"
    "iwater.prm"
    "mm2.prm"
    "mm3.prm"
    "mm3pro.prm"
    "mmff.prm"
    "oplsaal.prm"
    "oplsaa.prm"
    "oplsua.prm"
    "smoothaa.prm"
    "smoothua.prm"
    "tiny.prm"
    "water.prm"))


;;;
;;; Is used to test ideas:
;;;
(define (tinker-test)
  (for-each
   ;; (lambda (file)
   ;;   (let ((content (tinker-slurp (string-append "/home/alexei/devel/tinker/params/" file))))
   ;;     (for-each (lambda (row)
   ;;                 (display (second row))
   ;;                 (newline))
   ;;               (filter (lambda (row) (equal? (car row) 'text))
   ;;                       content)))
   (lambda (file)
     (with-output-to-file file
       (lambda ()
         (pretty-print
          (tinker-slurp (string-append "/home/alexei/devel/tinker/params/"
                                       file))))))
   *files*))

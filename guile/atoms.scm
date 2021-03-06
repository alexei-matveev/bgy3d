;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;; Atomic data
;;;
(define-module (guile atoms)
  #:use-module (srfi srfi-1)            ; list manipulation
  #:use-module (srfi srfi-2)            ; and-let*
  #:use-module (ice-9 pretty-print)
  #:export
  (canonical-name
   covalent-radius))

;;;
;;; So far only radii (copied from ParaTools).
;;;
(define *periodic-table*
  '(("H"  0.3200)
    ("He" 0.9300)
    ("Li" 1.2300)
    ("Be" 0.9000)
    ("B"  0.8200)
    ("C"  0.7700)
    ("N"  0.7500)
    ("O"  0.7300)
    ("F"  0.7200)
    ("Ne" 0.7100)
    ("Na" 1.5400)
    ("Mg" 1.3600)
    ("Al" 1.1800)
    ("Si" 1.1100)
    ("P"  1.0600)
    ("S"  1.0200)
    ("Cl" 0.9900)
    ("Ar" 0.9800)
    ("K"  2.0300)
    ("Ca" 1.7400)
    ("Sc" 1.4400)
    ("Ti" 1.3200)
    ("V"  1.2200)
    ("Cr" 1.1800)
    ("Mn" 1.1700)
    ("Fe" 1.1700)
    ("Co" 1.1600)
    ("Ni" 1.1500)
    ("Cu" 1.1700)
    ("Zn" 1.2500)
    ("Ga" 1.2600)
    ("Ge" 1.2200)
    ("As" 1.2000)
    ("Se" 1.1600)
    ("Br" 1.1400)
    ("Kr" 1.1200)
    ("Rb" 2.1600)
    ("Sr" 1.9100)
    ("Y"  1.6200)
    ("Zr" 1.4500)
    ("Nb" 1.3400)
    ("Mo" 1.3000)
    ("Tc" 1.2700)
    ("Ru" 1.2500)
    ("Rh" 1.2500)
    ("Pd" 1.2800)
    ("Ag" 1.3400)
    ("Cd" 1.4800)
    ("In" 1.4400)
    ("Sn" 1.4100)
    ("Sb" 1.4000)
    ("Te" 1.3600)
    ("I"  1.3300)
    ("Xe" 1.3100)
    ("Cs" 2.3500)
    ("Ba" 1.9800)
    ("La" 1.6900)
    ("Ce" 1.6500)
    ("Pr" 1.6500)
    ("Nd" 1.6400)
    ("Pm" 1.6300)
    ("Sm" 1.6200)
    ("Eu" 1.8500)
    ("Gd" 1.6100)
    ("Tb" 1.5900)
    ("Dy" 1.5900)
    ("Ho" 1.5800)
    ("Er" 1.5700)
    ("Tm" 1.5600)
    ("Yb" 1.5600)
    ("Lu" 1.5600)
    ("Hf" 1.4400)
    ("Ta" 1.3400)
    ("W"  1.3000)
    ("Re" 1.2800)
    ("Os" 1.2600)
    ("Ir" 1.2700)
    ("Pt" 1.3000)
    ("Au" 1.3400)
    ("Hg" 1.4900)
    ("Tl" 1.4800)
    ("Pb" 1.4700)
    ("Bi" 1.4600)
    ("Po" 1.4600)
    ("At" 1.4500)
    ("Rn" 1.0000)
    ("Fr" 1.0000)
    ("Ra" 1.0000)
    ("Ac" 1.0000)
    ("Th" 1.6500)
    ("Pa" 1.0000)
    ("U"  1.4200)
    ("Np" 1.0000)
    ("Pu" 1.0000)
    ("Am" 1.0000)
    ("Cm" 1.0000)
    ("Bk" 1.0000)
    ("Cf" 1.0000)
    ("Es" 1.0000)
    ("Fm" 0.8000)
    ("Md" 1.0000)
    ("No" 1.0000)
    ("Lr" 1.0000)
    ))

;;;
;;; List  of  site   names  for  each  atom  that   deviate  from  the
;;; default. Hm, not a good idea. I've seen H1-H8.
;;;
(define *site-names*
  '(("H" "HW")
    ("O" "OH" "OW" "OU")))

;;;
;;; Hash map from aliases to canonical names:
;;;
(define *canonical-names*
  (let ((hash (make-hash-table)))
    (for-each (lambda (row)
                (let ((name (car row))
                      (aliases (cdr row)))
                  (for-each (lambda (alias)
                              (hash-set! hash alias name))
                            aliases)))
              *site-names*)
    hash))

(define (canonical-name name)
  (or (hash-ref *canonical-names* name)
      name))

;;;
;;; FIXME: Solute  database contains atoms with "weird"  names such as
;;; "H8"  (yes, really). We  should not  make code  to guess  ... Self
;;; energy of such solutes will be wrong. These radii are only used to
;;; "guess topology" and  omit computing non-bonded interaction within
;;; a species.
;;;
(define (covalent-radius name)
  (let ((atom (canonical-name name)))
    (or (and-let* ((row (assoc-ref *periodic-table* atom)))
          (first row))
        0.0)))                          ; FIXME: a better solution?

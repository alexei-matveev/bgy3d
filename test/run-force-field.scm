;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;;  ../bgy3d -L ../ -s ./run-force-field.scm
;;;
(use-modules
 (guile bgy3d)          ; ccap0, ccap1, etc.
 (guile utils)          ; dfridr
 (ice-9 pretty-print))

(define RMIN 0.5)                       ; see bgy3d-solutes.c
(define h 1.0e-7)
(define RMIN+ (+ RMIN h))
(define RMIN- (- RMIN h))

(define (test name f0 f1)
  (let ((rs (list RMIN- RMIN RMIN+)))
    (pretty-print name)
    (for-each
     (lambda (r)
       (let* ((pair (dfridr f0 r h))    ; num diff
              (f1-num (car pair))
              (f1-err (cdr pair)))
         (pretty-print (list r (f0 r) (f1 r) f1-num))))
     rs)))

;;; FIXME: RMIN here and in bgy3d-solutes.c may diverge:
(unless (equal? (ccap0 RMIN+)
                (/ 1.0 RMIN+))
  (error "Capped Coulomb is wrong for r > RMIN"))

(when (equal? (ccap0 RMIN-)
              (/ 1.0 RMIN-))
  (error "Capped Coulomb is wrong for r < RMIN"))

(define primitives
  (list (list "lj" lj0 lj1)
        (list "ljcap" ljcap0 ljcap1)
        (list "cl" cl0 cl1)
        (list "cs" cs0 cs1)
        (list "cscap" cscap0 cscap1)
        (list "ccap" ccap0 ccap1)))

(for-each
 (lambda (args)
   (apply test args))
 primitives)

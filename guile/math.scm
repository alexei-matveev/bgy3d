;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
;;; Export a few functions from libm.
;;;
(define-module (guile math)
  #:use-module (system foreign)
  #:export
  (erf
   erfc
   j0))

(define libm (dynamic-link "libm"))

(define erf
  (pointer->procedure double
                      (dynamic-func "erf" libm)
                      (list double)))
(define erfc
  (pointer->procedure double
                      (dynamic-func "erfc" libm)
                      (list double)))
(define j0
  (pointer->procedure double
                      (dynamic-func "j0" libm)
                      (list double)))



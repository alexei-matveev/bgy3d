;;;
;;; Copyright (c) 2013 Alexei Matveev
;;;
(define-module (guile compat))

(cond-expand
 ((not guile-2)
  ;;
  ;; Code to backport selected features to 1.8 goes here:
  ;;
  (use-syntax (ice-9 syncase))

  ;;
  ;; Thanks to taylanub at #guile:
  ;;
  (define-syntax define-syntax-rule
    (syntax-rules ()
      ((define-syntax-rule (keyword . pattern) template)
       (define-syntax keyword
         (syntax-rules ()
           ((keyword . pattern) template))))))

  ;;
  ;; Not in the define-module form:
  ;;
  (export-syntax define-syntax-rule)

  ;;
  ;; When/uless were added in 2.0:
  ;;
  (define-syntax-rule (when test stmt stmt* ...)
    (if test (begin stmt stmt* ...)))

  (export-syntax when)

  (define-syntax-rule (unless test stmt stmt* ...)
    (if (not test) (begin stmt stmt* ...)))

  (export-syntax unless)

  ;;
  ;; Re-export all symbols from (ice-9 syncase). Credit to mark_weaver
  ;; at #guile. Otherwise the importer of this module, (guile compat),
  ;; still  needs to  explicitly import (ice-9  syncase) in  order for
  ;; define-syntax-rule to work.
  ;;
  (let ((names (module-map (lambda (name var) name)
                           (resolve-interface '(ice-9 syncase)))))
    (module-re-export! (current-module) names)))
 (else)
 ;;
 ;; Nothing here.
 ;;
 )

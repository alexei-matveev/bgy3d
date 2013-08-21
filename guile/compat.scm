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
  ;; Re-export all symbols from (ice-9 symcase). Credit to mark_weaver
  ;; at #guile. Otherwise the importer of this module, (guile compat),
  ;; still needs to explicitly import (ice-9 syncase) too in order for
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

#!/users/alexei/darcs/bgy3d/bgy3d -s
!#
;;;
;;; The usual trick  with #!/usr/bin/env bgy3d does not  work as guile
;;; needs a  -s argument to  work as an  interpeter.  Use the  path to
;;; your own bgy3d interpreter, if necessary.
;;;
;;; To find the  custom modules (and solute tables)  extend the module
;;; search path to contain a path to the BGY3d repository:
;;;
(set! %load-path (cons "/users/alexei/darcs/bgy3d" %load-path))

(use-modules ((guile bgy3d)
              #:select
              (old-main
               new-main)))

;;;
;;; We are  trying to  emulate behaviour of  old executable  unless we
;;; find a better interface:
;;;
(let ((argv (command-line)))
  (if (string-prefix? "--" (cadr argv)) ; what if no options supplied?
      (old-main argv)
      (new-main (cdr argv)))) ; otherwise the first argument should be
                              ; a subcommand


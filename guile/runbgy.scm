#!/users/alexei/darcs/bgy3d/bgy3d -s
!#
;;;
;;; Copyright (c) 2013 Alexei Matveev
;;;
;;; The usual trick  with #!/usr/bin/env bgy3d does not  work as guile
;;; needs a  -s argument to  work as an  interpeter.  Use the  path to
;;; your own bgy3d interpreter, if necessary.
;;;
;;; To find the  custom modules (and solute tables)  extend the module
;;; search path to contain a path to the BGY3d repository. This should
;;; work with  both 1.8  and 2.0.  Use  absolute paths  when extending
;;; %load-path --- there is no globbing for ~/ at compile time.
;;;
(cond-expand
 ((not guile-2) (use-modules (ice-9 syncase))) ; eval-when for 1.8
 (else))                                ; nothing

(eval-when
 (eval load compile)      ; extend %load-path at compile time with 2.0
 (set! %load-path (cons "/users/alexei/darcs/bgy3d" %load-path)))

(use-modules (guile bgy3d))             ; bgy3d-main


(bgy3d-main (command-line))


;;; -*- mode: scheme; -*- vim: set syntax=scheme:
;;;
;;; Copyright (c) 2014 Alexei Matveev
;;;
(define ss '("0w"
             "4w" "5w" "6w"
             ))
(define cs '("KH"
             "PSE2" "PSE3"
             ))

(define (pg-input nw)
  (string-append nw ",gp.scm"))

(define (pg-output nw)
  (string-append "o." nw ",gp"))

(define (make-cmd nw closure)
  (let ((closure-flag (string-append "--closure=" closure))
        (solute-name (string-append "uranyl, " nw ", pcm"))
        (input (pg-input nw)))
    `("mpiexec" "-n" "16" "/users/alexei/darcs/bgy3d-wheezy/test-qm/uranyl/sym/runqmmm"
      "--solvent" "water, PR-SPC/E"
      "--solute" ,solute-name
      "--norm-tol=1e-14"
      "--dielectric=78.4"
      "--rho=0.0333295"
      "--beta=1.6889"
      "--L=20"
      "--N=256"
      "--hnc"
      ,closure-flag
      ,input)))

(define (gxfile nw closure)
  (string-append nw "," closure ".gx"))

(define (sh* . args)
  (apply system* args))

(define (run nw clo)
  (sh* "rm" "-f" "hesse.dat" "saved_scfstate.dat")
  (sh* "rm" "-rf" (pg-output nw))
  (sh* "cp" (gxfile nw clo) "gxfile")
  (apply sh* (make-cmd nw clo))
  (sh* "mv" "gxfile" (gxfile nw clo))
  (sh* "mv" (pg-output nw) (string-append "o." nw "," clo)))

(for-each (lambda (clo)
            (for-each (lambda (nw)
                        (run nw "KH"))
                      ss))
          cs)


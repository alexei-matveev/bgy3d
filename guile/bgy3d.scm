;;;
;;; Scheme  interface  to  BGY3d  code.   Not to  pollute  the  global
;;; namespace  we put bgy3d-*  functions into  this module.
;;;
(define-module (guile bgy3d)
  #:use-module (guile molecule)         ; site representation
  #:use-module (guile punch-file)       ; write-punch-file
  #:use-module (srfi srfi-1)            ; list manipulation
  #:use-module (srfi srfi-2)            ; and-let*
  #:use-module (srfi srfi-11)           ; let-values
  #:use-module (ice-9 match)            ; match-lambda
  #:use-module (ice-9 pretty-print)     ; pretty-print
  #:use-module (ice-9 getopt-long)      ; getopt-long, option-ref
  #:use-module (guile bgy3d internal)   ; see bgy3d-guile.c
  #:re-export                           ; from (guile bgy3d internal)
  (state-make
   state-destroy
   vec-make
   vec-make-complex
   vec-destroy
   vec-save
   vec-set-random
   vec-dot
   vec-fft
   vec-ifft
   vec-fft-interp
   vec-map1
   vec-map2
   vec-moments
   vec-shift!
   genpts
   f64dst
   f64+
   f64*
   f64scale
   rism-solvent
   rism-solute
   hnc3d-run-solvent
   hnc3d-run-solute
   bgy3d-run-solvent
   bgy3d-run-solute
   bgy3d-restart-destroy)
  #:export
  (new-main
   old-main
   vec-print
   vec-norm
   bgy3d-test
   solvent/solvent
   solute/solvent))

;;;
;;; The list of the procedures defined in bgy3d-guile.c includes
;;;
;;;   hnc3d-run-solvent
;;;   hnc3d-run-solute
;;;   bgy3d-run-solvent
;;;   bgy3d-run-solute
;;;   bgy3d-pot-interp
;;;   bgy3d-pot-destroy
;;;   bgy3d-restart-destroy
;;;   vec-make
;;;   vec-make-complex
;;;   vec-destroy
;;;   vec-save
;;;   vec-load
;;;   vec-length
;;;   vec-ref
;;;   vec-set-random
;;;   vec-dot
;;;   vec-fft
;;;   vec-ifft
;;;   vec-fft-interp
;;;   vec-map1
;;;   vec-map2
;;;   vec-moments
;;;   vec-shift!
;;;   bgy3d-rank
;;;   bgy3d-size
;;;   state-make
;;;   state-destroy
;;;
;;; and  possibly  more, depending  on  the  compilation options.  The
;;; definitions are  put into (guile bgy3d  internal) module, imported
;;; here  and  re-exported  selectively  (see use-modules  form).   An
;;; earlier  version actively  invoked the  procedure to  define those
;;; symbols. However, this caused annoying warnings with Guile 2 which
;;; complains about "possibly undefined symbols" at compile time.
;;;

;;;
;;; So far C-pointers  are represented by integers in  Scheme. In case
;;; this is ever going to change use this as the null pointer literal:
;;;
(define NULL 0)

;;;
;;; This  may become  define-syntax some  time, so  do not  assume the
;;; arguments are evaluated:
;;;
(define (maybe-print . args)
  (if (zero? (bgy3d-rank))
      (begin
        (apply pretty-print args)
        (force-output))))

;;;
;;; Settings are handled as an  association list. The human input as a
;;; sequence of 2-lists.  This is the  input as it would appear in the
;;; file.  These defaults are targeting water at normal conditions.
;;;
;;; OUTDATED?  To get  butanoic  acid and  hexane  (two largest  ones)
;;; converge when treated by QM one needs about 1500 iterations (well,
;;; we should make the solver better).
;;;
(define *defaults*
  '((solvent "water")                   ; solvent name
    (L 10.0)                            ; [-L, L] gives the box size
    (N 64)                              ; grid dimension
    (rho 0.033427745)                   ; solvent density
    (beta 1.6889)                       ; inverse temperature
    (norm-tol 1.0e-7)                   ; convergence threshold
    (max-iter 1500)                     ; max number of iterations
    (damp-start 1.0)                    ; scaling factor?
    (lambda 0.02)))                     ; not the scheme lambda


;;;
;;; This compacts  association list  with entries coming  later having
;;; precedence. FIXME: thus any prior acons may be ignored.
;;;
(define (overwrite alist)
  (let loop ((merged '())
             (alist alist))
    (if (null? alist)
        merged
        (let ((head (car alist))
              (tail (cdr alist)))
          (loop (assoc-set! merged (car head) (cdr head))
                tail)))))


;;;
;;; The  input comes  from reading  (in  the Scheme  sense) the  human
;;; editable file.   To keep it simple  for a human we  do not require
;;; him/her  to properly  specify pairs  by Scheme  syntax as  in (key
;;; .  value).  Instead  a  key/value shall  be  a 2-list  as in  (key
;;; value), e.g.:
;;;
;;;   ((solute "butanoic acid")
;;;    (closure KH))
;;;
;;; This function is supposed to process the input and return a proper
;;; association list of (key . value) pairs:
;;;
(define (input->settings input)
  (overwrite (map (lambda (x) (cons (first x) (second x)))
                  input)))


;;;
;;; FIXME: must die, use *defaults* directly:
;;;
(define *default-molecule* (assoc-ref (input->settings *defaults*) 'solvent))

;;;
;;; Only uses the length of g1 list to generate file names:
;;;
(define (g1-file-names g1)
  (map (lambda (i) (format #f "g~A.bin" i))
       (iota (length g1))))

;;;
;;; Solute parameters are  currently represented by a name  and a list
;;; of sites. See (guile molecule) for accessors/constructors:
;;;
;;; ("water"
;;;  (("O" (-0.2929 0.0 0.0) 3.1506 0.1521 -0.834)
;;;   ("OH" (0.2929 0.757 0.0) 0.4 0.046 0.417)
;;;   ("OH" (0.2929 -0.757 0.0) 0.4 0.046 0.417)))
;;;

;;;
;;; Update site force field parameters by replacing them with the data
;;; from a  table identified  by string table-name.   The name  of the
;;; site serves as a key to lookup table entries.  The table is stored
;;; in  guile/solutes.scm as  a fake  solute with  all  possible sites
;;; listed once.
;;;
(define (update-sites table-name sites)
  (let* ((table                     ; fake solute with site-parameters
          (molecule-sites (find-molecule table-name)))
         (update-one
          (lambda (site)
            (let* ((name (site-name site))
                   (table-site (assoc name table)))
              (make-site name                      ; original
                         (site-position site)      ; original
                         (site-sigma table-site)   ; from table
                         (site-epsilon table-site) ; from table
                         (site-charge site))))))   ; original
    (map update-one sites)))


;;;
;;; L2-norm   of  a   vec.  The   funciton  vec-dot   is   defined  in
;;; bgy3d-guile.c:
;;;
(define (vec-norm v)
  (sqrt (vec-dot v v)))

;;;
;;; Will work inefficiently with  distributed vectors.  The problem is
;;; vec-ref  is collective.  The  Petsc primitive  VecGetValues() only
;;; works for local section.
;;;
(define (vec-fold kons knil vec)
  "A (left) fold of a (Petsc) vector.  E.g. (vec-fold + 0.0 vec)
computes the sum of all vector elements."
  (let ((len (vec-length vec))
        (vec (lambda (i) (vec-ref vec i))))
    (let loop ((knil knil)
               (i 0))
      (if (< i len)
          (loop (kons (vec i) knil)
                (+ 1 i))
          knil))))

(define (vec-for-each f vec)
  (vec-fold
   (lambda (x seed) (f x))              ; seed not used
   (if #f #f)                           ; seed is #<unspecified>
   vec))

;;;
;;; Returns an "iterator" --- here  a function that takes a callback f
;;; and invokes that callback for each vector element:
;;;
(define (make-vec-iter vec)
  (lambda (f) (vec-for-each f vec)))

;;
;; FIXME: Guile has problems with denormal numbers, replace small ones
;; by zeros. Otherwise Guile will show #.# when printing it:
;;
(define (normalize x)
  (if (> (abs x) 2.0e-308)
      x
      0.0))

(define (vec-print v)
  (vec-for-each
   (lambda (x) (format #t "~a\n" (normalize x)))
   v))

;;;
;;; For   each  position  x   in  the   list  xs   evaluate  potential
;;; v(x). Potential is a BGY3D interator object, not a function:
;;;
(define (potential-map v xs)
  (map (lambda (x) (bgy3d-pot-interp v x)) xs)) ; so far this accepts only one position

;;;
;;; These hooks,  solvent/solvent and solute/solvent,  are called from
;;; the "runqm" script used for QM solutes:
;;;
(define (solvent/solvent input)
  (let* ((settings      (input->settings (append *defaults* input)))
         (solvent-name  (assoc-ref settings 'solvent)) ; string or #f
         (solvent       (find-molecule solvent-name))  ; fails for #f
         (closure       (assoc-ref settings 'closure)) ; symbol or #f
         (run-solvent   (if closure
                            hnc3d-run-solvent
                            bgy3d-run-solvent))
         (solvent-info-file (if closure
                                "x00-fft.bin"
                                "g00.bin"))
         (always-run-solvent #t))
    ;;
    ;; At the moment the function bound to run-solvent echoes settings
    ;; as is, the output is written to disk instead. FIXME: a hack to
    ;; save some time re-running pure solvent calculations.  Set
    ;; always-run-solvent to #f. Then the calulation will only run if
    ;; the file solvent-info-file does not exist. It may still
    ;; correspond to a different solvent or settings.
    ;;
    (if (or always-run-solvent
            (not (file-exists? solvent-info-file)))
        (run-solvent solvent settings)))) ; writes solvent info to disk

;;;
;;; Input  comes   from  the  text  file   (see  ../test-qm/*.scm  and
;;; ../test-qm/runqm) but  sites, funptr and restart come  from the QM
;;; code.  The force  field parameters of sites are  to be replaced by
;;; something meaningful as the QM code knows little about that.
;;;
(define (solute/solvent input sites funptr restart)
  (let* ((settings      (input->settings (append *defaults* input)))
         (solvent-name  (assoc-ref settings 'solvent)) ; string or #f
         (solute-name   (assoc-ref settings 'solute))  ; string or #f
         (solvent       (find-molecule solvent-name))
         (solute        (make-molecule solute-name
                                       (update-sites solute-name
                                                     sites)))
         (closure       (assoc-ref settings 'closure)) ; symbol or #f
         (run-solute    (if closure
                            hnc3d-run-solute
                            bgy3d-run-solute)))
    ;;
    ;; Extend settings by an  entry with the funciton pointer that can
    ;; be used to compute additional solute charge density:
    ;;
    (set! settings (acons 'qm-density   ; key
                          funptr        ; value
                          settings))    ; alist
    ;; Print on master only:
    (maybe-print (list 'SETTINGS: settings))
    (maybe-print (list 'SOLVENT: solvent))
    (maybe-print (list 'SOLUTE: solute))
    (maybe-print (list 'RESTART: restart))
    ;;
    ;; The function  bound to run-solute allocates and  returns a list
    ;; of  Petsc Vecs, a  descriptor for electrostatic  potential, and
    ;; some  restart info  to be passed  back on the  next invokations
    ;; (returned as multiple values). It is the callers responsibility
    ;; to destroy  all of them ---  do not assume any of  them will be
    ;; garbage-collected. As of now it expects the solvent description
    ;; such as g??.bin or x??-fft.bin to exist as files:
    ;;
    (let-values (((g1 potential restart)
                  (run-solute solute solvent settings restart)))
      ;;
      ;; Save g1-files to disk:
      ;;
      (map vec-save (g1-file-names g1) g1)
      ;;
      ;; Dont forget to destroy them after use:
      ;;
      (map vec-destroy g1)
      ;;
      ;; Return (a) the (iterator over) potential (the caller must
      ;; bgy3d-pot-destroy it) and (b) the restart info that, when
      ;; passed to this function next time, could help saving a few
      ;; BGY iterations:
      ;;
      (cons potential restart))))

;;;
;;; Specifications of command line  flags common for old- and new-main
;;; for  use with getopt-long.  Most of  these happen  to be  real- or
;;; integer numbers:
;;;
(define option-spec-base
  (quasiquote
   ((verbosity          (value #t)      (predicate ,string->number)) ; dont use -v
    (N                  (value #t)      (predicate ,string->number))
    (rho                (value #t)      (predicate ,string->number))
    (beta               (value #t)      (predicate ,string->number))
    (norm-tol           (value #t)      (predicate ,string->number))
    (max-iter           (value #t)      (predicate ,string->number))
    (L                  (value #t)      (predicate ,string->number))
    (damp-start         (value #t)      (predicate ,string->number))
    (lambda             (value #t)      (predicate ,string->number))
    (solvent            (value #t)      (predicate ,find-molecule)) ; a string
    (solute             (value #t)      (predicate ,find-molecule)) ; a string
    (dielectric         (value #t)      (predicate ,string->number))
    (from-radial-g2     (value #f))
    (save-guess         (value #f))
    (load-guess         (value #f))
    (snes-solver        (value #t)
                        (predicate ,(lambda (x)
                                    (member x '("jager" "newton" "picard" "trial")))))
    (verbose            (single-char #\v)
                        (value #f)) ; use --verbosity num instead
    (no-cage            (value #f)) ; turns off metallic cage boundary
    (no-hacks           (value #f)))))  ; turns off ugly hacks



(define option-spec-new
  (quasiquote
   ((save-binary        (value #f))
    (unquote-splicing option-spec-base)))) ; common options

;;;
;;; FIXME: ./runbgy.scm invokes old-main at the moment, so that it has
;;; to accept  the options of  the new-main. This leads  to unhelpfull
;;; error messages:
;;;
(define option-spec-all
  (quasiquote
   ((bgy        (value #f))
    (hnc        (value #f))
    (rism       (value #f))
    (closure    (value #t)
                (predicate ,(lambda (x)
                              (member x '("HNC" "PY" "KH")))))
    (unquote-splicing option-spec-base)))) ; common options

;;;
;;; Returns  new   settings  with   updated  fields  taken   from  the
;;; getopt-long options:
;;;
(define (update-settings settings options)
  (let ((update-pair                    ; a function ...
         (match-lambda
          ((key . value)                ; that takes a pair ...
           (let ((op (option-ref options key "")))
             (cons key                     ; and returns updated pair.
                   (or (string->number op) ; converts "" -> #f
                       value)))))))
    (map update-pair settings)))


;;;
;;; For testing, primarily, evaluate  potential at positions of solute
;;; sites and the corresponding total energy:
;;;
(define (maybe-print-potentials solute potential)
  (let* ((sites         (molecule-sites solute))
         (positions     (map site-position sites))
         (potentials    (potential-map potential positions))
         (charges       (map site-charge sites))
         (energy        (apply + (map * charges potentials))))
    ;;
    ;; Only print on master:
    ;;
    (maybe-print (cons 'core-potentials potentials))
    (maybe-print (cons 'core-energy energy))))


;;;
;;; FIXME: at  the moment this  function only emulates the  minimum of
;;; the   functionality   of  the   original   executable.   The   new
;;; functionality is in  flux.  Note that at variance  with the legacy
;;; code   the  function  find-molecule   uses  on-disk   database  in
;;; ./solutes.scm and not the compiled in set from bgy3d-solutes.c and
;;; bgy3d-solvents.h:
;;;
(define (old-main argv)
  (let* ((opts (getopt-long argv option-spec-all))
         (solvent (and-let* ((name (option-ref opts 'solvent #f)))
                            (find-molecule name))) ; Maybe solvent
         (solute (and-let* ((name (option-ref opts 'solute #f)))
                           (find-molecule name)))) ; Maybe solute
    (maybe-print (list 'solvent: solvent))
    (maybe-print (list 'solute: solute))
    (let-values
        (((method run-solvent run-solute) ; three method-dependent values:
          (cond
           ((option-ref opts 'hnc #f)  (values 'hnc hnc3d-run-solvent hnc3d-run-solute))
           ((option-ref opts 'bgy #f)  (values 'bgy bgy3d-run-solvent bgy3d-run-solute))
           ((option-ref opts 'rism #f) (values 'rism rism-solvent rism-solute)))))
      (case method
        ;;
        ;; 3d HNC/BGY.  The functions bound to run-solvent and
        ;; run-solute share the interface. The code that calls them is
        ;; the same:
        ;;
        ((hnc bgy)
         (if solute                     ; either #f or real solute
             ;;
             ;; Solute with solvent. Supply NULL as the restart
             ;; parameter --- we cannot offer anything here:
             ;;
             (let-values (((g1 potential restart)
                           (run-solute solute solvent '() NULL)))
               ;;
               ;; Evaluate and print potential at positions of solute
               ;; sites and the corresponding total energy:
               ;;
               (maybe-print-potentials solute potential)
               ;;
               ;; Write g?.bin files:
               ;;
               (map vec-save (g1-file-names g1) g1)
               ;;
               ;; Then destroy the objects returned:
               ;;
               (bgy3d-restart-destroy restart) ; ignores NULL
               (bgy3d-pot-destroy potential)
               (map vec-destroy g1))
             ;;
             ;; Pure solvent:
             ;;
             (run-solvent solvent '())))
        ;;
        ;; 1d-RISM, very experimental:
        ;;
        ((rism)
         (if solute                     ; either #f or real solute
             ;;
             ;; Solute with solvent:
             ;;
             (run-solute solute solvent '())
             ;;
             ;; Pure solvent:
             ;;
             (run-solvent solvent '()))) ; use Petsc env
        ;;
        ;; Fall through to the new variant:
        ;;
        (else
         (new-main argv))))))


;;;
;;; Linear mesh.  Note that the last point  is max - dr  < max. FIXME:
;;; division by zero for n = 0.
;;;
(define (mesh/lin min max n)
  (let ((dr (/ (- max min) n)))
    (map (lambda (i) (+ min (* dr i)))
         (iota n))))

;;;
;;; Logarithmic mesh.  Similar gotchas  as with mesh/lin. Moreover one
;;; needs min > 0:
;;;
(define (mesh/log min max n)
  (map exp (mesh/lin (log min) (log max) n)))


;;;
;;; Should work for any number of columns:
;;;
(define (print-columns . columns)
  (apply for-each
         (lambda row
           (display (string-join (map number->string row) " "))
           (newline))
         columns))


;;;
;;; Act  according  to the  subcommand  in  (car  argv). With  cmd  ==
;;; "solutes" interprete each argument as the name of the solute. Note
;;; that you may  need to first run a solvent  calculation with cmd ==
;;; "solvent":
;;;
(define (new-main argv)
  (let ((cmd            (car argv))     ; should be non-empty
        (options        (getopt-long argv option-spec-new)))
    (let ((args         (option-ref options '() '())) ; positional arguments (solute names)
          (solvent      (find-molecule (option-ref options 'solvent *default-molecule*)))
          (solute       (find-molecule (option-ref options 'solute *default-molecule*)))
          (save-binary  (option-ref options 'save-binary #f))
          (settings     (update-settings (input->settings *defaults*) options))) ; defaults updated from command line
      (match cmd
        ("solvent"
         ;;
         ;; Only then run pure solvent, if --solvent was present in the
         ;; command line:
         ;;
         (bgy3d-run-solvent solvent settings))
        ;;
        ((or "solute" "solutes")
         ;;
         ;; Check  if we can find  the solutes by names  early, typos are
         ;; common:
         ;;
         (let ((solutes (map find-molecule args)))
           (map (lambda (solute)
                  (let-values (((g1 ve restart)
                                (bgy3d-run-solute solute solvent settings NULL)))
                    ;;
                    ;; Save distributions if requested from command
                    ;; line. FIXME: the file names do not relate to
                    ;; solute, so that when processing more than one
                    ;; solute in a row files will get overwritten:
                    ;;
                    (if save-binary
                        (map vec-save (g1-file-names g1) g1))
                    ;;
                    ;; Dont forget to destroy them:
                    ;;
                    (bgy3d-pot-destroy ve) ; not yet used
                    (bgy3d-restart-destroy restart)
                    (map vec-destroy g1)))
                solutes)))
        ("punch"
         ;;
         ;; Use g1 vectors to produce a *.pun file for visualization:
         ;;
         (let ((g1 (map vec-load args)))
            (write-punch-file solute
                              (map vec-length g1)
                              (map make-vec-iter g1)
                              settings)
            (map vec-destroy g1)))
        ;;
        ("rdf"
         ;;
         ;; Print Gnuplot-ready RDF table. E.g. by executing this:
         ;;
         ;; mpirun guile/runbgy.scm rdf --N 96 --L 10 0 0 0 test/*.bin
         ;;
         (let* ((center (map string->number (take args 3)))
                (args (drop args 3))
                (angular-order 110)     ; FIXME: literals here ...
                (rmin 0.75)             ; cannot be zero for log-mesh
                (rmax (assoc-ref settings 'L)) ; a choice ...
                (npts (assoc-ref settings 'N)) ; another choice ...
                (mesh (mesh/log rmin rmax npts))
                (domain (state-make settings))
                (rdfs (map (lambda (path)
                             (let* ((g (vec-load path))
                                    (rdf (rism-rdf domain g center mesh angular-order)))
                               (vec-destroy g)
                               rdf))
                           args)))
           (state-destroy domain)
           (if (zero? (bgy3d-rank))
               (begin
                 (format #t "# center = ~A\n" center)
                 (apply print-columns mesh rdfs)))))
        ;;
        ("dump"
         ;;
         ;; Dump each Vec from a *.bin file to tty:
         ;;
         (for-each
          (lambda (path)
            (let ((v (vec-load path)))
              (vec-print v)
              (vec-destroy v)))
          args))
        ;;
        ("print-molecule"
         (for-each
             (lambda (name)
               (print-molecule/xyz (find-molecule name)))
           args))
        ;;
        ("find-molecule"
         (for-each
             (lambda (name)
               (let* ((mol (find-molecule name))
                      (charge-sum (molecule-charge mol))
                      (dipole-vec (molecule-dipole mol))
                      (dipole-mag (sqrt (apply + (map * dipole-vec dipole-vec)))))
                 ;;
                 ;; Abuse debugging helper to print commented expressions:
                 ;;
                 (pk 'CHARGE-SUM: charge-sum)
                 (pk 'DIPOLE-VEC: dipole-vec)
                 (pk 'DIPOLE-MAG: dipole-mag)
                 ;;
                 ;; Actual data:
                 ;;
                 (pretty-print mol)))
           args))))))


(define (bgy3d-test)
  (pretty-print "does nothing"))
;;;
;;; Scheme  interface  to  BGY3d  code.   Not to  pollute  the  global
;;; namespace  we put bgy3d-*  functions into  this module.
;;;
(define-module (guile bgy3d)
  #:use-module (guile compat)           ; define-syntax-rule for 1.8
  #:use-module (guile molecule)         ; site representation
  #:use-module (guile punch-file)       ; write-punch-file
  #:use-module (guile utils)            ; memoize, ddd
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
   vec-scale!
   genpts
   f64dst
   f64+
   f64*
   f64scale
   rism-solvent
   rism-solute
   rism-self-energy
   hnc3d-run-solvent
   hnc3d-run-solute
   bgy3d-run-solvent
   bgy3d-run-solute
   bgy3d-restart-destroy)
  #:export
  (bgy3d-main
   vec-print
   vec-norm
   solvent/solvent
   solute/solvent))

;;;
;;; The list of the procedures defined in bgy3d-guile.c includes
;;;
;;;   rism-solvent
;;;   rism-solute
;;;   rism-self-energy
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
;;;   vec-scale!
;;;   comm-rank
;;;   comm-size
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
;;; This  form can  be  used for  code  to be  evaluated  only on  one
;;; (master)  worker. In  particular no  collective procedures  can be
;;; used from this context. Do not assume the forms are evaluated:
;;;
(define-syntax-rule (begin/serial e ...)
  (begin
    (if (zero? (comm-rank))
        (begin e ...))
    (if #f #f)))

(define (pretty-print/serial . args)
  (begin/serial
   (apply pretty-print args)
   (force-output)))

;;;
;;; Settings are handled as an  association list. The human input as a
;;; sequence of 2-lists.  This is the  input as it would appear in the
;;; file.  These defaults are targeting water at normal conditions.
;;;
;;; The  defaults, updated  from the  command  line, will  be used  to
;;; prepare   the   reall   association   list  with   settings.   See
;;; parse-command-line below.
;;;
(define *defaults*
  '((solvent #f)                      ; uninitialized string, actually
    (solute #f)                       ; uninitialized string
    (L 10.0)                          ; [-L, L] gives the box size
    (N 64)                            ; grid dimension
    (rho 0.033427745)                 ; solvent density (1 g/cm^3)
    (beta 1.6889)                     ; inverse temperature
    (norm-tol 1.0e-7)                 ; convergence threshold
    (max-iter 1500)                   ; max number of iterations
    (damp-start 1.0)                  ; scaling factor?
    (lambda 0.02)                     ; not the scheme lambda
    (bond-length-thresh 2.0) ; lower bound for non-bonding interactions, in A.
    (closure HNC)                     ; HNC, KH or PY
    (hnc #f)                          ; FIXME: exclusive
    (bgy #f)                          ; FIXME: exclusive
    (rism #f)                         ; FIXME: exclusive
    (derivatives #f)                  ; #t or #f
    ))


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
    (pretty-print/serial (list 'SETTINGS: settings))
    (pretty-print/serial (list 'SOLVENT: solvent))
    (pretty-print/serial (list 'SOLUTE: solute))
    (pretty-print/serial (list 'RESTART: restart))
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
;;; FIXME: ./runbgy.scm invokes old-main at the moment, so that it has
;;; to accept  the options of  the new-main. This leads  to unhelpfull
;;; error messages sometimes:
;;;
(define option-spec
  (quasiquote
   ((bgy                (value #f))
    (hnc                (value #f))
    (rism               (value #f))
    (verbosity          (value #t)      (predicate ,string->number)) ; dont use -v
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
    (closure            (value #t)      (predicate ,(lambda (x)
                                                      (member x '("HNC" "PY" "KH")))))
    (dielectric         (value #t)      (predicate ,string->number))
    (comb-rule          (value #t)      (predicate ,string->number))
    (solvent-1d         (value #f)) ; use solvent susceptibility from 1D RISM
    (from-radial-g2     (value #f))
    (save-guess         (value #f))
    (load-guess         (value #f))
    (derivatives        (value #f))
    (snes-solver        (value #t)
                        (predicate ,(lambda (x)
                                    (member x '("jager" "newton" "picard" "trial")))))
    (verbose            (single-char #\v)
                        (value #f)) ; use --verbosity num instead
    (rbc                (value #f)) ; add repulsive bridge correction
    (no-cage            (value #f)) ; turns off metallic cage boundary
    (no-hacks           (value #f)))))  ; turns off ugly hacks


;;;
;;; Returns  new   settings  with   updated  fields  taken   from  the
;;; getopt-long options. The values from the command line (provided as
;;; strings) are converted to the proper type derived from the default
;;; value of the setting.
;;;
(define (update-settings settings options)
  ;;
  ;; Here "old" is the default and should have the proper type (e.g. a
  ;; number, where numbers are expected), "new" is either a string or
  ;; boolean (for pure or missing options). FIXME: there is an ugly
  ;; exception though --- we use old == #f to indicate no meanigful
  ;; default for solute and solvent strings.
  ;;
  (define (update old new)
    (cond
     ((and (boolean? old) (boolean? new)) new) ; for pure options
     ((not new) old) ; keep the default if not specified in command line
     ((string? old) new)
     ((number? old) (string->number new))
     ((symbol? old) (string->symbol new))
     ((not old) new))) ; FIXME: if the default is #f use new, whatever that is
  ;;
  ;; ... otherwise someting went wrong and the return value of (update
  ;; old new) will be unspecified.
  ;;
  (define update-pair                   ; a function ...
    (match-lambda
     ((key . val)                              ; that takes a pair ...
      (let ((str (option-ref options key #f))) ; string or #f
        (cons key (update val str)))))) ; and returns updated pair.
  ;;
  ;; FIXME: this updates the  defaults, but does not add other entries
  ;; specified in the command line. Augment the defaults.
  ;;
  (map update-pair settings))


;;;
;;; Try  to  keep the  command  line parsing  at  a  single place  for
;;; re-use. This converts at least some known options to numbers. Note
;;; that  as  implemented,  the  function  update-settings,  does  not
;;; introduce new key/value pairs  not present in defaults. FIXME: The
;;; flags  present  in  option-spec  but  not  in  defaults  are  this
;;; inaccessible (unless one uses PETSC command line accessors).
;;;
(define (parse-command-line argv)
  (let* ((options (getopt-long argv option-spec))
         (settings (update-settings (input->settings *defaults*) options)))
    (acons '() (option-ref options '() '()) settings))) ; positional args


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
    (pretty-print/serial (cons 'core-potentials potentials))
    (pretty-print/serial (cons 'core-energy energy))))


;;;
;;; FIXME: at  the moment this  function only emulates the  minimum of
;;; the   functionality   of  the   original   executable.   The   new
;;; functionality is in  flux.  Note that at variance  with the legacy
;;; code   the  function  find-molecule   uses  on-disk   database  in
;;; ./solutes.scm and not the compiled in set from bgy3d-solutes.c and
;;; bgy3d-solvents.h:
;;;
(define (old-main argv)
  (let* ((settings (parse-command-line argv))
         (solvent (and-let* ((name (assoc-ref settings 'solvent)))
                    (find-molecule name))) ; Maybe solvent
         (solute (and-let* ((name (assoc-ref settings 'solute)))
                   (find-molecule name)))) ; Maybe solute
    (pretty-print/serial (list 'solvent: solvent))
    (pretty-print/serial (list 'solute: solute))
    (let-values
        (((method run-solvent run-solute) ; three method-dependent values:
          (cond
           ((assoc-ref settings 'hnc)  (values 'hnc hnc3d-run-solvent hnc3d-run-solute))
           ((assoc-ref settings 'bgy)  (values 'bgy bgy3d-run-solvent bgy3d-run-solute))
           ((assoc-ref settings 'rism) (values 'rism rism-solvent rism-solute)))))
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
             (let ((res (run-solute solute solvent '())))
               (pretty-print/serial res))
             ;;
             ;; Pure solvent:
             ;;
             (let ((res (run-solvent solvent '())))
               (pretty-print/serial res))))
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
;;; Construct a PES as a  function of solute (or solvent) geometry and
;;; return that together with the initial geometry as multiple values.
;;; The solvent  susceptibility χ is precomuted lazily  and forced the
;;; first time the PES is used (if at all).
;;;
(define (make-pes solute solvent settings)
  (let ((chi (delay (let ((dct (rism-solvent solvent settings))) ; solvent run here!
                      (assoc-ref dct 'susceptibility)))) ; extract susceptibility
        (section (if solute 'solute 'solvent)))          ; alist key
    (let* ((m (or solute solvent))      ; molecule to be "moved"
           (x0 (molecule-positions m))  ; unperturbed geometry
           (rism (memoize               ; FIXME: leaking much?
                  (lambda (x s)         ; (x, settings) -> alist
                    (let* ((m' (move-molecule m x))
                           (d (if solute
                                  (rism-solute m' solvent s (force chi))
                                  (rism-solvent m' s))))
                      (assoc-ref d section))))) ; u or v
           ;; PES (f x) -> energy:
           (f (lambda (x)
                (assoc-ref (rism x settings) 'XXX)))
           ;; Make sure to request evaluation of gradients for (g x)
           ;; but not for (f x). Using two different settings impedes
           ;; memoization:
           (settings (if (assoc-ref settings 'derivatives)
                         settings       ; for the sake of memoization
                         (acons 'derivatives #t settings)))
           ;; PES tuple (fg x) -> energy, gradients:
           (fg (lambda (x)
                 (let ((dict (rism x settings)))
                   (values (assoc-ref dict 'XXX)
                           (assoc-ref dict 'XXX-GRADIENTS))))))
      ;;
      ;; Return initial geometry, PES function f, and its Taylor
      ;; expansion fg as multiple values:
      ;;
      (values x0 f fg))))

;;;
;;; "Gas-phase" PES. Hm,  you cannot search for a  minimum on this PES
;;; without  introducing intra-molecular  interactions or  making each
;;; species rigid.
;;;
(define (make-pes/gp solute settings)
  (let* ((rmax (assoc-ref settings 'bond-length-thresh))
         (species (molecule-species solute rmax))
         (x0 (molecule-positions solute))
         (f (lambda (x)                 ; x -> energy
              (let ((solute' (move-molecule solute x)))
                (rism-self-energy solute' species))))
         (g (lambda (x)                 ; x -> gradients
              (ddd f x)))
         (fg (lambda (x)
               (values (f x) (g x)))))
    (values x0 f fg)))


;;;
;;; Act   according   to   the   subcommand  (the   first   positional
;;; argument). With  cmd == "solutes" interprete each  argument as the
;;; name of the solute. Note that  you may need to first run a solvent
;;; calculation with cmd == "solvent":
;;;
(define (new-main argv)
  (let* ((settings (parse-command-line argv)) ; argv[0] is ignored
         (args (assoc-ref settings '())) ; positional arguments
         (cmd (car args))                ; first the command ...
         (args (cdr args)))              ; ... then the real args
    (let ((solvent (and-let* ((name (assoc-ref settings 'solvent)))
                     (find-molecule name))) ; Maybe solvent
          (solute (and-let* ((name (assoc-ref settings 'solute)))
                    (find-molecule name))) ; Maybe solute
          (save-binary (assoc-ref settings 'save-binary)))
      ;;
      (match cmd
        ;;
        ((or "energy" "gradients" "taylor")
         (let-values (((x f fg) (if solvent
                                    (make-pes solute solvent settings)
                                    (make-pes/gp solute settings))))
           (match cmd
             ("energy"
              (let ((e (f x)))
                (pretty-print/serial e)))
             ("gradients"
              (let-values (((e g) (fg x)))
                (pretty-print/serial g)))
             ("taylor"
              (let-values (((e g) (fg x)))
                (pretty-print/serial (cons e g)))))))
        ;;
        ;; Start a server that communicates with a client via two
        ;; named pipes for input/output:
        ;;
        ("server"
         (let ((finp (first args))  ; coordinates or #f read from here
               (fout (second args))) ; resulting energy is written here
           ;;
           ;; Make a function of 3d geometry. The initial geometry is
           ;; not used:
           ;;
           (let-values (((x0 f fg) (if solvent
                                       (make-pes solute solvent settings)
                                       (make-pes/gp solute settings))))
             (let loop ()
               ;;
               ;; Read the input pipe. This will block utill the
               ;; client writes the geometry from the other side. Note
               ;; that only the first s-expression is read:
               ;;
               (let ((x (with-input-from-file finp read)))
                 ;;
                 ;; Convention is when the input is #f, then
                 ;; terminate. Otherwise evaluate PES at this point
                 ;; and write the resulting energy and gradients to
                 ;; our side of the output pipe as a single
                 ;; s-expression. This will block utill the client
                 ;; reads the results from his end:
                 ;;
                 (when x
                   (let-values (((e g) (fg x))) ; get energy, gradient
                     (with-output-to-file fout
                       (lambda () (write (cons e g))))) ; or (list e g)?
                   (loop)))))))
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
        ("moments"
         ;;
         ;; Print a few first moments (with respect to the grid
         ;; center) of each distribution supplied in the command line.
         ;;
         ;; mpirun guile/runbgy.scm moments --N 96 --L 10 test/*.bin
         ;;
         (let ((domain (state-make settings)))
           (for-each (lambda (path)
                       (let* ((h (vec-load path))
                              (h (vec-shift! h -1.0))
                              (h (vec-scale! h -1.0))
                              (moments (vec-moments domain h)))
                         (begin/serial
                           (format #t ";;; Moments from ~A\n" path)
                           (pretty-print moments))
                         (vec-destroy h)))
                     args)
           (state-destroy domain)))
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
           (begin/serial
            (format #t "# center = ~A\n" center)
            (apply print-columns mesh rdfs))))
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
        ("print-species"
         (let ((rmax (assoc-ref settings 'bond-length-thresh)))
          (for-each
           (lambda (name)
             (let ((mol (find-molecule name)))
               (pretty-print/serial
                (map cons
                     (map site-name (molecule-sites mol))
                     (molecule-species mol rmax)))))
           args)))
        ;;
        ("self-energy"
         (for-each
          (lambda (name)
            (let-values (((x0 f fg) (make-pes/gp (find-molecule name) settings)))
              (let ((e (f x0)))
                (pretty-print/serial e))))
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


;;;
;;; We are  trying to  emulate behaviour of  old executable  unless we
;;; find a better  interface. The old behaviour is  triggered only for
;;; --bgy, --hnc and --rism options:
;;;
(define (bgy3d-main argv)
  (if (or (member "--bgy" argv)
          (member "--hnc" argv)
          (member "--rism" argv))
      (old-main argv)
      (new-main argv)))

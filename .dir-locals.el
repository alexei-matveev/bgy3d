;; Per-directory local variables for GNU Emacs 23 and later.
((scheme-mode
  . ((indent-tabs-mode . nil)
     (eval . (put 'and-let* 'scheme-indent-function 1))
     (eval . (put 'match 'scheme-indent-function 1))
     (eval . (put 'unless 'scheme-indent-function 1))
     (eval . (put 'when 'scheme-indent-function 1))
     (eval . (put 'with-fluids 'scheme-indent-function 1)))))

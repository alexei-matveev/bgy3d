module lisp
  !
  ! Copyright (c) 2013 Alexei Matveev
  !
  use iso_c_binding
  implicit none
  private

  type, bind (c) :: obj
     private
     integer (c_intptr_t) :: ptr
  end type obj

  ! FIXME: this is for Guile 2.0:
  type (obj), parameter :: false = obj (z"004")       ! #f
  type (obj), parameter :: nil = obj (z"304")         ! ()
  type (obj), parameter :: true = obj (z"404")        ! #t
  type (obj), parameter :: unspecified = obj (z"804") ! *unspecified*
  type (obj), parameter :: undefined = obj (z"904")   ! for SCM_UNBNDP()
  type (obj), parameter :: eof = obj (z"a04")         ! for eof-object?

  interface values
     function scm_values (list) result (vals) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: list
       type (obj) :: vals
     end function scm_values
  end interface values

  interface car
     function scm_car (pair) result (car) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: pair
       type (obj) :: car
     end function scm_car
  end interface car

  interface cdr
     function scm_cdr (pair) result (cdr) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: pair
       type (obj) :: cdr
     end function scm_cdr
  end interface cdr

  interface cddr
     function scm_cddr (list) result (item) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: list
       type (obj) :: item
     end function scm_cddr
  end interface cddr

  interface cdar
     function scm_cdar (list) result (item) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: list
       type (obj) :: item
     end function scm_cdar
  end interface cdar

  interface cadr
     function scm_cadr (list) result (item) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: list
       type (obj) :: item
     end function scm_cadr
  end interface cadr

  interface caar
     function scm_caar (list) result (item) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: list
       type (obj) :: item
     end function scm_caar
  end interface caar

  interface assoc
     function scm_assoc (key, alist) result (pair) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: key, alist
       type (obj) :: pair
     end function scm_assoc
  end interface assoc

  interface cons
     function scm_cons (car, cdr) result (pair) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: car, cdr
       type (obj) :: pair
     end function scm_cons
  end interface cons

  interface flonum
     function scm_from_double (d) result (inexact) bind (c)
       import
       implicit none
       real (c_double), intent (in), value :: d
       type (obj) :: inexact
     end function scm_from_double

     function scm_to_double (inexact) result (d) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: inexact
       real (c_double) :: d
     end function scm_to_double
  end interface flonum

  interface int
     function scm_from_int32 (i) result (exact) bind (c)
       import
       implicit none
       integer (c_int32_t), intent (in), value :: i
       type (obj) :: exact
     end function scm_from_int32

     function scm_from_int64 (i) result (exact) bind (c)
       import
       implicit none
       integer (c_int64_t), intent (in), value :: i
       type (obj) :: exact
     end function scm_from_int64

     function scm_to_int32 (exact) result (i) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: exact
       integer (c_int32_t) :: i
     end function scm_to_int32
  end interface int

  interface bool
     module procedure from_bool
     module procedure to_bool
  end interface bool

  interface
     function scm_from_locale_symboln (str, len) result (symbol) bind (c)
       import
       implicit none
       character (kind=c_char), intent (in) :: str(*)
       integer (c_size_t), intent (in), value :: len
       type (obj) :: symbol
     end function scm_from_locale_symboln

     function scm_from_locale_stringn (str, len) result (string) bind (c)
       import
       implicit none
       character (kind=c_char), intent (in) :: str(*)
       integer (c_size_t), intent (in), value :: len
       type (obj) :: string
     end function scm_from_locale_stringn

     function scm_to_bool (bool) result (int) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: bool
       integer (c_int) :: int
     end function scm_to_bool

     function scm_display (object, port) result (undef) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: object, port
       type (obj) :: undef
     end function scm_display

     function scm_newline (port) result (undef) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: port
       type (obj) :: undef
     end function scm_newline

     function scm_variable_ref (var) result (val) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: var
       type (obj) :: val
     end function scm_variable_ref

     function scm_c_resolve_module (path) result (mod) bind (c)
       import
       implicit none
       character (kind=c_char), intent (in) :: path
       type (obj) :: mod
     end function scm_c_resolve_module

     function scm_c_module_lookup (mod, key) result (var) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: mod
       character (kind=c_char), intent (in) :: key
       type (obj) :: var
     end function scm_c_module_lookup
  end interface

  interface not
     function scm_not (object) result (bool) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: object
       type (obj) :: bool
     end function scm_not
  end interface

  interface list
     function scm_list_1 (e1) result (list) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: e1
       type (obj) :: list
     end function scm_list_1

     function scm_list_2 (e1, e2) result (list) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: e1, e2
       type (obj) :: list
     end function scm_list_2

     function scm_list_3 (e1, e2, e3) result (list) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: e1, e2, e3
       type (obj) :: list
     end function scm_list_3

     function scm_list_4 (e1, e2, e3, e4) result (list) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: e1, e2, e3, e4
       type (obj) :: list
     end function scm_list_4

     function scm_list_5 (e1, e2, e3, e4, e5) result (list) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: e1, e2, e3, e4, e5
       type (obj) :: list
     end function scm_list_5
 end interface list

 interface funcall
    function scm_call_0 (proc) result (res) bind (c)
      import
      implicit none
      type (obj), intent (in), value :: proc
      type (obj) :: res
    end function scm_call_0

    function scm_call_1 (proc, arg1) result (res) bind (c)
      import
      implicit none
      type (obj), intent (in), value :: proc, arg1
      type (obj) :: res
    end function scm_call_1

    function scm_call_2 (proc, arg1, arg2) result (res) bind (c)
      import
      implicit none
      type (obj), intent (in), value :: proc, arg1, arg2
      type (obj) :: res
    end function scm_call_2

    function scm_call_3 (proc, arg1, arg2, arg3) result (res) bind (c)
      import
      implicit none
      type (obj), intent (in), value :: proc, arg1, arg2, arg3
      type (obj) :: res
    end function scm_call_3

    function scm_call_4 (proc, arg1, arg2, arg3, arg4) result (res) bind (c)
      import
      implicit none
      type (obj), intent (in), value :: proc, arg1, arg2, arg3, arg4
      type (obj) :: res
    end function scm_call_4

    function scm_call_5 (proc, arg1, arg2, arg3, arg4, arg5) result (res) bind (c)
      import
      implicit none
      type (obj), intent (in), value :: proc, arg1, arg2, arg3, arg4, arg5
      type (obj) :: res
    end function scm_call_5

    !
    ! FIXME: extend up to 9 if necessary.
    !
 end interface

  public :: obj                              ! type
  public :: flonum, int, symbol, string, bool ! constructors
  public :: is_true, true, false
  public :: cons, nil
  public :: car, cdr
  public :: cddr, cdar, cadr, caar
  public :: acons, assoc, not
  public :: display, newline
  public :: values, list
  public :: funcall
  public :: lookup

contains

  subroutine display (object, port)
    implicit none
    type (obj), intent (in) :: object, port
    optional :: port
    ! *** end of interface **

    type (obj) :: ignore

    if (present (port)) then
       ignore = scm_display (object, port)
    else
       ignore = scm_display (object, undefined)
    endif
  end subroutine display


  subroutine newline (port)
    implicit none
    type (obj), intent (in), optional :: port
    ! *** end of interface **

    type (obj) :: ignore

    if (present (port)) then
       ignore = scm_newline (port)
    else
       ignore = scm_newline (undefined)
    endif
  end subroutine newline


  function symbol (fstr) result (symb)
    implicit none
    character (len=*), intent (in) :: fstr
    type (obj) :: symb
    ! *** end of interface **

    integer (c_size_t) :: slen

    slen = len (fstr)
    symb = scm_from_locale_symboln (fstr, slen)
  end function symbol


  function string (fstr) result (lstr)
    implicit none
    character (len=*), intent (in) :: fstr
    type (obj) :: lstr
    ! *** end of interface **

    integer (c_size_t) :: slen

    slen = len (fstr)
    lstr = scm_from_locale_stringn (fstr, slen)
  end function string


  function from_bool (flag) result (bool)
    implicit none
    logical, intent (in) :: flag
    type (obj) :: bool
    ! *** end of interface **

    ! SCM scm_from_bool(int) is a C-macro!
    if (flag) then
       bool = true
    else
       bool = false
    endif
  end function from_bool


  function to_bool (bool) result (flag)
    implicit none
    type (obj), intent (in) :: bool
    logical :: flag
    ! *** end of interface **

    ! This one will choke on anything  that is not #t or #f, the usual
    ! convention that any object except  #f is "trueish" does not hold
    ! here:
    flag = (scm_to_bool (bool) /= 0)
  end function to_bool


  function is_true (x) result (y)
    implicit none
    type (obj), intent (in) :: x
    logical :: y
    ! *** end of interface **

    y = .not. bool (not (x))
  end function is_true


  function acons (key, val, old) result (new)
    !
    ! (acons k v lst) => (cons (cons k v) lst)
    !
    implicit none
    type (obj), intent (in) :: key, val, old
    type (obj) :: new
    ! *** end of interface **

    new = cons (cons (key, val), old)
  end function acons


  function lookup (module, name) result (val)
    !
    ! Lookup a name in a module.
    !
    use iso_c_binding, only: C_NULL_CHAR
    implicit none
    character (len=*), intent (in) :: module, name
    type (obj) :: val
    ! *** end of interface ***

    val = scm_variable_ref (scm_c_module_lookup &
         (scm_c_resolve_module (module // C_NULL_CHAR), name // C_NULL_CHAR))
  end function lookup

end module lisp

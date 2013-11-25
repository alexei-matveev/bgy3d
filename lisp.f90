module lisp

  use iso_c_binding
  implicit none
  private

  type, bind (c) :: obj
     private
     integer (c_intptr_t) :: ptr
  end type obj

  ! FIXME: this is for Guile 2.0:
  type (obj), parameter :: nil = obj (z"304")

  interface cons
     function scm_cons (car, cdr) result (pair) bind (c)
       import
       implicit none
       type (obj), intent (in), value :: car, cdr
       type (obj) :: pair
     end function scm_cons
  end interface cons

  interface float
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
  end interface float

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
  end interface

  public :: obj                        ! type
  public :: float, int, symbol, string ! constructors
  public :: cons, nil
  public :: acons

contains

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

end module lisp

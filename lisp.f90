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

  interface num
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

     function scm_from_double (d) result (inexact) bind (c)
       import
       implicit none
       real (c_double), intent (in), value :: d
       type (obj) :: inexact
     end function scm_from_double
  end interface num

  interface
     function scm_from_utf8_symboln (str, len) result (symbol) bind (c)
       import
       implicit none
       character (kind=c_char), intent (in) :: str(*)
       integer (c_size_t), intent (in), value :: len
       type (obj) :: symbol
     end function scm_from_utf8_symboln

     function scm_from_utf8_stringn (str, len) result (string) bind (c)
       import
       implicit none
       character (kind=c_char), intent (in) :: str(*)
       integer (c_size_t), intent (in), value :: len
       type (obj) :: string
     end function scm_from_utf8_stringn
  end interface

  public :: obj                 ! type
  public :: num, sym, str       ! constructors
  public :: cons, nil

contains

  function sym (fstr) result (symb)
    implicit none
    character (len=*), intent (in) :: fstr
    type (obj) :: symb
    ! *** end of interface **

    integer (c_size_t) :: slen

    slen = len (fstr)
    symb = scm_from_utf8_symboln (fstr, slen)
  end function sym

  function str (fstr) result (lstr)
    implicit none
    character (len=*), intent (in) :: fstr
    type (obj) :: lstr
    ! *** end of interface **

    integer (c_size_t) :: slen

    slen = len (fstr)
    lstr = scm_from_utf8_stringn (fstr, slen)
  end function str

end module lisp

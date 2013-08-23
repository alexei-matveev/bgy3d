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

  public :: obj
  public :: cons, nil
  public :: num

end module lisp

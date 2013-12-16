module options
  !
  ! Copyright (c) 2013 Alexei Matveev
  !
  use kinds, only: rk
  implicit none
  private

  interface getopt
     module procedure getenv_test
     module procedure getenv_bool
     module procedure getenv_real
     module procedure getenv_int

     module procedure getopt_test
     module procedure getopt_int
     module procedure getopt_real
     module procedure getopt_string
  end interface getopt

  public :: getopt

contains

  !
  ! These are  used to inquire settings organized  into an association
  ! list supplied as the first argument:
  !
  function getenv_test (env, key) result (ok)
    use lisp, only: obj, assoc, symbol, not, bool
    implicit none
    type (obj), intent (in) :: env
    character (len=*), intent (in) :: key
    logical :: ok
    ! *** end of interface ***

    type (obj) :: pair

    ! "not not" converts any value to true or false:
    pair = assoc (symbol (key), env)
    ok = .not. bool (not (pair))
  end function getenv_test


  function getenv_bool (env, key, val) result (ok)
    use lisp, only: obj, assoc, symbol, not, bool, cdr
    implicit none
    type (obj), intent (in) :: env
    character (len=*), intent (in) :: key
    logical, intent (inout) :: val
    logical :: ok
    ! *** end of interface ***

    type (obj) :: pair

    pair = assoc (symbol (key), env)
    ok = .not. bool (not (pair))
    if (ok) val = bool (cdr (pair))
  end function getenv_bool


  function getenv_real (env, key, val) result (ok)
    use lisp, only: obj, assoc, symbol, not, bool, cdr, flonum
    implicit none
    type (obj), intent (in) :: env
    character (len=*), intent (in) :: key
    real (rk), intent (inout) :: val
    logical :: ok
    ! *** end of interface ***

    type (obj) :: pair

    pair = assoc (symbol (key), env)
    ok = .not. bool (not (pair))
    if (ok) val = flonum (cdr (pair))
  end function getenv_real


  function getenv_int (env, key, val) result (ok)
    use lisp, only: obj, assoc, symbol, not, bool, cdr, int
    implicit none
    type (obj), intent (in) :: env
    character (len=*), intent (in) :: key
    integer, intent (inout) :: val
    logical :: ok
    ! *** end of interface ***

    type (obj) :: pair

    pair = assoc (symbol (key), env)
    ok = .not. bool (not (pair))
    if (ok) val = int (cdr (pair))
  end function getenv_int

  !
  ! These are to inquire PETSC environment:
  !
  function op (key) result (str)
    use iso_c_binding, only: C_NULL_CHAR
    implicit none
    character (len=*), intent (in) :: key
    character (len = len (key) + 3) :: str
    ! *** end of interface ***

    str = "--" // key // C_NULL_CHAR
  end function op

  function getopt_test (key) result (ok)
    use foreign, only: bgy3d_getopt_test
    implicit none
    character (len=*), intent (in) :: key
    logical :: ok
    ! *** end of interface ***

    ok = bgy3d_getopt_test (op (key))
  end function getopt_test

  function getopt_int (key, val) result (ok)
    use foreign, only: bgy3d_getopt_int
    implicit none
    character (len=*), intent (in) :: key
    integer, intent (inout) :: val
    logical :: ok
    ! *** end of interface ***

    ok = bgy3d_getopt_int (op (key), val)
  end function getopt_int

  function getopt_real (key, val) result (ok)
    use foreign, only: bgy3d_getopt_real
    implicit none
    character (len=*), intent (in) :: key
    real (rk), intent (inout) :: val
    logical :: ok
    ! *** end of interface ***

    ok = bgy3d_getopt_real (op (key), val)
  end function getopt_real

  function getopt_string (key, val) result (ok)
    use foreign, only: bgy3d_getopt_string
    implicit none
    character (len=*), intent (in) :: key
    character (len=*), intent (inout) :: val
    logical :: ok
    ! *** end of interface ***

    ok = bgy3d_getopt_string (op (key), len (val), val)
  end function getopt_string

end module options

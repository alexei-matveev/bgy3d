module options
  !
  ! Copyright (c) 2013 Alexei Matveev
  !
  use kinds, only: rk
  implicit none
  private

  interface getopt
     module procedure getopt_test
     module procedure getopt_int
     module procedure getopt_real
     module procedure getopt_string
     module procedure getopt_obj
     module procedure getopt_bool
  end interface getopt

  public :: getopt

contains

  !
  ! These are to  inquire dynamic environment. It is  assumed that the
  ! dynvar   *settings*    in   bgy3d.scm   is    set   to   something
  ! meaningful. E.g. as specified in the command line.
  !
  function op (key) result (str)
    use iso_c_binding, only: C_NULL_CHAR
    implicit none
    character (len=*), intent (in) :: key
    character (len = len (key) + 1) :: str
    ! *** end of interface ***

    str = key // C_NULL_CHAR
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

  function getopt_obj (key, val) result (ok)
    use foreign, only: bgy3d_getopt_scm
    use lisp, only: obj
    implicit none
    character (len=*), intent (in) :: key
    type (obj), intent (inout) :: val
    logical :: ok
    ! *** end of interface ***

    ok = bgy3d_getopt_scm (op (key), val)
  end function getopt_obj

  function getopt_bool (key, val) result (ok)
    use foreign, only: bgy3d_getopt_scm
    use lisp, only: obj, bool
    implicit none
    character (len=*), intent (in) :: key
    logical, intent (inout) :: val
    logical :: ok
    ! *** end of interface ***

    type (obj) :: v

    ok = getopt_obj (key, v)
    if (ok) val = bool (v)
  end function getopt_bool

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

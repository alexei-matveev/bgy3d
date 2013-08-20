module options
  use kinds, only: rk
  implicit none
  private

  interface getopt
     module procedure getopt_test
     module procedure getopt_int
     module procedure getopt_real
     module procedure getopt_string
  end interface getopt

  public :: getopt

contains

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

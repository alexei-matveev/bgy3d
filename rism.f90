program rism
  use iso_c_binding, only: rk => c_double
  implicit none

  interface
     subroutine rism_dst (n, out, in) bind (c)
       !
       ! See rism-dst.c
       !
       use iso_c_binding, only: c_size_t, c_double
       implicit none
       integer(c_size_t), intent(in), value :: n
       real(c_double), intent(out) :: out (n)
       real(c_double), intent(in) :: in (n)
     end subroutine rism_dst
  end interface

  integer :: i

  do i = 1, 10**4
     call test_dst (2**12)
  end do

contains

  function dst (a) result (b)
    use iso_c_binding, only: c_size_t
    implicit none
    real(rk), intent(in) :: a(:)
    real(rk) :: b(size (a))
    ! *** end of interface ***

    integer(c_size_t) :: n

    ! cast to size_t
    n = size (a)

    call rism_dst (n, b, a)
  end function dst

  subroutine test_dst (n)
    implicit none
    integer, intent(in) :: n
    ! *** end of interface ***

    real(RK) :: a(n)

    call random_number (a)

    ! RODFT11 (DST-IV) is self inverse up to this normalization factor:
    if (maxval (a - dst (dst (a)) / (2 * n)) > 1.0e-14) stop "does not match"
  end subroutine test_dst
end program

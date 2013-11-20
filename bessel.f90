module bessel
  !
  ! One of  the many equivalent expresssions for  the spherical Bessel
  ! functions:
  !
  !               n  / 1  d \ n  sin(x)
  !   j (x) = (-x)  (  - --- )   ------
  !    n             \ x dx /      x
  !
  ! In particular:
  !
  !                               2         4
  !   j (x) = sin(x) / x  ~  1 - x / 6 + o(x )
  !    0
  !
  !                                                    3
  !   j (x) = (sin(x) / x - cos(x)) / x  ~  x / 3 + o(x )
  !    1
  !
  ! Note that
  !
  !    d
  !   --- j (x) = -j (x)
  !    dx  0        1
  !
  use kinds, only: rk
  implicit none
  private

  public :: j0, j1

contains

  elemental function j0 (x) result (f)
    !
    ! sin(x) / x
    !
    implicit none
    real (rk), intent (in) :: x
    real (rk) :: f
    ! *** end of interface ***

    real (rk), parameter :: c1 = -1 / 6.0d0
    real (rk), parameter :: c2 =  1 / 120.0d0
    real (rk), parameter :: c3 = -1 / 5040.0d0
    real (rk), parameter :: c4 =  1 / 362880.0d0
    real (rk), parameter :: c5 = -1 / 39916800.0d0
    real (rk), parameter :: c6 =  1 / 6227020800.0d0 ! > 2^32, use doubles

    if (abs (x) < 0.5) then
      block
         real (rk) :: y

         y = x * x
         f =  1 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * c6)))))
      end block
    else
       f = sin (x) / x
    endif
  end function j0


  elemental function j1 (x) result (f)
    !
    ! (sin(x) / x - cos(x)) / x
    !
    implicit none
    real (rk), intent (in) :: x
    real (rk) :: f
    ! *** end of interface ***

   real (rk), parameter :: c1 = -1 / 10.0d0
   real (rk), parameter :: c2 =  1 / 280.0d0
   real (rk), parameter :: c3 = -1 / 15120.0d0
   real (rk), parameter :: c4 =  1 / 1330560.0d0
   real (rk), parameter :: c5 = -1 / 172972800.0d0

   if (abs (x) < 0.25) then
      block
         real (rk) :: y, sum

         y = x * x
         sum = 1 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * c5))))
         f = x / 3 * sum
      end block
   else
      f = (sin (x) / x - cos (x)) / x
   endif
  end function j1

end module bessel

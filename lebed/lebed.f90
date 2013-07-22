module lebed
  use iso_c_binding, only: ik => c_int, rk => c_double
  implicit none
  private

  public :: genpts

  integer, parameter, public :: order(32) =  &
       [   6,   14,   26,   38,   50,   74,   86,  110,  146,  170, &
         194,  230,  266,  302,  350,  434,  590,  770,  974, 1202, &
        1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, &
        5294, 5810]

  abstract interface
      subroutine ldxxxx (x, y, z, w, n)
       double precision x(*)
       double precision y(*)
       double precision z(*)
       double precision w(*)
       integer :: n
     end subroutine ldxxxx
  end interface

  ! See Lebedev-Laikov.F  for the implementation,  these are interface
  ! declarations only:
  procedure (ldxxxx) :: &
       LD0006, LD0014, LD0026, LD0038, LD0050, &
       LD0074, LD0086, LD0110, LD0146, LD0170, &
       LD0194, LD0230, LD0266, LD0302, LD0350, &
       LD0434, LD0590, LD0770, LD0974, LD1202, &
       LD1454, LD1730, LD2030, LD2354, LD2702, &
       LD3074, LD3470, LD3890, LD4334, LD4802, &
       LD5294, LD5810

contains

  function genpts (m, pts, wts) result (n) bind (c)
    implicit none
    integer (ik), intent (in), value :: m
    real (rk), intent (out) :: pts (3, m), wts (m)
    integer (ik) :: n
    ! *** end of interface ***

    integer :: i

    ! The largest n <= m for which a Lebedev quadrature does exist:
    n = -1
    do i = 1, size (order)
       if (order(i) <= m) then
          n = order(i)
       else
          exit ! the loop
       endif
    end do

    ! FIXME: copy-out here:
    if (n > 0 ) then
       call ld (pts (1, :n), pts (2, :n), pts (3, :n), wts(:n), n)
    else
       pts = huge (pts)
       wts = huge (wts)
    endif
  end function genpts

  subroutine ld (x, y, z, w, n)
    implicit none
    integer, intent (in) :: n
    real (rk), intent (out) :: x(n), y(n), z(n), w(n)
    ! *** end of interface ***

    select case (n)
    case (0006)
       call LD0006 (x, y, z, w, n)
    case (0014)
       call LD0014 (x, y, z, w, n)
    case (0026)
       call LD0026 (x, y, z, w, n)
    case (0038)
       call LD0038 (x, y, z, w, n)
    case (0050)
       call LD0050 (x, y, z, w, n)
    case (0074)
       call LD0074 (x, y, z, w, n)
    case (0086)
       call LD0086 (x, y, z, w, n)
    case (0110)
       call LD0110 (x, y, z, w, n)
    case (0146)
       call LD0146 (x, y, z, w, n)
    case (0170)
       call LD0170 (x, y, z, w, n)
    case (0194)
       call LD0194 (x, y, z, w, n)
    case (0230)
       call LD0230 (x, y, z, w, n)
    case (0266)
       call LD0266 (x, y, z, w, n)
    case (0302)
       call LD0302 (x, y, z, w, n)
    case (0350)
       call LD0350 (x, y, z, w, n)
    case (0434)
       call LD0434 (x, y, z, w, n)
    case (0590)
       call LD0590 (x, y, z, w, n)
    case (0770)
       call LD0770 (x, y, z, w, n)
    case (0974)
       call LD0974 (x, y, z, w, n)
    case (1202)
       call LD1202 (x, y, z, w, n)
    case (1454)
       call LD1454 (x, y, z, w, n)
    case (1730)
       call LD1730 (x, y, z, w, n)
    case (2030)
       call LD2030 (x, y, z, w, n)
    case (2354)
       call LD2354 (x, y, z, w, n)
    case (2702)
       call LD2702 (x, y, z, w, n)
    case (3074)
       call LD3074 (x, y, z, w, n)
    case (3470)
       call LD3470 (x, y, z, w, n)
    case (3890)
       call LD3890 (x, y, z, w, n)
    case (4334)
       call LD4334 (x, y, z, w, n)
    case (4802)
       call LD4802 (x, y, z, w, n)
    case (5294)
       call LD5294 (x, y, z, w, n)
    case (5810)
       call LD5810 (x, y, z, w, n)
    end select
  end subroutine ld
end module lebed

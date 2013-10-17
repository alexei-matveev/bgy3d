module fft
  use kinds, only: rk
  implicit none
  private

  real (rk), parameter, public :: pi = 4 * atan (1.0_rk)

  !
  ! These should  obey: ab  = (2π)³, a²b  = (2π)³. The  first equality
  ! guarantees that the forward- and backward FT are mutually inverse,
  ! whereas  the second  equation  will make  the convolution  theorem
  ! factor-less as in FT(f * g) = FT(f) FT(g):
  !
  real (rk), parameter, public :: FT_FW = 1           ! == a
  real (rk), parameter, public :: FT_BW = (2 * pi)**3 ! == b

  !
  ! Functions operated upon by fourier(), integrate() and their fellow
  ! siblings depend on radius "r" only and are represented by an array
  ! such as f(1:n) where
  !
  !   f  = f(dr * [2 * i - 1] / 2),  1 <= i <= n
  !    i
  !
  ! and  dr  =  1  for  the  purposes  of  the  current  module.   The
  ! corresponding Fourier  transform, g =  FT(f), is represented  on a
  ! similarly spaced k-grid,
  !
  !   g  = g(dk * [2 * k - 1] / 2),  1 <= k <= n
  !    k
  !
  ! albeit with dk related to dr by the following expression:
  !
  !   dr * dk = 2π / 2n
  !
  public :: fourier_rows        ! f(*, 1:n) -> g(*, 1:n)
! public :: fourier_cols        ! f(1:n, *) -> g(1:n, *)
  public :: integrate           ! f(1:n) -> scalar
  !
  ! *** END OF INTERFACE ***
  !

  !
  ! This is a concrete function, implemented in C:
  !
  interface
     subroutine rism_dst (n, out, in) bind (c)
       !
       ! Performs DST. In FFTW terms this is RODFT11 (or DST-IV) which
       ! is self inverse up to a normalization factor.
       !
       ! void rism_dst (size_t n, double out[n], const double in[n])
       !
       ! See ./rism-dst.c
       !
       use iso_c_binding, only: c_size_t, c_double
       implicit none
       integer (c_size_t), intent (in), value :: n
       real (c_double), intent (out) :: out(n)
       real (c_double), intent (in) :: in(n)
     end subroutine rism_dst

     subroutine rism_dst_columns (m, n, buf) bind (c)
       !
       ! Performs DST of size n in-place for each of the m columns. In
       ! FFTW terms this is RODFT11  (or DST-IV) which is self inverse
       ! up to a normalization factor.
       !
       ! See ./rism-dst.c
       !
       use iso_c_binding, only: c_int, c_double
       implicit none
       integer (c_int), intent (in), value :: m, n
       real (c_double), intent (inout) :: buf(n, m)
     end subroutine rism_dst_columns

     subroutine rism_dst_rows (n, m, buf) bind (c)
       !
       ! Performs DST  of size n in-place  for each of the  m rows. In
       ! FFTW terms this is RODFT11  (or DST-IV) which is self inverse
       ! up to a normalization factor.
       !
       ! See ./rism-dst.c
       !
       use iso_c_binding, only: c_int, c_double
       implicit none
       integer (c_int), intent (in), value :: m, n
       real (c_double), intent (inout) :: buf(m, n)
     end subroutine rism_dst_rows
  end interface

contains

  function fourier_rows (f) result (g)
    implicit none
    real (rk), intent (in) :: f(:, :, :)
    real (rk) :: g(size (f, 1), size (f, 2), size (f, 3))
    ! *** end of interface ***

    integer :: p, i, j, n
    real (rk) :: fac
    real (rk), parameter :: one = 1

    n = size (f, 3)
    !
    ! We use  RODFT11 (DST-IV) that is  "odd around j =  -0.5 and even
    ! around j  = n - 0.5".   Here we use integer  arithmetics and the
    ! identity (2 * j - 1) / 2 == j - 0.5.
    !
    !$omp parallel do private(fac, p, i, j)
    do p = 1, n
       fac = 2 * n * (2 * p - 1)
       do j = 1, size (f, 2)
          do i = 1, size (f, 1)
             g(i, j, p) =  f(i, j, p) * fac
          enddo
       enddo
    enddo
    !$omp end parallel do

    call dst_rows (g)

    !$omp parallel do private(fac, p, i, j)
    do p = 1, n
       fac = one / (2 * p - 1)
       do j = 1, size (g, 2)
          do i = 1, size (g, 1)
             g(i, j, p) = g(i, j, p) * fac
          enddo
       enddo
    enddo
    !$omp end parallel do
  end function fourier_rows


  function fourier_cols (f) result (g)
    implicit none
    real (rk), intent (in) :: f(:, :, :)
    real (rk) :: g(size (f, 1), size (f, 2), size (f, 3))
    ! *** end of interface ***

    integer :: p, i, j, n

    n = size (f, 1)
    !
    ! We use  RODFT11 (DST-IV) that is  "odd around j =  -0.5 and even
    ! around j  = n - 0.5".   Here we use integer  arithmetics and the
    ! identity (2 * j - 1) / 2 == j - 0.5.
    !
    !$omp parallel workshare private(p, i, j)
    forall (p = 1:n, i = 1:size (f, 2), j = 1:size (f, 3))
       g(p, i, j) =  f(p, i, j) * (2 * n * (2 * p - 1))
    end forall
    !$omp end parallel workshare

    call dst_columns (g)

    !$omp parallel workshare private(p, i, j)
    forall (p = 1:n, i = 1:size (g, 2), j = 1:size (g, 3))
       g(p, i, j) = g(p, i, j) / (2 * p - 1)
    end forall
    !$omp end parallel workshare
  end function fourier_cols


  function fourier (f) result (g)
    implicit none
    real (rk), intent (in) :: f(:)
    real (rk) :: g(size (f))
    ! *** end of interface ***

    integer :: i, n

    n = size (f)

    !
    ! We use  RODFT11 (DST-IV) that is  "odd around j =  -0.5 and even
    ! around j  = n - 0.5".   Here we use integer  arithmetics and the
    ! identity (2 * j - 1) / 2 == j - 0.5.
    !
    forall (i = 1:n)
       g(i) = f(i) * (2 * i - 1)
    end forall

    g = 2 * n * dst (g)

    forall (i = 1:n)
       g(i) = g(i) / (2 * i - 1)
    end forall
  end function fourier


  pure function integrate (f) result (g)
    !
    ! Approximates 4π  ∫f(r)r²dr.  Should  be "related" to  the k  = 0
    ! component of the  FT. This is also the  only reason the function
    ! is put into this module.
    !
    implicit none
    real (rk), intent (in) :: f(:)
    real (rk) :: g
    ! *** end of interface ***

    integer :: i, n

    n = size (f)

    g = 0.0
    do i = n, 1, -1
       g = g + f(i) * (2 * i - 1)**2 / 4
    enddo
    g = 4 * pi * g
  end function integrate


  subroutine dst_columns (f)
    use iso_c_binding, only: c_int
    implicit none
    real (rk), intent (inout) :: f(:, :, :) ! ~ (n, m)
    ! *** end of interface ***

    integer (c_int) :: m, n

    ! cast to c_int
    n = size (f, 1)
    m = size (f) / n

    call rism_dst_columns (m, n, f)
  end subroutine dst_columns


  subroutine dst_rows (f)
    use iso_c_binding, only: c_int
    implicit none
    real (rk), intent (inout) :: f(:, :, :) ! ~ (m, n)
    ! *** end of interface ***

    integer (c_int) :: m, n

    ! cast to c_int
    n = size (f, 3)
    m = size (f) / n

    call rism_dst_rows (n, m, f)
  end subroutine dst_rows


  function dst (a) result (b)
    use iso_c_binding, only: c_size_t
    implicit none
    real (rk), intent (in) :: a(:)
    real (rk) :: b(size (a))
    ! *** end of interface ***

    integer (c_size_t) :: n

    ! cast to size_t
    n = size (a)

    call rism_dst (n, b, a)
  end function dst


  subroutine test_dst (n)
    implicit none
    integer, intent (in) :: n
    ! *** end of interface ***

    real (rk) :: a(n)

    call random_number (a)

    ! RODFT11 (DST-IV) is self inverse up to a normalization factor:
    if (maxval (abs (a - dst (dst (a)) / (2 * n))) > 1.0e-10) then
       print *, "diff=", maxval (abs (a - dst (dst (a)) / (2 * n)))
       stop "unnormalized DST does not match"
    endif
  end subroutine test_dst


  subroutine test_ft (rmax, n)
    implicit none
    integer, intent (in) :: n
    real (rk), intent (in) :: rmax
    ! *** end of interface ***

    real (rk) :: r(n), k(n), f(n), g(n), h(n)
    real (rk) :: dr, dk

    integer :: i

    !
    ! The 3d analytical unitary FT of the function
    !
    !   f = exp(-r²/2)
    !
    ! is the function
    !
    !   g = exp(-k²/2), FIXME!
    !
    dr = rmax / n
    dk = pi / rmax
    forall (i = 1:n)
       r(i) = (i - 0.5) * dr
       k(i) = (i - 0.5) * dk
    end forall

    ! Gaussian:
    f = exp (-r**2 / 2)

    ! Unit "charge", a "fat" delta function:
    f = f / (integrate (f) * dr**3)

    ! Forward transform:
    g = fourier (f) * (dr**3 / FT_FW)

    ! Backward transform:
    h = fourier (g) * (dk**3 / FT_BW)

    print *, "# norm (f )^2 =", integrate (f**2) * dr**3
    print *, "# norm (g )^2 =", integrate (g**2) * dk**3
    print *, "# norm (f')^2 =", integrate (h**2) * dr**3

    print *, "# |f - f'| =", maxval (abs (f - h))

    print *, "# int (f') =", integrate (h) * dr**3
    print *, "# int (f ) =", integrate (f) * dr**3

    ! This should correspond  to the convolution (f *  f) which should
    ! be again a gaussian, twice as "fat":
    h = g * g
    h = fourier (h) * (dk**3 / FT_BW)

    print *, "# int (h ) =", integrate (h) * dr**3

    ! Compare width as <r^2>:
    print *, "# sigma (f) =", integrate (r**2 * f) * dr**3
    print *, "# sigma (h) =", integrate (r**2 * h) * dr**3

    ! print *, "# n=", n
    ! print *, "# r, f, k, g, h = (f*f)"
    ! do i = 1, n
    !    print *, r(i), f(i), k(i), g(i), h(i)
    ! enddo
  end subroutine test_ft
end module fft

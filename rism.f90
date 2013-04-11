module rism
  use iso_c_binding, only: rk => c_double
  implicit none
  private

  public :: rism_main           ! ()
  real(rk), parameter, public :: pi = 4 * atan (1.0_rk)
  ! *** end of interface ***

  !
  ! This defines two  iterator interfaces x -> dx  for use in fixpoint
  ! non-linear problems.  One is  the type of Fortran closures another
  ! is the corresponding type of C closures. Such an iterator function
  ! is supposed to return array of zeros (ok, rather small numbers) at
  ! convergence.
  !
  abstract interface
     function f_iterator (x) result (dx)
       import rk
       implicit none
       real(rk), intent(in) :: x(:)
       real(rk) :: dx(size (x))
     end function f_iterator

     subroutine c_iterator (ctx, n, x, dx) bind (c)
       use iso_c_binding, only: c_ptr, c_int, c_double
       implicit none
       type(c_ptr), intent(in), value :: ctx
       integer(c_int), intent(in), value :: n
       real(c_double), intent(in) :: x(n)
       real(c_double), intent(out) :: dx(n)
     end subroutine c_iterator
  end interface

  !
  ! These are concrete functions, implemented in C:
  !
  interface
     subroutine rism_dst (n, out, in) bind (c)
       !
       ! Performs DST. In FFTW terms this is RODFT11 (or DST-IV) which
       ! is is self inverse up to a normalization factor.
       !
       ! void rism_dst (size_t n, double out[n], const double in[n])
       !
       ! See ./rism-dst.c
       !
       use iso_c_binding, only: c_size_t, c_double
       implicit none
       integer(c_size_t), intent(in), value :: n
       real(c_double), intent(out) :: out (n)
       real(c_double), intent(in) :: in (n)
     end subroutine rism_dst

     subroutine rism_snes (ctx, f, n, x) bind (c)
       !
       ! void rism_snes (void *ctx, ArrayFunc f,
       !                 int n, real x_[n])
       !
       ! using procedure(c_iterator) or in C-lang ArrayFunc:
       !
       ! typedef void (*ArrayFunc) (void *ctx,
       !                            int n, const real x[n], real r[n]);
       !
       ! See ./bgy3d-snes.c
       !
       use iso_c_binding, only: c_ptr, c_int, c_double
       implicit none
       type(c_ptr), intent (in), value :: ctx
       procedure(c_iterator) :: f
       integer(c_int), intent (in), value :: n
       real(c_double), intent (inout) :: x(n)
     end subroutine rism_snes
  end interface

  !
  ! A pointer to to the structure of this type will serve as a closure
  ! context to  be passed to  the C-world. Upon callback  the fourtran
  ! sub  implementing a  procedure(c_iterator) can  use  the procedure
  ! pointer to perform the actual  work. I wish closures could be made
  ! simpler than that.
  !
  type context
     procedure (f_iterator), pointer, nopass :: f
  end type context

contains

  subroutine rism_main () bind (c)
    implicit none
    ! *** end of interface ***

    integer :: i
    do i = 1, 10**0
       call test_dst (2**12)
    end do
    call test_ft (rmax=10.d0, n=2**14)
    ! stop

    ! FIXME: why 4π?
    ! call rism1d (rho=1.054796d0, beta=0.261610d0, rmax=10.0d0, n=256)
    call rism1d (rho=(4 * pi * 1.054796d0), beta=0.261610d0, rmax=20.0d0, n=2**14)
  end subroutine rism_main

  subroutine rism1d (rho, beta, rmax, n)
    implicit none
    real(rk), intent(in) :: rho, beta, rmax
    integer, intent(in) :: n
    ! *** end of interface ***

    real(rk), parameter :: half = 0.5, alpha = 0.02
    real(rk) :: r(n), k(n), v(n), t(n), c(n), g(n)
    real(rk) :: dr, dk
    integer :: i

    call print_info (rho=rho, beta=beta)

    ! dr * dk = 2π/2n:
    dr = rmax / n
    dk = pi / rmax
    forall (i = 1:n)
       r(i) = (i - half) * dr
       k(i) = (i - half) * dk
    end forall

    ! LJ potential, sigma=1, epsilon=1:
    v = lj (r)

    ! Intirial guess:
    t = 0.0

    ! Find t such that iterate_t  (t) == 0. FIXME: passing an internal
    ! function  as a  callback is  an F2008  feature. GFortran  4.3 on
    ! Debian Lenny does not support that:
    call snes_default (iterate_t, t)

    ! Do not assume c has a meaningfull value, it was overwritten with
    ! c(k):
    c = closure_hnc (beta, v, t)
    g = 1 + c + t

    print *, "# rho=", rho, "beta=", beta, "n=", n
    print *, "# r, v, t, c, g"
    do i = 1, n
       print *, r(i), v(i), t(i), c(i), g(i)
    enddo

  contains

    function iterate_t (t) result (dt)
      !
      ! Closure over  host variables: r, k,  dr, dk, v,  c, beta, rho,
      ! ... Implements procedure(f_iterator).
      !
      implicit none
      real(rk), intent(in) :: t(:)
      real(rk) :: dt(size (t))
      ! *** end of interface ***

      c = closure_hnc (beta, v, t)

      ! Forward FT via DST:
      c = dst (r * c) / k * (dr / sqrt (2 * pi))

      ! OZ equation, involves "convolutions", take care of the
      ! normalization here:
      dt = oz_equation_c_t (rho, c)

      ! Inverse FT via DST:
      dt = dst (k * dt) / r * (dk / sqrt (2 * pi))

      ! Return the increment that vanishes at convergence:
      dt = dt - t
    end function iterate_t
  end subroutine rism1d


  subroutine snes_picard (f, x)
    !
    ! Simple Picard iteration
    !
    !   x    = x + α f(x )
    !    n+1    n       n
    !
    implicit none
    procedure(f_iterator) :: f  ! (x) -> dx
    real(rk), intent(inout) :: x(:)
    ! *** end of interface ***

    real(rk), parameter :: alpha = 0.02
    real(rk) :: dx(size (x))
    real(rk) :: diff
    integer :: iter
    logical :: converged

    iter = 0
    converged = .false.
    do while (.not. converged)
       iter = iter + 1

       dx = f (x)

       diff = maxval (abs (dx))
       converged = (diff < 1.0e-12)

       x = x + alpha * dx
       print *, "# iter=", iter, "diff=", diff
    end do
  end subroutine snes_picard


  subroutine snes_default (f, x)
    !
    ! Delegates the actual work to Petsc by way of C-func rism_snes().
    !
    use iso_c_binding
    implicit none
    procedure(f_iterator) :: f  ! (x) -> dx
    real(rk), intent(inout) :: x(:)
    ! *** end of interface ***

    type(context), target :: f_ctx
    type(c_ptr) :: ctx

    f_ctx % f => f
    ctx = c_loc (f_ctx)

    call rism_snes (ctx, iterator, size (x), x)
  end subroutine snes_default


  subroutine iterator (ctx, n, x, dx) bind (c)
    !
    ! Implements  procedure(c_iterator)  and  will  be passed  to  the
    ! rism_snes()   together   with    the   suitable   context.   See
    ! snes_default().
    !
    use iso_c_binding
    implicit none
    type(c_ptr), intent(in), value :: ctx
    integer(c_int), intent(in), value :: n
    real(rk), intent(in) :: x(n)
    real(rk), intent(out) :: dx(n)
    ! *** end of interface ***

    type(context), pointer :: f_ctx

    ! We  dont  so any  work  ourselves,  just  extract a  pointer  to
    ! procedure(f_iterator) an let it do the rest:
    call c_f_pointer (ctx, f_ctx)

    dx = f_ctx % f (x)
  end subroutine iterator

  !
  ! 1)  Hypernetted Chain  (HNC)  closure relation  to compute  direct
  ! correlation function c  in real space.  See OZ  equation below for
  ! the   second  relation   between  two   unknowns.    The  indirect
  ! correlation γ = h - c is denoted by latin "t" in other sources. We
  ! will  use  that to  avoid  greek  identifiers  and confusion  with
  ! distribution functions:
  !
  !   c := exp (-βv + γ) - 1 - γ
  !
  elemental function closure_hnc (beta, v, t) result (c)
    implicit none
    real(rk), intent(in) :: beta, v, t
    real(rk) :: c
    ! *** end of interface ***

    ! exp (-beta * v + t) - 1.0 - t:
    c = exp (-beta * v + t) - 1 - t
  end function closure_hnc


  !
  ! Use the k-representation of Ornstein-Zernike (OZ) equation
  !
  !   h = c + ρ c * h
  !
  ! to compute γ =  h - c form c:
  !
  !                -1                -1   2
  !   γ =  (1 - ρc)   c - c = (1 - ρc)  ρc
  !
  ! If you scale c by h3 beforehand  or pass rho' = rho * h3 and scale
  ! the result  by h3 in addition,  you will compute  exactly what the
  ! older version of the function did:
  !
  elemental function oz_equation_c_t (rho, c) result (t)
    implicit none
    real(rk), intent(in) :: rho, c
    real(rk) :: t
    ! *** end of interface ***

    ! Alternatively: t = c / (1 - rho * c) - c
    t = rho * (c * c) / (1 - rho * c)
  end function oz_equation_c_t


  elemental function lj (r) result (f)
    !
    ! To be called as in eps * lj (r / sigma)
    !
    implicit none
    real(rk), intent(in) :: r   ! r / sigma, in general
    real(rk) :: f
    ! *** end of interfce ***

    real(rk) :: sr6

    sr6 = 1 / r**6

    f = 4 * sr6 * (sr6 - 1)
  end function lj


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


  function ndst (a) result (b)
    implicit none
    real(rk), intent(in) :: a(:)
    real(rk) :: b(size (a))
    ! *** end of interface ***

    real(rk) :: norm

    norm = 2 * size (b)         ! cast to real
    b = dst (a) / sqrt (norm)
  end function ndst


  subroutine test_dst (n)
    implicit none
    integer, intent(in) :: n
    ! *** end of interface ***

    real(rk) :: a(n)

    call random_number (a)
    a = lj (a + 0.5)

    ! RODFT11 (DST-IV) is self inverse up to a normalization factor:
    if (maxval (abs (a - dst (dst (a)) / (2 * n))) > 1.0e-10) then
       print *, "diff=", maxval (abs (a - dst (dst (a)) / (2 * n)))
       stop "unnormalized DST does not match"
    endif

    if (maxval (abs (a - ndst (ndst (a)))) > 1.0e-10) then
       print *, "diff=", maxval (abs (a - ndst (ndst (a))))
       stop "normalized DST does not match"
    endif
  end subroutine test_dst


  subroutine test_ft (rmax, n)
    implicit none
    integer, intent(in) :: n
    real(rk), intent(in) :: rmax
    ! *** end of interface ***

    real(rk) :: r(n), k(n), f(n), g(n), h(n)
    real(rk) :: dr, dk
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

    ! L2-normalized:
    f = f / sqrt (sum ((r * f)**2) * 4 * pi * dr)

    ! Unitary forward transform:
    g = dst (r * f) / k * (dr / sqrt (2 * pi))

    ! Unitary backward transform:
    h = dst (k * g) / r * (dk / sqrt (2 * pi))

    print *, "# int(f)=", sum ((r * f)**2) * 4 * pi * dr
    print *, "# int(g)=", sum ((k * g)**2) * 4 * pi * dk
    print *, "# int(h)=", sum ((r * h)**2) * 4 * pi * dr
    print *, "# |f - h|=", maxval (abs (f - h))
    print *, "# sigma(f)=", sum (r**2 * (r * f)**2) * 4 * pi * dr
    print *, "# sigma(g)=", sum (k**2 * (k * g)**2) * 4 * pi * dk

    ! print *, "# n=", n
    ! print *, "# r, f, k, g"
    ! do i = 1, n
    !    print *, r(i), f(i), k(i), g(i)
    ! enddo
  end subroutine test_ft


  subroutine print_info (rho, beta)
    !
    ! van der Hoef (Ref. [6] Eqs. 25 and 26):
    !
    !          -1/4                          2         3         4         5
    ! ρ      =β    [0.92302-0.09218β+0.62381β -0.82672β +0.49124β -0.10847β ]
    !  solid
    !          -1/4                          2         3         4         5
    ! ρ      =β    [0.91070-0.25124β+0.85861β -1.08918β +0.63932β -0.14433β ]
    !  liquid
    !
    ! Mastny and de Pablo (Ref [7] Eqs. 20 and 21):
    !
    !          -1/4                             2          3          4          5
    ! ρ      =β    [0.908629-0.041510β+0.514632β -0.708590β +0.428351β -0.095229β ]
    !  solid
    !          -1/4                          2         3         4         5
    ! ρ      =β    [0.90735-0.27120β+0.91784β -1.16270β +0.68012β -0.15284β ]
    !  liquid
    !
    implicit none
    real(rk), intent(in) :: rho, beta
    ! *** end of interface ***

    real(rk), parameter :: psol_hoef(6) = [0.92302d0, -0.09218d0, +0.62381d0, -0.82672d0, +0.49124d0, -0.10847d0]
    real(rk), parameter :: pliq_hoef(6) = [0.91070d0, -0.25124d0, +0.85861d0, -1.08918d0, +0.63932d0, -0.14433d0]
    real(rk), parameter :: psol_mast(6) = [0.908629d0, -0.041510d0, +0.514632d0, -0.708590d0, +0.428351d0, -0.095229d0]
    real(rk), parameter :: pliq_mast(6) = [0.90735d0, -0.27120d0, +0.91784d0, -1.16270d0, +0.68012d0, -0.15284d0]

    print *, "# rho =", rho, "rs =", (4 * pi * rho / 3)**(-1d0/3d0)
    print *, "# beta =", beta, "T =", 1 / beta
    print *, "#"
    print *, "# At this temperature ..."
    print *, "#"
    print *, "# According to Hoef"
    print *, "# rho solid  =", beta**(-0.25d0) * poly (psol_hoef, beta)
    print *, "# rho liquid =", beta**(-0.25d0) * poly (pliq_hoef, beta)
    print *, "#"
    print *, "# According to Mastny and de Pablo"
    print *, "# rho solid  =", beta**(-0.25d0) * poly (psol_mast, beta)
    print *, "# rho liquid =", beta**(-0.25d0) * poly (pliq_mast, beta)
    print *, "#"
  end subroutine print_info

  function poly (p, x) result (y)
    implicit none
    real(rk), intent(in) :: p(0:), x
    real(rk) :: y
    ! *** end of interface ***

    integer :: n

    y = 0.0
    do n = 0, size (p) - 1
       y = y + p(n) * x**n
    enddo
  end function poly
end module rism


program rism_prog
  use rism, only: rism_main
  implicit none

  call rism_main ()
end program rism_prog

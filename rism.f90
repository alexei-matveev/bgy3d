module rism
  use iso_c_binding, only: c_int, c_double, c_char
  implicit none
  private

  integer, parameter :: ik = c_int, rk = c_double

  public :: rism_main
  public :: bgy3d_problem_data

  real (rk), parameter, public :: pi = 4 * atan (1.0_rk)

  ! Keep this in sync with bgy3d-solutes.h:
  type, public, bind (c) :: site
     character (kind=c_char, len=5) :: name ! atom types. What are they used for?
     real (c_double) :: x(3)                ! coordinates
     real (c_double) :: sigma               ! sigma for LJ
     real (c_double) :: epsilon             ! epsilon for LJ
     real (c_double) :: charge              ! charge
  end type site

  ! Keep this in sync with bgy3d.h:
  type, public, bind (c) :: problem_data
     real (c_double) :: interval(2) ! min and max of the domain: 3d-box*/
     real (c_double) :: h(3)        ! mesh width
     real (c_double) :: beta        ! inverse temperature, 1/kT
     real (c_double) :: rho         ! solvent density
     integer (c_int) :: N(3), N3    ! global Grid size

     ! Other staff  that was retrieved by the  solvers themselves from
     ! the (Petsc) environment:
     real (c_double) :: lambda   ! Mixing parameter.
     real (c_double) :: damp     ! Scaling factor.
     integer (c_int) :: max_iter ! Maximal number of iterations.
     real (c_double) :: norm_tol ! Convergence threshold.
     real (c_double) :: zpad     ! FIXME: ???
  end type problem_data

  !
  ! *** END OF INTERFACE ***
  !

  !                              2
  ! These should obey: ab = 2π, a b = 1:
  !
  real (rk), parameter :: a = 1 / (2 * pi), b = (2 * pi)**2

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
       real (rk), intent (in) :: x(:, :, :)
       real (rk) :: dx(size (x, 1), size (x, 2), size (x, 3))
     end function f_iterator

     subroutine c_iterator (ctx, n, x, dx) bind (c)
       use iso_c_binding, only: c_ptr, c_int, c_double
       implicit none
       type (c_ptr), intent (in), value :: ctx
       integer (c_int), intent (in), value :: n
       real (c_double), intent (in) :: x(n)
       real (c_double), intent (out) :: dx(n)
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
       integer (c_size_t), intent (in), value :: n
       real (c_double), intent (out) :: out(n)
       real (c_double), intent (in) :: in(n)
     end subroutine rism_dst

     subroutine rism_snes (ctx, f, n, x) bind (c)
       !
       ! void rism_snes (void *ctx, ArrayFunc f,
       !                 int n, real x_[n])
       !
       ! using procedure (c_iterator) or in C-lang ArrayFunc:
       !
       ! typedef void (*ArrayFunc) (void *ctx,
       !                            int n, const real x[n], real r[n]);
       !
       ! See ./bgy3d-snes.c
       !
       use iso_c_binding, only: c_ptr, c_int, c_double
       implicit none
       type (c_ptr), intent (in), value :: ctx
       procedure (c_iterator) :: f
       integer (c_int), intent (in), value :: n
       real (c_double), intent (inout) :: x(n)
     end subroutine rism_snes

     function bgy3d_problem_data () result (pd) bind (c)
       import problem_data
       implicit none
       type (problem_data) :: pd
     end function bgy3d_problem_data
  end interface

  !
  ! A pointer to to the structure of this type will serve as a closure
  ! context to  be passed to  the C-world. Upon callback  the fourtran
  ! sub  implementing a  procedure(c_iterator) can  use  the procedure
  ! pointer to perform the actual  work. I wish closures could be made
  ! simpler than that.
  !
  type context
     integer :: shape(3)
     procedure (f_iterator), pointer, nopass :: f
  end type context

contains

  subroutine rism_main (pd, m, sites) bind (c)
    use iso_c_binding, only: c_int
    implicit none
    type (problem_data), intent (in) :: pd ! no VALUE!
    integer (c_int), intent (in), value :: m
    type (site), intent (in) :: sites(m)
    ! *** end of interface ***

    integer :: i, n
    real (rk) :: rmax
    real (rk) :: rho(m)

    rho = pd % rho              ! all site densities are the same
    rmax = 0.5 * (pd % interval(2) - pd % interval(1))
    n = maxval (pd % n)

    print *, "# rho=", pd % rho
    print *, "# beta=", pd % beta
    print *, "# L=", rmax
    print *, "# n=", n
    print *, "# Sites:"
    do i = 1, m
       print *, "#", i, &
            &        pad (sites(i) % name), &
            &        sites(i) % sigma, &
            &        sites(i) % epsilon, &
            &        sites(i) % charge, &
            &        rho(i)
    enddo

    ! This is applicable to LJ only, and should take reduced
    ! density/temperature:
    ! call print_info (rho = pd % rho, beta = pd % beta)

    call rism1d (n, rmax, beta= pd % beta, rho= rho, sites= sites)
  end subroutine rism_main


  subroutine rism1d (n, rmax, beta, rho, sites)
    implicit none
    integer, intent (in) :: n            ! grid size
    real (rk), intent (in) :: rmax       ! cell size
    real (rk), intent (in) :: beta       ! inverse temp
    real (rk), intent (in) :: rho(:)     ! (m)
    type (site), intent (in) :: sites(:) ! (m)
    ! *** end of interface ***

    real (rk), parameter :: half = 0.5

    ! Pair quantities. FIXME: they are symmetric, one should use that:
    real (rk), dimension (n, size (sites), size (sites)) :: v, t, c, g

    ! Radial grids:
    real (rk) :: r(n), dr, dk

    integer :: i, m

    m = size (sites)

    ! dr * dk = 2π/2n:
    dr = rmax / n
    dk = pi / rmax
    forall (i = 1:n)
       r(i) = (i - half) * dr
    end forall

    ! Tabulate pairwise potential v() on the r-grid:
    call force_field (sites, r, v)

    ! Intitial guess:
    t = 0.0

    ! Find t such that iterate_t  (t) == 0. FIXME: passing an internal
    ! function  as a  callback is  an F2008  feature. GFortran  4.3 on
    ! Debian Lenny does not support that:
    call snes_default (iterate_t, t)

    ! Do not assume c has a meaningfull value, it was overwritten with
    ! c(k):
    c = closure_hnc (beta, v, t)
    g = 1 + c + t

    ! Done with it, print results:
    block
       integer :: p, i, j
       print *, "# rho=", rho, "beta=", beta, "n=", n
       print *, "# r, and v, t, c, g, each for m * (m + 1) / 2 pairs"
       do p = 1, n
          write (*, *) r(p), &
               &     ((v(p, i, j), i=1,j), j=1,m), &
               &     ((t(p, i, j), i=1,j), j=1,m), &
               &     ((c(p, i, j), i=1,j), j=1,m), &
               &     ((g(p, i, j), i=1,j), j=1,m)

       enddo
    end block

  contains

    function iterate_t (t) result (dt)
      !
      ! Closure over  host variables: r, k,  dr, dk, v,  c, beta, rho,
      ! ... Implements procedure(f_iterator).
      !
      implicit none
      real (rk), intent (in) :: t(:, :, :) ! (n, m, m)
      real (rk) :: dt(size (t, 1), size (t, 2), size (t, 3))
      ! *** end of interface ***

      integer :: i, j

      c = closure_hnc (beta, v, t)

      ! Forward FT via DST:
      do j = 1, m
         do i = 1, m
            c(:, i, j) = fourier (c(:, i, j)) * (dr**3 * (n/pi) / a)
         enddo
      enddo

      ! OZ equation, involves "convolutions", take care of the
      ! normalization here:
      dt = oz_equation_c_t (rho, c)

      ! Inverse FT via DST:
      do j = 1, m
         do i = 1, m
            dt(:, i, j) = fourier (dt(:, i, j)) * (dk**3 * (n/pi) / b)
         enddo
      enddo

      ! Return the increment that vanishes at convergence:
      dt = dt - t
    end function iterate_t
  end subroutine rism1d


  subroutine force_field (sites, r, v)
    implicit none
    type (site), intent (in) :: sites(:) ! (m)
    real (rk), intent (in) :: r(:)       ! (n)
    real (rk), intent (out) :: v(:, :, :) ! (n, m, m)
    ! *** end of inteface ***

    real (rk) :: epsilon, sigma, charge
    integer :: i, j, m
    type (site) :: a, b

    m = size (sites)

    ! LJ potential:
    do j = 1, m
       b = sites(j)
       do i = 1, m
          a = sites(i)

          epsilon = sqrt (a % epsilon * b % epsilon)
          sigma = (a % sigma + b % sigma) / 2
          charge = a % charge * b % charge

          if (sigma /= 0.0) then
             v(:, i, j) = epsilon * lj (r / sigma)
          else
             v(:, i, j) = 0.0
          endif
       enddo
    enddo
  end subroutine force_field


  subroutine snes_picard (f, x)
    !
    ! Simple Picard iteration
    !
    !   x    = x + α f(x )
    !    n+1    n       n
    !
    implicit none
    procedure (f_iterator) :: f  ! (x) -> dx
    real (rk), intent (inout) :: x(:, :, :)
    ! *** end of interface ***

    real (rk), parameter :: alpha = 0.02
    real (rk) :: dx(size (x, 1), size (x, 2), size (x, 3))
    real (rk) :: diff
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
    use iso_c_binding, only: c_ptr, c_loc
    implicit none
    procedure (f_iterator) :: f  ! (x) -> dx
    real (rk), intent (inout) :: x(:, :, :)
    ! *** end of interface ***

    type (context), target :: f_ctx
    type (c_ptr) :: ctx

    f_ctx % f => f
    f_ctx % shape = shape (x)
    ctx = c_loc (f_ctx)

    call rism_snes (ctx, iterator, size (x), x)
  end subroutine snes_default


  subroutine iterator (ctx, n, x, dx) bind (c)
    !
    ! Implements  procedure(c_iterator)  and  will  be passed  to  the
    ! rism_snes()   together   with   the   suitable   context.    See
    ! snes_default().
    !
    use iso_c_binding, only: c_f_pointer, c_ptr, c_int
    implicit none
    type (c_ptr), intent (in), value :: ctx
    integer (c_int), intent (in), value :: n
    real (rk), intent (in) :: x(n)
    real (rk), intent (out) :: dx(n)
    ! *** end of interface ***

    type (context), pointer :: f_ctx
    integer :: shape(3)

    ! We  dont  so any  work  ourselves,  just  extract a  pointer  to
    ! procedure(f_iterator) an let it do the rest:
    call c_f_pointer (ctx, f_ctx)

    shape = f_ctx % shape

    block
       real (rk) :: y(shape(1), shape(2), shape(3))
       real (rk) :: dy(shape(1), shape(2), shape(3))

       y = reshape (x, shape)

       dy = f_ctx % f (y)

       dx = reshape (dy, [n])
    end block
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
    real (rk), intent (in) :: beta, v, t
    real (rk) :: c
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
  function oz_equation_c_t (rho, C) result (T)
    implicit none
    real (rk), intent (in) :: rho(:), C(:, :, :) ! (m), (n, m, m)
    real (rk) :: T(size (C, 1), size (C, 2), size (C, 3))
    ! *** end of interface ***

    integer :: i

    if (size (rho) == 1) then
       do i = 1, size (C, 1)
          T(i, 1, 1) = oz_equation_c_t_1x1 (rho(1), C(i, 1, 1))
       enddo
    else
       do i = 1, size (C, 1)
          T(i, :, :) = oz_equation_c_t_MxM (rho(:), C(i, :, :))
       enddo
    endif
  end function oz_equation_c_t


  elemental function oz_equation_c_t_1x1 (rho, c) result (t)
    implicit none
    real (rk), intent (in) :: rho, c
    real (rk) :: t
    ! *** end of interface ***

    !
    ! The actual (matrix) expression
    !
    !                  -1
    !   t = (1 - ρ * c)  * c - c
    !
    ! simplifies in the 1x1 case to this:
    !
    t = rho * (c * c) / (1 - rho * c)
  end function oz_equation_c_t_1x1


  function oz_equation_c_t_MxM (rho, C) result (T)
    !
    ! So far  rho is the same for  all sites, we may  have mixed left-
    ! and right-side multiplies.
    !
    implicit none
    real (rk), intent (in) :: rho(:)  ! (m)
    real (rk), intent (in) :: C(:, :) ! (m, m)
    real (rk) :: T(size (rho), size (rho)) ! (m, m)
    ! *** end of interface ***

    real (rk), dimension (size (rho), size (rho)) :: H
    integer :: i, j, m

    m = size (rho)

    ! T :=  1 - ρC, H  := C. The  latter will be owerwritten  with the
    ! real H after  solving the linear equations. The  output matrix T
    ! is used here as a free work array:
    do j = 1, m
       do i = 1, m
          H(i, j) = C(i, j)
          T(i, j) = delta (i, j) - rho(i) * C(i, j)
       enddo
    enddo

    !
    ! Solving the linear equation makes  H have the literal meaning of
    ! the total correlation matrix (input T is destroyed):
    !
    !       -1                -1
    ! H := T   * H == (1 - ρC)   * C
    !
    call sles (m, T, H)

    !
    ! The  same  effect  is  achieved  in  1x1  version  of  the  code
    ! differently:
    !
    ! T := H - C
    !
    T = H - C

  contains

    function delta (i, j) result (d)
      implicit none
      integer, intent (in) :: i, j
      integer :: d
      ! *** end of interface ***

      if (i == j) then
         d = 1
      else
         d = 0
      endif
    end function delta

    subroutine sles (m, a, b)
      !
      ! Solve linear equations  A X = B. As in LAPACK  the matrix A is
      ! destroyed and the result is returned in B.
      !
      implicit none
      integer, intent (in) :: m
      real (rk), intent (inout) :: a(m, m), b(m, m)
      ! *** end of interface ***

      integer ipiv(m), info

      ! B will  be overwriten with  the result, A will  be overwritten
      ! with its factorization:
      call dgesv (m, m, a, m, ipiv, b, m, info)
      if (info /= 0) stop "dgesv failed"
    end subroutine sles
  end function oz_equation_c_t_MxM


  elemental function lj (r) result (f)
    !
    ! To be called as in eps * lj (r / sigma)
    !
    implicit none
    real (rk), intent (in) :: r   ! r / sigma, in general
    real (rk) :: f
    ! *** end of interfce ***

    real (rk) :: sr6

    sr6 = 1 / r**6

    f = 4 * sr6 * (sr6 - 1)
  end function lj


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


  function ndst (a) result (b)
    implicit none
    real (rk), intent (in) :: a(:)
    real (rk) :: b(size (a))
    ! *** end of interface ***

    real (rk) :: norm

    norm = 2 * size (b)         ! cast to real
    b = dst (a) / sqrt (norm)
  end function ndst


  subroutine test_dst (n)
    implicit none
    integer, intent (in) :: n
    ! *** end of interface ***

    real (rk) :: a(n)

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
    do i = 1, n
       g(i) = f(i) * (2 * i - 1)
    enddo

    g = dst (g)

    do i = 1, n
       g(i) = g(i) / (2 * i - 1)
    enddo
  end function fourier


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
    f = f / (sum (r**2 * f) * 4 * pi * dr)

    ! Forward transform:
    g = fourier (f) * (dr**3 * (n/pi) / a)

    ! Backward transform:
    h = fourier (g) * (dk**3 * (n/pi) / b) ! a * b = 2pi

    print *, "# norm (f )^2 =", sum ((r * f)**2) * 4 * pi * dr
    print *, "# norm (g )^2 =", sum ((k * g)**2) * 4 * pi * dk
    print *, "# norm (f')^2 =", sum ((r * h)**2) * 4 * pi * dr

    print *, "# |f - f'| =", maxval (abs (f - h))

    print *, "# int (f') =", sum (r**2 * h) * 4 * pi * dr
    print *, "# int (f ) =", sum (r**2 * f) * 4 * pi * dr

    ! This should correspond  to the convolution (f *  f) which should
    ! be again a gaussian, twice as "fat":
    h = g * g
    h = fourier (h) * (dk**3 * (n/pi) / b) ! a*a*b = 1

    print *, "# int (h ) =", sum (r**2 * h) * 4 * pi * dr

    ! Compare width as <r^2>:
    print *, "# sigma (f) =", sum (r**4 * f) * 4 * pi * dr
    print *, "# sigma (h) =", sum (r**4 * h) * 4 * pi * dr

    ! print *, "# n=", n
    ! print *, "# r, f, k, g, h = (f*f)"
    ! do i = 1, n
    !    print *, r(i), f(i), k(i), g(i), h(i)
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
    real (rk), intent (in) :: rho, beta
    ! *** end of interface ***

    real (rk), parameter :: &
         psol_hoef(6) = [0.92302d0, -0.09218d0, +0.62381d0, -0.82672d0, +0.49124d0, -0.10847d0], &
         pliq_hoef(6) = [0.91070d0, -0.25124d0, +0.85861d0, -1.08918d0, +0.63932d0, -0.14433d0], &
         psol_mast(6) = [0.908629d0, -0.041510d0, +0.514632d0, -0.708590d0, +0.428351d0, -0.095229d0], &
         pliq_mast(6) = [0.90735d0, -0.27120d0, +0.91784d0, -1.16270d0, +0.68012d0, -0.15284d0]

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
    real (rk), intent (in) :: p(0:), x
    real (rk) :: y
    ! *** end of interface ***

    integer :: n

    y = 0.0
    do n = 0, size (p) - 1
       y = y + p(n) * x**n
    enddo
  end function poly


  function pad (s) result (t)
    use iso_c_binding, only: c_null_char
    implicit none
    character (len=*), intent (in) :: s
    character (len=len (s)) :: t
    ! *** end of interface ***

    integer :: i

    t = s
    do i = 1, len (t)
       if (t(i:i) /= c_null_char) cycle
       t(i:) = " "
       exit
    enddo
  end function pad
end module rism


program rism_prog
  use rism, only: rism_main, site, problem_data, bgy3d_problem_data
  implicit none

  type (problem_data) :: pd
  type (site), parameter :: a(1) = &
       [site ("lj", [0.0d0, 0.0d0, 0.0d0], 1.0d0, 1.0d0, 0.0d0)]

  pd = bgy3d_problem_data()

  call rism_main (pd, size (a), a)
end program rism_prog

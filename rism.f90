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

  !
  ! These should  obey: ab  = (2π)³, a²b  = (2π)³. The  first equality
  ! guarantees  that that the  forward- and  backward FT  are mutually
  ! inverse,  whereas the  second equation  will make  the convolution
  ! theorem factor-less as in FT(f * g) = FT(f) FT(g):
  !
  real (rk), parameter :: FT_FW = 1           ! == a
  real (rk), parameter :: FT_BW = (2 * pi)**3 ! == b

  !
  ! The interaction energy of two unit charges separated by 1 A is
  !                  -1
  !   E = 1 * 1 / 1 A   = 0.529 au = 332 kcal [/ mol]
  !
  ! The next parameter appears to have the meaning of this interaction
  ! energy of such two unit charges:
  !
  !   EPSILON0INV = 1 / ε₀
  !
  ! and  is  used  to  covert electrostatic  interaction  energies  to
  ! working units.   It has  to be consistent  with other  force field
  ! parameters defined in bgy3d-solvents.h, notably with Lennard-Jones
  ! parameters σ and ε (FIXME: so maybe it was not a good idea to move
  ! it here). These are the original comments:
  !
  !   You have: e^2/4/pi/epsilon0/angstrom, you want: kcal/avogadro/mol
  !
  !   => 331.84164
  !
  ! Again, the code appears to use the IT-calorie to define 1/ε₀. Keep
  ! this in sync with bgy3d.h:
  !
  real (rk), parameter :: EPSILON0INV = 331.84164d0

  !
  ! Inverse range parameter for  separation of the Coulomb into short-
  ! and long range components. The  inverse of this number 1/α has the
  ! dimension of length. FIXME: with the hardwired parameter like this
  ! the  short-range  Coulomb is  effectively  killed alltogether  for
  ! water-like particles  --- a water-like  particle has a  typical LJ
  ! radius σ = 3.16 A. There are no visible changes in the RDF with or
  ! without short range Coulomb and the charges of the order ±1.
  !
  real (rk), parameter :: ALPHA = 1.2

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
  ! context to be passed to the C-world. Upon callback the Fortran sub
  ! implementing a procedure(c_iterator) can use the procedure pointer
  ! to perform the actual work.  I wish closures could be made simpler
  ! than that.
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
    real (rk), dimension (n, size (sites), size (sites)) :: &
         v, vk, t, c, g

    ! Radial grids:
    real (rk) :: r(n), dr, dk
    real (rk) :: k(n)

    integer :: i, m

    m = size (sites)

    ! dr * dk = 2π/2n:
    dr = rmax / n
    dk = pi / rmax
    forall (i = 1:n)
       r(i) = (i - half) * dr
       k(i) = (i - half) * dk
    end forall

    ! Tabulate short-range  pairwise potential  v() on the  r-grid and
    ! long-range pairwise potential vk() on the k-grid:
    call force_field (sites, r, k, v, vk)

    ! Intitial guess:
    t = 0.0

    ! Find t such that iterate_t  (t) == 0. FIXME: passing an internal
    ! function as a callback is an F08 feature. GFortran 4.3 on Debian
    ! Lenny does not support that:
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
            c(:, i, j) = fourier (c(:, i, j)) * (dr**3 / FT_FW)
         enddo
      enddo

      !
      ! The  real-space representation  encodes  only the  short-range
      ! part  of  the  direct   corrlation.  The  (fixed)  long  range
      ! contribution is added here:
      !
      !   C := C  - βV
      !         S     L
      !
      c = c - beta * vk

      !
      ! OZ equation, involves "convolutions", take care of the
      ! normalization here:
      !
      dt = oz_equation_c_t (rho, c)

      !
      ! Since we plugged  in the Fourier transform of  the full direct
      ! correlation including the long range part into the OZ equation
      ! what we get out is the full indirect correlation including the
      ! long-range part.  The menmonic is  C + T is short range.  Take
      ! it out:
      !
      !   T  := T - βV
      !    S          L
      !
      dt = dt - beta * vk

      ! Inverse FT via DST:
      do j = 1, m
         do i = 1, m
            dt(:, i, j) = fourier (dt(:, i, j)) * (dk**3 / FT_BW)
         enddo
      enddo

      ! Return the increment that vanishes at convergence:
      dt = dt - t
    end function iterate_t
  end subroutine rism1d


  subroutine force_field (sites, r, k, vr, vk)
    !
    ! Force field  is represented by two contributions,  the first one
    ! is (hopefully) of short range and is returned in array vr on the
    ! r-grid. The other term is  then of long range and, naturally, is
    ! represented by array vk on the k-grid.
    !
    implicit none
    type (site), intent (in) :: sites(:)   ! (m)
    real (rk), intent (in) :: r(:)         ! (n)
    real (rk), intent (in) :: k(:)         ! (n)
    real (rk), intent (out) :: vr(:, :, :) ! (n, m, m)
    real (rk), intent (out) :: vk(:, :, :) ! (n, m, m)
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

          ! Short range on the r-grid:
          if (sigma /= 0.0) then
             vr(:, i, j) = epsilon * lj (r / sigma) + &
                  EPSILON0INV * charge * coulomb_short (r, ALPHA)
          else
             vr(:, i, j) = &
                  EPSILON0INV * charge * coulomb_short (r, ALPHA)
          endif

          ! Long range on the k-grid:
          vk(:, i, j) = &
               EPSILON0INV * charge * coulomb_long_fourier (k, ALPHA)
       enddo
    enddo

  contains

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
  end subroutine force_field


  !
  ! In the most general case
  !
  !   coulomb_long (r, a) = (1/ε₀) erf (a r) / r
  !
  ! The corresponding Fourier transform is
  !
  !   coulomb_long_fourier (k, a) = (4π/ε₀) exp (- k² / 4a²) / k²
  !
  ! Note that the long-range Coulomb is regular in real space:
  !
  !   erf(x) / x = 2 / √π - 2x² / 3√π + O(x⁴)
  !
  ! So that for a typical value of a = 1.2 the long range potential at r
  ! = 0  is 331.84164 * 1.2 *  1.12837916709551 ~ 449 kcal  for two unit
  ! charges. By doing a finite-size FFT of the k-gitter appriximation of
  ! the  above  Fourier representation  you  will  likely get  something
  ! different.
  !
  elemental function coulomb_long_fourier (k, alpha) result (f)
    implicit none
    real (rk), intent (in) :: k, alpha
    real (rk) :: f
    ! *** end of interface ***

    if (k == 0.0) then
       f = 0.0                  ! 1/k² is undefined
    else
       f = 4 * pi * exp (-k**2 / (4 * alpha**2)) / k**2
    endif
  end function coulomb_long_fourier


  elemental function coulomb_long (r, alpha) result (f)
    !
    ! It is unlikely that you  need this function as as the long-range
    ! interactions are  best specified  on the k-grid  and not  on the
    ! r-grid. See coulomb_long_fourier() instead.
    !
    implicit none
    real (rk), intent (in) :: r, alpha
    real (rk) :: f
    ! *** end of interface ***

    ! ERF() is F08 and later:
    if (r == 0.0) then
       f = alpha * 2 / sqrt (pi)
    else
       f = erf (alpha * r) / r
    endif
  end function coulomb_long


  elemental function coulomb_short (r, alpha) result (f)
    implicit none
    real (rk), intent (in) :: r, alpha
    real (rk) :: f
    ! *** end of interface ***

    ! This will return infinity for r = 0. ERFC() is F08 and later:
    f = erfc (alpha * r) / r
  end function coulomb_short


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
    real (rk), intent (in), target :: x(n)
    real (rk), intent (out), target :: dx(n)
    ! *** end of interface ***

    type (context), pointer :: f_ctx

    ! We  dont  so any  work  ourselves,  just  extract a  pointer  to
    ! procedure(f_iterator) an let it do the rest:
    call c_f_pointer (ctx, f_ctx)

    block
       real (rk), pointer :: y(:, :, :), dy(:, :, :)
       integer :: n(3)

       n = f_ctx % shape

       ! F03 pointer association with contiguous array:
       y(1:n(1), 1:n(2), 1:n(3)) => x
       dy(1:n(1), 1:n(2), 1:n(3)) => dx

       ! The warning by GF 4.6 is incorrect:
       ! http://gcc.gnu.org/bugzilla/show_bug.cgi?id=55855
       dy = f_ctx % f (y)
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

    g = 2 * n * dst (g)

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
    g = fourier (f) * (dr**3 / FT_FW)

    ! Backward transform:
    h = fourier (g) * (dk**3 / FT_BW)

    print *, "# norm (f )^2 =", sum ((r * f)**2) * 4 * pi * dr
    print *, "# norm (g )^2 =", sum ((k * g)**2) * 4 * pi * dk
    print *, "# norm (f')^2 =", sum ((r * h)**2) * 4 * pi * dr

    print *, "# |f - f'| =", maxval (abs (f - h))

    print *, "# int (f') =", sum (r**2 * h) * 4 * pi * dr
    print *, "# int (f ) =", sum (r**2 * f) * 4 * pi * dr

    ! This should correspond  to the convolution (f *  f) which should
    ! be again a gaussian, twice as "fat":
    h = g * g
    h = fourier (h) * (dk**3 / FT_BW)

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

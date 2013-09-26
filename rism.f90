module rism
  use kinds, only: rk
  implicit none
  private

  public :: rism_solvent
  public :: rism_solute
  !
  ! *** END OF INTERFACE ***
  !

  real (rk), parameter :: pi = 4 * atan (1.0_rk)

  !
  ! Working units *are*  angstroms and kcals.  Still when  you wish to
  ! make dimensions  explicit (and  you should) use  combinations like
  ! 0.5 * angstrom, or 1.2 * angstrom**(-1), or 332 kcal * angstrom:
  !
  real (rk), parameter :: ANGSTROM = 1
  real (rk), parameter :: KCAL = 1

  ! ITC calorie here, see also bgy3d.h:
  real (rk), parameter :: KJOULE = KCAL / 4.1868d0 ! only for output

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
  real (rk), parameter :: EPSILON0INV = 331.84164d0 * KCAL * ANGSTROM

  !
  ! Inverse range parameter for  separation of the Coulomb into short-
  ! and long range components. The  inverse of this number 1/α has the
  ! dimension of length. FIXME: with the hardwired parameter like this
  ! the  short-range  Coulomb is  effectively  killed alltogether  for
  ! water-like particles  --- a water-like  particle has a  typical LJ
  ! radius σ = 3.16 A. There are no visible changes in the RDF with or
  ! without short range Coulomb and the charges of the order ±1.
  !

  real (rk), parameter :: ALPHA = 1.2d0 * ANGSTROM**(-1)

  interface gnuplot
     module procedure gnuplot1
     module procedure gnuplot3
  end interface gnuplot

contains

  function rism_nrad (pd) result (nrad) bind (c)
    !
    ! Needs to be consistent with ./rism.h
    !
    use iso_c_binding, only: c_int
    use foreign, only: problem_data
    implicit none
    type (problem_data), intent (in) :: pd ! no VALUE
    integer (c_int) :: nrad
    ! *** end of interface ***

    nrad = maxval (pd % n)
  end function rism_nrad


  function rism_rmax (pd) result (rmax) bind (c)
    !
    ! Needs to be consistent with ./rism.h
    !
    use iso_c_binding, only: c_double
    use foreign, only: problem_data
    implicit none
    type (problem_data), intent (in) :: pd ! no VALUE
    real (c_double) :: rmax
    ! *** end of interface ***

    rmax = maxval (pd % L) / 2
  end function rism_rmax


  function rism_upscale (pd) result (xpd) bind (c)
    !
    ! Needs to be consistent with ./rism.h
    !
    use foreign, only: problem_data
    implicit none
    type (problem_data), intent (in) :: pd ! no VALUE
    type (problem_data) :: xpd
    ! *** end of interface ***

    xpd = pd

    ! FIXME: get rid of redundancies:
    xpd % n = xpd % n * 16
    xpd % L = xpd % L * 4
    xpd % h = xpd % L / xpd % n
  end function rism_upscale


  subroutine rism_solvent (pd, m, solvent, t_buf, x_buf, ptr) bind (c)
    !
    ! When either t_buf or x_buf is not NULL it should be a pointer to
    ! a buffer  with sufficient space  to store nrad  * m *  m doubles
    ! where nrad = maxval (pd % n). If this is the case for t_buf then
    ! the real space represented  of indirect correlation t is written
    ! there.  If x_buf is not  NULL then the Fourier representation of
    ! the solvent susceptibility (not offset by one) is written there.
    !
    ! A  NULL  C-pointer  will   be  cast  by  c_f_pointer()  into  an
    ! unassociated   Fortran  pointer   which  is   isomorphic   to  a
    ! non-present optional argument down the call chain.
    !
    ! Needs to be consistent with ./rism.h
    !
    use iso_c_binding, only: c_int, c_ptr, c_f_pointer
    use foreign, only: problem_data, site, bgy3d_problem_data_print
    use lisp, only: obj
    implicit none
    type (problem_data), intent (in) :: pd ! no VALUE
    integer (c_int), intent (in), value :: m
    type (site), intent (in) :: solvent(m)
    type (c_ptr), intent (in), value :: t_buf ! double[m][m][nrad] or NULL
    type (c_ptr), intent (in), value :: x_buf ! double[m][m][nrad] or NULL
    type (c_ptr), intent (in), value :: ptr   ! SCM* or NULL
    ! *** end of interface ***

    integer :: nrad
    real (rk), pointer :: t(:, :, :), x(:, :, :)
    type (obj), pointer :: dict

    nrad = rism_nrad (pd)

    call c_f_pointer (t_buf, t, shape = [nrad, m, m])
    call c_f_pointer (x_buf, x, shape = [nrad, m, m])
    call c_f_pointer (ptr, dict)

    ! Fill supplied  storage with solvent susceptibility.  If NULL, it
    ! is interpreted as not present() in main() and thus ignored:
    call main (pd, solvent, t=t, x=x, dict=dict)
  end subroutine rism_solvent


  subroutine rism_solute (pd, n, solute, m, solvent, x_buf, ptr) bind (c)
    !
    ! A  NULL  C-pointer  will   be  cast  by  c_f_pointer()  into  an
    ! unassociated   Fortran  pointer   which  is   isomorphic   to  a
    ! non-present optional argument down the call chain.
    !
    ! Needs to be consistent with ./rism.h
    !
    use iso_c_binding, only: c_int, c_ptr, c_f_pointer
    use foreign, only: problem_data, site, bgy3d_problem_data_print
    use lisp, only: obj
    implicit none
    type (problem_data), intent (in) :: pd ! no VALUE
    integer (c_int), intent (in), value :: n, m
    type (site), intent (in) :: solute(n)
    type (site), intent (in) :: solvent(m)
    type (c_ptr), intent (in), value :: x_buf ! double[m][m][nrad] or NULL
    type (c_ptr), intent (in), value :: ptr ! SCM* or NULL
    ! *** end of interface ***

    integer :: nrad
    real (rk), pointer :: x(:, :, :)
    type (obj), pointer :: dict

    nrad = rism_nrad (pd)

    call c_f_pointer (x_buf, x, shape = [nrad, m, m])
    call c_f_pointer (ptr, dict)

    call main (pd, solvent, solute, x=x, dict=dict)
  end subroutine rism_solute


  subroutine main (pd, solvent, solute, t, x, dict)
    !
    ! This one does not need to be interoperable.
    !
    use foreign, only: problem_data, site, bgy3d_problem_data_print, &
         verbosity, comm_rank
    use lisp, only: obj, nil, cons, sym
    implicit none
    type (problem_data), intent (in) :: pd
    type (site), intent (in) :: solvent(:)
    type (site), optional, intent (in) :: solute(:)
    real (rk), optional, intent (out) :: t(:, :, :) ! (nrad, m, m)
    real (rk), optional, intent (out) :: x(:, :, :) ! (nrad, m, m)
    type (obj), intent (out), optional :: dict
    ! *** end of interface ***

    logical :: vv, uv
    integer :: nrad
    real (rk) :: rmax

    ! We do  vv-caclulation if there is  no solute, or  if the solvent
    ! susceptibility was not supplied for solute/solvent calculation.
    uv = present (solute)
    vv = .not. uv .or. uv .and. .not. present (x)

    rmax = rism_rmax (pd)
    nrad = rism_nrad (pd)

    if (verbosity > 0 .and. comm_rank() == 0) then
       print *, "# L =", rmax, "(for 1d)"
       print *, "# N =", nrad, "(for 1d)"

       call bgy3d_problem_data_print (pd)

       call show_sites ("Solvent", solvent)

       if (present (solute)) then
          call show_sites ("Solute", solute)
       endif

       print *, "# ε(RISM) =", epsilon_rism (pd % beta, pd % rho, solvent)
    endif

    ! This is applicable to LJ only, and should take reduced
    ! density/temperature:
    ! call print_info (rho = pd % rho, beta = pd % beta)

    block
       ! Indirect correlation γ aka t = h - c:
       real (rk) :: gam(nrad, size (solvent), size (solvent))

       ! Solvent susceptibility χ = ω + ρh:
       real (rk) :: chi_fft(nrad, size (solvent), size (solvent))

       ! Procedures  down  the  call  chain  will set  this  to  valid
       ! assiciation lists:
       type (obj) :: vdict, udict

       if (vv) then
          call rism_vv (pd % closure, nrad, rmax, pd % beta, pd % rho, &
               solvent, gam, chi_fft, vdict)

          ! Output, if requested:
          if (present (t)) then
             t = gam
          endif

          if (present (x)) then
             x = chi_fft
          endif
       else
          if (.not. present (x)) error stop "no way to get chi!"
          chi_fft = x
          vdict = nil
       endif

       if (uv) then
          call rism_uv (pd % closure, nrad, rmax, pd % beta, pd % rho, &
               solvent, chi_fft, solute, udict)
       endif

       if (present (dict)) then
          dict = cons (cons (sym ("solvent"), vdict), nil)
          if (uv) then
             dict = cons (cons (sym ("solute"), udict), dict)
          endif
       endif
    end block
  end subroutine main


  subroutine rism_vv (method, nrad, rmax, beta, rho, sites, gam, chi, dict)
    use fft, only: fourier_many, FT_FW, FT_BW, integrate
    use snes, only: snes_default
    use foreign, only: site, comm_rank
    use options, only: getopt
    use lisp, only: obj
    implicit none
    integer, intent (in) :: method          ! HNC, KH, or PY
    integer, intent (in) :: nrad            ! grid size
    real (rk), intent (in) :: rmax          ! cell size
    real (rk), intent (in) :: beta          ! inverse temp
    real (rk), intent (in) :: rho
    type (site), intent (in) :: sites(:)    ! (m)
    real (rk), intent (out) :: gam(:, :, :) ! (nrad, m, m)
    real (rk), intent (out) :: chi(:, :, :) ! (nrad, m, m)
    type (obj), intent (out) :: dict
    ! *** end of interface ***


    ! Pair quantities. FIXME: they are symmetric, one should use that:
    real (rk), dimension (nrad, size (sites), size (sites)) :: &
         v, vk, t, c, wk, xk

    ! Radial grids:
    real (rk) :: r(nrad), dr
    real (rk) :: k(nrad), dk

    integer :: i, m

    real (rk) :: eps        ! zero if not supplied in the command line
    real (rk) :: A          ! Cummings & Stell factor [CS81]


    ! FIXME: try  not to  proliferate use of  "environments", function
    ! behaviour is better controlled via its arguments:
    if (.not. getopt ("dielectric", eps)) then
       !
       ! Make sure it is not used anywhere if it was not supplied. The
       ! logic is  if eps  /= 0  then do DRISM,  otherwise ---  do the
       ! usual RISM with whatever dielectric constant comes out of it:
       !
       eps = 0.0
    endif

    ! Set  Cummings &  Stell  semi-empiric scaling  factor  A for  the
    ! long-range  site-site correlation.   This  conflicts with  DRISM
    ! approach.  FIXME: there is  currently no input control to switch
    ! to this method (that is why the hardwired no-op branch choice):
    if (.false.) then
       block
          real (rk) :: y

          y = dipole_density (beta=beta, rho=rho, sites=sites)

          ! FIXME: this factor  is infinite when y ->  0 (e.g. zero
          ! dipole) while keeping the target ε - 1 finite:
          A = dipole_factor (e=eps, y=y)
          if (comm_rank () == 0) then
             print *, "# Scaling the long-range asymptotics by A=", A, "for e=", eps
          endif
       end block
    else
       ! A = 1 has no effect:
       A = 1.0
    endif


    m = size (sites)

    ! dr * dk = 2π/2n:
    dr = rmax / nrad
    dk = pi / rmax
    forall (i = 1:nrad)
       r(i) = (2 * i - 1) * dr / 2
       k(i) = (2 * i - 1) * dk / 2
    end forall

    ! Tabulate short-range  pairwise potentials v() on  the r-grid and
    ! long-range pairwise potential vk() on the k-grid:
    call force_field (sites, sites, r, k, v, vk)

    ! Rigid-bond correlations on the k-grid:
    wk = omega_fourier (sites, k)

    !
    ! FIXME: the  dipole correction x(k) is meaningless  if the dipole
    ! vector is  zero.  Watch for the warning  issued by dipole_axes()
    ! in such cases:
    !
    if (eps /= 0.0) then
       xk = dipole_correction (beta, rho, eps, sites, k)
    else
       xk = 0.0                 ! maybe useful to avoid branches later
    endif

    ! Intitial guess:
    t = 0.0

    ! Find t such that iterate_t  (t) == 0. FIXME: passing an internal
    ! function as a callback is an F08 feature. GFortran 4.3 on Debian
    ! Lenny does not support that:
    call snes_default (iterate_t, t)

    ! Do not assume c has a meaningfull value, it was overwritten with
    ! c(k):
    c = closure (method, beta, v, t)

    !
    ! Return the  indirect correlation t(r)  aka γ(r) in real  rep and
    ! solvent susceptibility χ(k) in Fourier rep:
    !
    gam = t

    block
       real (rk) :: h(nrad, m, m)
       integer :: i, j

       h = fourier_many (c + t) * (dr**3 / FT_FW)

       do i = 1, m
          do j = 1, m
             chi(:, i, j) = wk(:, i, j) + rho * h(:, i, j)
          enddo
       enddo
    end block

    ! Done with it, print results. Here solute == solvent:
    call post_process (method, beta, rho, sites, sites, dr, dk, v, t, &
         A=A, eps=eps, dict=dict)

  contains

    function iterate_t (t) result (dt)
      !
      ! Closure over  host variables: r, k,  dr, dk, v,  c, beta, rho,
      ! ... Implements procedure(f_iterator).
      !
      implicit none
      real (rk), intent (in) :: t(:, :, :) ! (nrad, m, m)
      real (rk) :: dt(size (t, 1), size (t, 2), size (t, 3))
      ! *** end of interface ***

      c = closure (method, beta, v, t)

      ! Forward FT via DST:
      c = fourier_many (c) * (dr**3 / FT_FW)

      !
      ! The  real-space representation  encodes  only the  short-range
      ! part  of  the  direct   corrlation.  The  (fixed)  long  range
      ! contribution is added here:
      !
      !   C := C  - βV
      !         S     L
      !
      ! For dipole solvents Cummings and Stell implied that an overall
      ! scaling  factor A  in  the long  range  expression for  direct
      ! correlation is apropriate [CS81]:
      !
      c = c - (beta * A) * vk

      !
      ! OZ   equation,  involves  convolutions,   take  care   of  the
      ! normalization here --- this is one of a few places where it is
      ! assumed that convolution theorem is factor-less:
      !
      if (eps /= 0.0) then
         ! DRISM. This is going to break for  ρ = 0 as x(k) ~ 1/ρ², it
         ! will also break if dipole density y = 0:
         dt = oz_vv_equation_c_t (rho, c, wk + rho * xk) + xk
      else
         ! The branch is redundand since x(k)  = 0 in this case. It is
         ! here only not to waste CPU time:
         dt = oz_vv_equation_c_t (rho, c, wk)
      endif

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
      ! The factor A is due to Cummings and Stell.
      !
      dt = dt - (beta * A) * vk

      ! Inverse FT via DST:
      dt = fourier_many (dt) * (dk**3 / FT_BW)

      ! Return the increment that vanishes at convergence:
      dt = dt - t
    end function iterate_t
  end subroutine rism_vv


  subroutine rism_uv (method, nrad, rmax, beta, rho, solvent, chi, solute, dict)
    use snes, only: snes_default
    use foreign, only: site
    use fft, only: integrate
    use lisp, only: obj
    use options, only: getopt
    implicit none
    integer, intent (in) :: method         ! HNC, KH, or PY
    integer, intent (in) :: nrad           ! grid size
    real (rk), intent (in) :: rmax         ! cell size
    real (rk), intent (in) :: beta         ! inverse temp
    real (rk), intent (in) :: rho
    type (site), intent (in) :: solvent(:) ! (m)
    real (rk), intent(in) :: chi(:, :, :)  ! (nrad, m, m)
    type (site), intent (in) :: solute(:)  ! (n)
    type (obj), intent (out) :: dict
    ! *** end of interface ***

    ! Solute-solvent pair quantities:
    real (rk), dimension (nrad, size (solute), size (solvent)) :: &
         v, vk, t, c, expB      ! (nrad, n, m)

    ! Solute-solute pair quantities:
    real (rk), dimension (nrad, size (solute), size (solute)) :: &
         wk                     ! (nrad, m, m)

    ! Radial grids:
    real (rk) :: r(nrad), dr
    real (rk) :: k(nrad), dk

    integer :: i, m, n
    logical rbc

    n = size (solute)
    m = size (solvent)

    ! dr * dk = 2π/2n:
    dr = rmax / nrad
    dk = pi / rmax
    forall (i = 1:nrad)
       r(i) = (2 * i - 1) * dr / 2
       k(i) = (2 * i - 1) * dk / 2
    end forall

    ! Tabulate short-range  pairwise potentials v() on  the r-grid and
    ! long-range pairwise potential vk() on the k-grid:
    call force_field (solute, solvent, r, k, v, vk)

    ! Rigid-bond solute-solute correlations on the k-grid:
    wk = omega_fourier (solute, k)

    rbc = getopt ("rbc")

    if (rbc) then
       call bridge (solute, solvent, beta, r, dr, k, dk, expB)
    else
       ! Not used in this case:
       expB = 1.0
    endif


    ! Intitial guess:
    t = 0.0

    ! Find t such that iterate_t  (t) == 0. FIXME: passing an internal
    ! function as a callback is an F08 feature. GFortran 4.3 on Debian
    ! Lenny does not support that:
    call snes_default (iterate_t, t)

    ! Done with it, print results:
    call post_process (method, beta, rho, solvent, solute, dr, dk, v, t, &
         A=1.0d0, eps=0.0d0, dict=dict)

  contains

    function iterate_t (t) result (dt)
      !
      ! Closure over  host variables: r, k,  dr, dk, v,  c, beta, rho,
      ! ... Implements procedure(f_iterator).
      !
      use fft, only: fourier_many, FT_FW, FT_BW
      implicit none
      real (rk), intent (in) :: t(:, :, :) ! (nrad, m, m)
      real (rk) :: dt(size (t, 1), size (t, 2), size (t, 3))
      ! *** end of interface ***

      if (rbc) then
          c = closure_rbc (method, beta, v, t, expB)
      else
          c = closure (method, beta, v, t)
      endif

      ! Forward FT via DST:
      c = fourier_many (c) * (dr**3 / FT_FW)

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
      dt = oz_uv_equation_c_t (c, wk, chi)

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
      dt = fourier_many (dt) * (dk**3 / FT_BW)

      ! Return the increment that vanishes at convergence:
      dt = dt - t
    end function iterate_t
  end subroutine rism_uv


  subroutine bridge (solute, solvent, beta, r, dr, k, dk, expB)
    !
    ! Calculate repulsive bridge correction exp[-B(r)]:
    !
    !
    ! exp[-B  (r)] = Π ω(r)  * exp[-4βε  (σ  / r)¹²]
    !       ij      l≠j    il          lj  lj
    !
    use fft, only: fourier_many, FT_BW, FT_FW
    use foreign, only: site
    implicit none
    type (site), intent (in) :: solute(:)     ! (n)
    type (site), intent (in) :: solvent(:)    ! (m)
    real (rk), intent (in) :: beta            ! inverse temperature
    real (rk), intent (in) :: r(:)            ! (nrad)
    real (rk), intent (in) :: k(:)            ! (nrad)
    real (rk), intent (in) :: dr, dk
    real (rk), intent (out) :: expB(:, :, :) ! (nrad, n, m)
    ! *** end of interface ***

    integer :: m, n

    m = size (solvent)
    n = size (solute)

    block
      ! Solvent-solvent pair quantities:
      real (rk), dimension (size (r), m, m) :: w

      ! Storage for solute-solvent pair quantities (dont take the name
      ! g literally):
      real (rk), dimension (size (r), n, m) :: f, g

      ! Solvent  rigid-bond   correlations  ω  on   the  k-grid.   The
      ! self-correlation  is  a  δ-function   (or  1  in  the  Fourier
      ! representation).
      w = omega_fourier (solvent, k)

      ! Repulsive part of LJ potential goes to f(:, i, j). First index
      ! is that of a solute site  i, second index is of a solvent site
      ! j:
      call lj_repulsive (solute, solvent, r, f)

      ! A  Boltzman  factor  due   to  the  repulsive  branch  of  the
      ! potential:
      f = exp (-beta * f)

      ! Fourier transform of the  Boltzman factor exp(-βv). It is this
      ! factor   which   appears   in   the   convolution   with   the
      ! intra-molecular solvent-solvent site-site correlation ω:
      f = fourier_many (f) * (dr**3 / FT_FW)

      ! Compute expB(:,  i, j) as a  product over all  solvent sites l
      ! except j.  First, set initial value to 1.0:
      expB(:, :, :) = 1.0
      block
        integer :: i, j, l

        ! Loop over solvent sites indexing product terms:
        do l = 1, m
           ! Convolution  in the  Fourier space  is a  product  of two
           ! Fourier representations:
           do j = 1, m          ! solvent sites
              do i = 1, n       ! solute sites
                 g(:, i, j) =  f(:, i, l) * w(:, l, j)
              enddo
           enddo

           ! Transform convolutions to the real space:
           g = fourier_many (g) * (dk**3 / FT_BW)

           ! Here the  product is accumulated.  FIXME:  Note that even
           ! though the factors  with l == j are  computed above, they
           ! are excluded from the product:
           do j = 1, m          ! solvent sites
              do i = 1, n       ! solute sites
                 if (j == l) cycle
                 expB(:, i, j) = expB(:, i, j) * g(:, i, j)
              enddo
           enddo
        enddo
      end block
    end block
  end subroutine bridge


  function chempot_bridge (beta, rho, h, expB, r, dr) result (mu)
    !
    ! Computes   correction  to   the  chemical   potential,   βμ,  by
    ! integrating the  expression given by  thermodynamic perturbation
    ! theory over the volume:
    !
    !   βμ = 4πρ ∫ [exp(B) - 1] g(r) r²dr
    !
    ! Here dr == 1, scale the result by dr³ if that is not the case.
    !
    use fft, only: integrate
    implicit none
    real (rk), intent (in) :: beta, rho
    real (rk), intent (in) :: h(:, :, :)    ! (nrad, n, m)
    real (rk), intent (in) :: expB(:, :, :) ! (nrad, n, m)
    real (rk), intent (in) :: r(:)          ! (nrad)
    real (rk), intent (in) :: dr
    real (rk) :: mu
    ! *** end of interface ***

    integer :: i, j
    real (rk) :: density (size (r))

    ! According   to  thermodynamic   perturbation   theory,  chemical
    ! potential density to be integrated:
    density = 0.0
    do j = 1, size (h, 3)
       do i = 1, size (h, 2)
          density = (expB(:, i, j) - 1) * (1 + h(:, i, j))
       enddo
    enddo

    ! The integration  procedure assumes a  very specific radial  (i +
    ! 1/2) grid.  Multiply that by dr³ and divide by β to get the real
    ! number:
    mu = integrate (density) * (rho * dr**3 / beta)
  end function chempot_bridge


  subroutine lj_repulsive (asites, bsites, r, vr)
    !
    ! only calculate the repulsive part of LJ potential (r^-12 term)
    !
    use foreign, only: site
    implicit none
    type (site), intent (in) :: asites(:)       ! (n)
    type (site), intent (in) :: bsites(:)       ! (m)
    real (rk), intent (in) :: r(:)              ! (nrad)
    real (rk), intent (out) :: vr(:, :, :)      ! (nrad, n, m)
    ! ** end of interface ***

    real (rk) :: sigma, epsilon
    integer :: i, j
    type (site) :: a, b

    do j = 1, size (bsites)
       b = bsites(j)
       do i = 1, size (asites)
          a = asites(i)

          epsilon = sqrt (a % epsilon * b % epsilon)
          sigma = (a % sigma + b % sigma) / 2

          if (sigma /= 0.0) then
              vr (:, i, j) = epsilon * lj12 (r / sigma)
          else
              vr (:, i, j) = 0.0
          endif
        enddo
    enddo

  contains
    elemental function lj12 (r) result (f)
      !
      ! To be called as in eps * lj12 (r / sigma)
      !
      implicit none
      real (rk), intent (in) :: r   ! r / sigma, in general
      real (rk) :: f
      ! *** end of interfce ***

      real (rk) :: sr12

      sr12 = 1 / r**12

      f = 4 * sr12
    end function lj12
  end subroutine lj_repulsive


  subroutine post_process (method, beta, rho, solvent, solute, dr, dk, v, t, &
       A, eps, dict)
    !
    ! Prints some results.
    !
    use fft, only: fourier_many, FT_FW
    use linalg, only: polyfit
    use foreign, only: site, verbosity, comm_rank, &
         HNC => CLOSURE_HNC, KH => CLOSURE_KH, PY => CLOSURE_PY
    use lisp, only: obj, cons, nil, sym, num
    implicit none
    integer, intent (in) :: method         ! HNC, KH or PY
    real (rk), intent (in) :: beta         ! inverse temperature
    real (rk), intent (in) :: rho
    type (site), intent (in) :: solvent(:) ! (m)
    type (site), intent (in) :: solute(:)  ! (n)
    real (rk), intent (in) :: dr, dk       ! grid steps
    real (rk), intent (in) :: v(:, :, :)   ! (nrad, n, m)
    real (rk), intent (in) :: t(:, :, :)   ! (nrad, n, m)
    real (rk), intent (in) :: A            ! long-range scaling factor
    real (rk), value :: eps     ! requested dielectric constant, or 0
    type (obj), intent (out) :: dict
    ! *** end of interface ***

    integer :: nrad, n, m
    integer :: verb

    if (comm_rank () == 0) then
       verb = verbosity
    else
       verb = 0
    endif

    ! FIXME: to make debug  output somewhat meaningfull for plain RISM
    ! (not DRISM) calculations use a real eps:
    if (eps == 0.0) eps = 78.4d0 ! has value attribute

    nrad = size (t, 1)          ! number of radial point
    n = size (t, 2)             ! number of solute sites
    m = size (t, 3)             ! number of solvent sites

    block
       integer :: j
       integer, parameter :: one(3, 3) = &
            reshape ([1, 0, 0,  0, 1, 0,  0, 0, 1], [3, 3])
       real (rk) :: u(3, 3), x(3, size (solvent))
       type (site) :: sites(size (solvent))

       if (verb > 0) then
          sites = solvent
          call show_sites ("Before!", sites)
          print *, "# Dipole before =", dipole (sites)

          print *, "# Solvent axes:"
          print *, "# center =", center (sites)
          u = dipole_axes (sites)

          print *, "# diff =", maxval (abs (matmul (u, transpose (u)) - one))
          print *, "# diff =", maxval (abs (matmul (transpose (u), u) - one))

          do j = 1, 3
             print *, "# axis(", j, ") =", u(:, j)
          enddo

          x = local_coords (sites, u)
          do j = 1, size (sites)
             sites(j) % x = x(:, j)
          enddo
          call show_sites ("After!", sites)
          print *, "# Dipole after =", dipole (sites)
       endif
    end block

    block
       integer :: p, i, j
       real (rk) :: r(nrad), k(nrad)
       real (rk) :: c(nrad, n, m)
       real (rk) :: cl(nrad, n, m)
       real (rk) :: h(nrad, n, m)
       real (rk) :: g(nrad, n, m)
       real (rk) :: x(nrad), x0(nrad), x1(nrad), x2(nrad), xx(nrad) ! qhq in k-space
       real (rk) :: xd(nrad, m, m)

       ! Dont like to pass redundant  info, recompute r(:) from dr and
       ! k(:) from dk:
       forall (i = 1:nrad)
          r(i) = (2 * i - 1) * dr / 2
          k(i) = (2 * i - 1) * dk / 2
       end forall

       ! For the same reason recomute c and g:
       c = closure (method, beta, v, t)
       h = c + t
       g = 1 + h

       ! Real-space rep of the  long-range correlation. Note the extra
       ! scaling factor A:
       forall (p = 1:nrad, i = 1:n, j = 1:m)
          cl(p, i, j) = - (beta * A) * solute(i) % charge * solvent(j) % charge &
               * EPSILON0INV * coulomb_long (r(p), ALPHA)
       end forall

       if (verb > 0) then
          print *, "# rho =", rho, "beta =", beta, "n =", nrad
       endif

       ! Chemical potentials:
       block
          ! FIXME: as implemented, for method == PY the default energy
          ! functional is GF:
          integer :: methods(3) = [HNC, KH, PY]
          character(len=3) :: names(size (methods)) = ["HNC", "KH ", "GF "]
          real (rk) :: mu(size (methods))
          integer :: i

          ! Initialize intent (out) argument:
          dict = nil
          do i = 1, size (methods)
             mu(i) = chempot (methods(i), rho, h, c, cl) * (dr**3 / beta)

             ! Cons a key/value pair onto the list:
             dict = cons (cons (sym (trim (names(i))), num (mu(i))), dict)
          enddo

          if (verb > 0) then
             print *, "# Chemical potentials, default is marked with *:"
             do i = 1, size (methods)
                print *, "# MU =", mu(i) / KCAL, "kcal =", mu(i) / KJOULE, "kJ", &
                     " (", names(i), ")", merge ("*", " ", method == methods(i))
             enddo
          endif
       end block

       block
          integer :: i, j, p
          real (rk) :: q(size (solvent))
          real (rk) :: y, fac0, fac1, fac2
          real (rk) :: qhq(nrad)
          real (rk) :: hk(nrad, n, m)

          ! Small-k  behavior   of  qh(k)q  which   is  essential  for
          ! dielectric permittivity:
          hk = fourier_many (h) * (dr**3 / FT_FW)

          ! FIXME: dipole_density() assumes all solvent sites have the
          ! same number density:
          y = dipole_density (beta, rho, solvent)

          if (verb > 0) then
             print *, "# y = ", y, "e = 1 + 3y =", 1 + 3 * y
             print *, "# ε = ", epsilon_rism (beta, rho, solvent)
             print *, "# A(1 + 3 * y, y) = ", dipole_factor (1 + 3 * y, y)
             print *, "# A(ε, y) = ", dipole_factor (eps, y)
          endif

          !
          ! The coefficient in the  low-k expansion of the "dielectric
          ! susceptibility".    We   need    a    better   name    for
          ! that.  Ref.  [PP92b]  refers  to it  as  "weighted  second
          ! moment":
          !
          !                            k²
          !   Σ   z  ρ h  (k) ρ z  ~  --- [(ε - 1)/ε - 3y]
          !    ij  i    ij       j    4πβ
          !
          ! Cp.  to  e.g. Eq.  (35)  of Ref.  [HS76]  or  Eq. (15)  of
          ! Ref. [PP92b].
          !
          ! [HS76] Dielectric constant in terms of atom–atom
          !   correlation functions, Johan S.  Høye and G.  Stell, J.
          !   Chem.  Phys.  65, 18 (1976);
          !   http://dx.doi.org/10.1063/1.432793
          !
          ! [PP92b] A site–site theory for finite concentration saline
          !   solutions, John Perkyns and B. Montgomery Pettitt,
          !   J. Chem. Phys. 97, 7656 (1992);
          !   http://dx.doi.org/10.1063/1.463485
          !
          fac0 = ((eps - 1) / eps - 3 * y) / (4 * pi * beta)

          !
          ! When ε  = 1 + 3y  as predicted by  RISM, the corresponding
          ! expression for the factor becomes
          !
          !   fac1 = (-9 * y**2 / (1 + 3 * y)) / (4 * pi * beta)
          !
          ! Though we compute it in the same way:
          !
          associate (eps => 1 + 3 * y)
            fac1 = ((eps - 1) / eps - 3 * y) / (4 * pi * beta)
          end associate

          associate (eps0 => 1 + 3 * y)
            fac2 = (eps - eps0) / (4 * pi * beta)
          end associate

          if (verb > 0) then
             print *, "# [(ε - 1) / ε - 3y] / 4πβ = ", fac0, &
                  "e =", epsln (beta, y, fac0), &
                  "(target values)"
             print *, "#    -9y² / (1 + 3y) / 4πβ = ", fac1, &
                  "e =", epsln (beta, y, fac1), &
                  "(RISM limit)"
             print *, "# [ε  -  (1  +  3y)] / 4πβ = ", fac2, &
                  "(DRISM correction scale)"
          endif

          ! q = ρ * z:
          q(:) = rho * solvent(:) % charge

          xd = dipole_correction (beta, rho, eps, solvent, k)

          ! Note  the  extra  EPSILON0INV   factor,  it  seems  to  be
          ! necessary to get the kcal/A³ units right:
          qhq = 0.0
          x = 0.0
          xx = 0.0
          do i = 1, size (solvent)
             do j = 1, size (solvent)
                do p = 1, nrad
                   qhq(p) = qhq(p) + q(i) * h(p, i, j) * q(j) * EPSILON0INV
                   x(p) = x(p) + q(i) * hk(p, i, j) * q(j) * EPSILON0INV
                   xx(p) = xx(p) + q(i) * xd(p, i, j) * q(j) * EPSILON0INV
                enddo
             enddo
          enddo
          ! where (k < 2.0) x = x - fac * k**2
          x0 = fac0 * k**2
          x1 = fac1 * k**2
          x2 = fac2 * k**2

          block
             real (rk) :: a(0:2), c(0:2)
             integer :: i, n, npts(3) = [3, 5, 15]

             do n = 0, 2
                c(n) = moment (n, qhq, dr)
             enddo

             if (verb > 0) then
                print *, "# [(e - 1) / e - 3y] / 4πβ = ", -c(2) / 6, &
                     "e =", epsln (beta, y, -c(2) / 6), &
                     "(from second moment)"
                if (verb > 1) then
                   print *, "# first three moments = ", c(:)
                endif
             endif

             do i = 1, size (npts)
                n = npts(i)
                if (nrad < n) cycle

                a = polyfit (k(1:n), x(1:n), 2)

                if (verb > 0) then
                   print *, "# [(e - 1) / e - 3y] / 4πβ = ", a(2) , &
                        "e =", epsln (beta, y, a(2)), &
                        "(", n, "pts)"
                   if (verb > 1) then
                      print *, "# fitted coefficients = ", a(:)
                      print *, "# fitted k-range = ", k(1), "...", k(n), &
                           "with", n, "pts"
                   endif
                endif
             enddo
          end block
       end block


       ! This prints a lot of data on tty!
       if (verb > 1) then
          print *, "# r(i) then g(i), v(i), t(i), c(i), each for",  n, "x", m, "pairs"
          do p = 1, nrad
             write (*, *) r(p), &
                  &     ((g(p, i, j), i=1,n), j=1,m), &
                  &     ((v(p, i, j), i=1,n), j=1,m), &
                  &     ((t(p, i, j), i=1,n), j=1,m), &
                  &     ((c(p, i, j), i=1,n), j=1,m), &
                  &      k(p), x(p), x0(p), x1(p), x2(p), xx(p)
          enddo
       endif
    end block

  contains

    pure function epsln (beta, y, c) result (e)
      !
      ! Solves for e in
      !
      !   c = [(e - 1) / e - 3y] / 4πβ
      !
      implicit none
      real (rk), intent (in) :: beta, y, c
      real (rk) :: e
      ! *** end of interface ***

      associate (x => 4 * pi * beta * c + 3 * y)
        e = 1 / (1 - x)         ! x = (e - 1) / e
      end associate
    end function epsln


    pure function moment (n, f, dr) result (m)
      use fft, only: integrate
      implicit none
      integer, intent (in) :: n
      real (rk), intent (in) :: f(:) ! f(r) on the r-grid
      real (rk), intent (in) :: dr   ! grid spacing
      real (rk) :: m
      ! *** end of interface ***

      integer :: i
      real (rk) :: r(size (f))

      ! Recompute r(:) from dr:
      forall (i = 1:size (f))
         r(i) = (2 * i - 1) * dr / 2
      end forall

      m = integrate (r**n * f) * dr**3
    end function moment

  end subroutine post_process


  function omega_fourier (sites, k) result (wk)
    !
    ! The  intra-molecular force-field  amounts to  rigid  bonds. This
    ! function computes the Fourier representation
    !
    !   ω(k) = sinc (kR)
    !
    ! of the corresponding δ-like correlation functions:
    !
    !   ω(r) = δ(r - R) / 4πR²
    !
    use foreign, only: site
    use bessel, only: j0
    implicit none
    type (site), intent (in) :: sites(:)                  ! (m)
    real (rk), intent (in) :: k(:)                        ! (nrad)
    real (rk) :: wk(size (k), size (sites), size (sites)) ! (nrad, m, m)
    ! *** end of inteface ***

    real (rk) :: xa(3), xb(3), rab
    integer :: i, j, m

    m = size (sites)

    ! Site-site bond lengths are derived from site coordinates:
    do j = 1, m
       xb = sites(j) % x
       do i = 1, m
          xa = sites(i) % x

          ! Site-site distance vanishes (at least) for i == j, so that
          ! sinc(kr) == 1 in this case. NORM2() is an F08 intrinsic.
          rab = norm2 (xa - xb)

          ! Rigid bond correlation on a k-grid:
          wk(:, i, j) = j0 (k * rab)
       enddo
    enddo
  end function omega_fourier


  subroutine force_field (asites, bsites, r, k, vr, vk)
    !
    ! Force field  is represented by two contributions,  the first one
    ! is (hopefully) of short range and is returned in array vr on the
    ! r-grid. The other term is  then of long range and, naturally, is
    ! represented by array vk on the k-grid.
    !
    use foreign, only: site
    implicit none
    type (site), intent (in) :: asites(:)  ! (n)
    type (site), intent (in) :: bsites(:)  ! (m)
    real (rk), intent (in) :: r(:)         ! (nrad)
    real (rk), intent (in) :: k(:)         ! (nrad)
    real (rk), intent (out) :: vr(:, :, :) ! (nrad, n, m)
    real (rk), intent (out) :: vk(:, :, :) ! (nrad, n, m)
    ! *** end of inteface ***

    real (rk) :: epsilon, sigma, charge
    integer :: i, j
    type (site) :: a, b

    do j = 1, size (bsites)
       b = bsites(j)
       do i = 1, size (asites)
          a = asites(i)

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


  elemental function closure (method, beta, v, t) result (c)
    use foreign, only: HNC => CLOSURE_HNC, KH => CLOSURE_KH, PY => CLOSURE_PY
    implicit none
    integer, intent (in) :: method
    real (rk), intent (in) :: beta, v, t
    real (rk) :: c
    ! *** end of interface ***

    select case (method)
    case (HNC)
       c = closure_hnc (beta, v, t)
    case (KH)
       c = closure_kh (beta, v, t)
    case (PY)
       c = closure_py (beta, v, t)
    case default
       c = huge (c)            ! FIXME: cannot abort in pure functions
    end select
  end function closure


  elemental function expm1 (x) result (f)
    !
    ! Elemental version of libc expm1().
    !
    use foreign, only: c_expm1 => expm1
    implicit none
    real (rk), intent (in) :: x
    real (rk) :: f
    ! *** end of interface ***

    f = c_expm1 (x)
  end function expm1


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

    ! c = exp (-beta * v + t) - 1 - t
    c = expm1 (-beta * v + t) - t
  end function closure_hnc


  !
  ! For x <= 0 the same as exp(x) - 1, but does not grow exponentially
  ! for positive x:
  !
  elemental function lexpm1 (x) result (y)
    implicit none
    real (rk), intent (in) :: x
    real (rk) :: y
    ! *** end of interface ***

    if (x <= 0.0) then
       y = expm1 (x)
    else
       y = x
    endif
  end function lexpm1


  !
  ! 2) Kovalenko-Hirata (KH) closure.
  !
  elemental function closure_kh (beta, v, t) result (c)
    implicit none
    real (rk), intent (in) :: beta, v, t
    real (rk) :: c
    ! *** end of interface ***

    ! Note that lexpm1() /= expm1():
    c = lexpm1 (-beta * v + t) - t
  end function closure_kh


  ! 3)  Percus-Yevick  (PY)   closure  relation  between  direct-  and
  ! indirect correlation c and γ:
  !
  !   c := exp (-βv) [1 + γ] - 1 - γ
  elemental function closure_py (beta, v, t) result (c)
    implicit none
    real (rk), intent (in) :: beta, v, t
    real (rk) :: c
    ! *** end of interface ***

    c = exp (-beta * v) * (1 + t) - 1 - t
  end function closure_py

  elemental function closure_rbc (method, beta, v, t, expB) result (c)
    use foreign, only: HNC => CLOSURE_HNC, KH => CLOSURE_KH, PY => CLOSURE_PY
    implicit none
    integer, intent (in) :: method
    real (rk), intent (in) :: beta, v, t, expB
    real (rk) :: c
    ! *** end of interface ***

    select case (method)
    ! RBC only with HNC now
    case (HNC)
       c = closure_hnc_rbc (beta, v, t, expB)
    case default
       c = huge (c)            ! FIXME: cannot abort in pure functions
    end select
  end function closure_rbc

  !
  ! HNC closure with repulsive bridge correction:
  !
  !    c := exp (-βv + γ + B) - 1 - γ
  !
  elemental function closure_hnc_rbc (beta, v, t, expB) result (c)
    implicit none
    real (rk), intent (in) :: beta, v, t, expB
    real (rk) :: c
    ! *** end of interface ***

    c = exp (-beta * v + t) * expB - 1 - t
  end function closure_hnc_rbc

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
  function oz_vv_equation_c_t (rho, C, W) result (T)
    implicit none
    real (rk), intent (in) :: rho
    real (rk), intent (in) :: C(:, :, :)         ! (nrad, m, m)
    real (rk), intent (in) :: W(:, :, :)         ! (nrad, m, m)
    real (rk) :: T(size (C, 1), size (C, 2), size (C, 3))
    ! *** end of interface ***

    integer :: i

    ! There is  no reason to  handle the 1x1 case  differently, except
    ! clarity.  The MxM branch should be able to handle that case too.
    if (size (C, 3) == 1) then
       ! FIXME: it  is implied here  that W =  1. See comments  on the
       ! value of ω(k) for i == j in omega_fourier().
       do i = 1, size (C, 1)
          T(i, 1, 1) = oz_vv_equation_c_t_1x1 (rho, C(i, 1, 1))
       enddo
    else
       do i = 1, size (C, 1)
          T(i, :, :) = oz_vv_equation_c_t_MxM (rho, C(i, :, :), W(i, :, :))
       enddo
    endif
  end function oz_vv_equation_c_t


  elemental function oz_vv_equation_c_t_1x1 (rho, c) result (t)
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
  end function oz_vv_equation_c_t_1x1


  function oz_vv_equation_c_t_MxM (rho, C, W) result (T)
    !
    ! So far  rho is the same for  all sites, we may  have mixed left-
    ! and right-side multiplies.
    !
    ! RISM equation, here for h and c:
    !
    !   h = ω * c * ω + ω * c * ρ * h
    !
    ! or
    !                      -1
    !   h = (1 - ω * c * ρ)   * ω * c * ω
    !
    !                              -1
    !     = ω * c * (1 - ρ * ω * c)   * ω
    !
    !                                  -1
    !     = ω * c * ω * (1 - ρ * c * ω)
    !
    ! All three expressions are equivalent --- this can be verified by
    ! a series expansion of the inverse, given that ρ * ω = ω * ρ.  To
    ! use  this  code for  mixtures  of  otherwise uncorrelated  sites
    ! supply  unity matrix  for  ω. A  mixture  of molecular  solvents
    ! should  be represented  by  an ω-matrix  with  a suitable  block
    ! structure.
    !
    use linalg, only: sles
    implicit none
    real (rk), intent (in) :: rho
    real (rk), intent (in) :: C(:, :)      ! (m, m)
    real (rk), intent (in) :: W(:, :)      ! (m, m)
    real (rk) :: T(size (C, 1), size (C, 1)) ! (m, m)
    ! *** end of interface ***

    real (rk), dimension (size (C, 1), size (C, 1)) :: H
    integer :: i, j, m

    m = size (C, 1)

    ! H := WC, temporarily:
    H = matmul (W, C)

    ! T := 1  - WCρ. The output matrix  T is used here as  a free work
    ! array:
    forall (i = 1:m, j = 1:m)
       T(i, j) = delta (i, j) - H(i, j) * rho
    end forall

    ! H := WCW.  Still temporarily --- it will be overwritten with the
    ! real H after solving the linear equations.
    H = matmul (H, W)

    !
    ! Solving the linear equation makes  H have the literal meaning of
    ! the total correlation matrix (input T is destroyed):
    !
    !       -1                 -1
    ! H := T   * H == (1 - WCρ)   * WCW
    !
    call sles (m, T, H)         ! FIXME: copy in/out.

    !
    ! The  same  effect  is  achieved  in  1x1  version  of  the  code
    ! differently:
    !
    ! T := H - C
    !
    T = H - C

  contains

    pure function delta (i, j) result (d)
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

  end function oz_vv_equation_c_t_MxM


  function oz_uv_equation_c_t (cuv, wuu, xvv) result (tuv)
    !
    ! RISM equation, here for h and c:
    !
    !   h   =  ω  * c   * [ω + ρ h  ]
    !    uv     u    uv     v     vv
    !
    ! The term  is square brackets  is the solvent property  (fixed by
    ! the assumption  of the infinite  dilution) and is passed  as the
    ! solvent susceptibility
    !
    !   χ   =  ω + ρ h
    !    vv     v     vv
    !
    ! The returned value is the indirect correlation
    !
    !   t   =  h  -  c
    !    uv     uv    uv
    !
    implicit none
    real (rk), intent (in) :: cuv(:, :, :)         ! (nrad, n, m)
    real (rk), intent (in) :: wuu(:, :, :)         ! (nrad, n, n)
    real (rk), intent (in) :: xvv(:, :, :)         ! (nrad, m, m)
    real (rk) :: tuv(size (cuv, 1), size (cuv, 2), size (cuv, 3))
    ! *** end of interface ***

    integer :: p

    ! Many associative matrix multiplies: NxN * (NxM * MxM) - NxM
    do p = 1, size (cuv, 1)
       tuv(p, :, :) = matmul (wuu(p, :, :), matmul (cuv(p, :, :), xvv(p, :, :))) - cuv(p, :, :)
    enddo
  end function oz_uv_equation_c_t


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

  contains

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
  end subroutine print_info

  function chempot_density (method, rho, h, cs, cl) result (mu)
    !
    ! Returns the β-scaled "density" of the chemical potential, βμ(r).
    ! To  get the  excess  chemical potential  integrate  it over  the
    ! volume and divide by β, cf.:
    !
    !   βμ = 4πρ ∫ [½h²(r) - c(r) - ½h(r)c(r)] r²dr
    !
    ! The first h²-term in the definition of the chemical potential is
    ! of  short  range  and  should  not present  any  difficulty  for
    ! integration.
    !
    ! Assuming that the long-range component of the direct correlation
    ! c₁₂(r)  of two sites,  1 and  2, is  proportional to  the charge
    ! product q₁ and  q₂ as in -βq₁q₂v(r) one may  deduce that ρ₁c₁₁ +
    ! ρ₂c₂₁ vanishes identically for  *neutral* systems.  In this case
    ! the Coulomb tails of the direct correlation in the second c-term
    ! do not  contribute to  the chemical potential.   This is  a good
    ! thing, because otherwise an integral 4π∫r²dr/r would diverge.
    !
    ! The third  hc-term is of the  short range due to  h. However the
    ! long-range  Coulomb-type correction  to the  direct correlation,
    ! since weighted  by a pair  specific total correlation  h₁₂, does
    ! not vanish in such a sum.
    !
    ! To get an idea about the order of magnitude of such a long-range
    ! correction: for  SPC/E water  with σ(H) =  1.0 A, ε(H)  = 0.0545
    ! kcal  [1],  and using  the  short-range  correlation with  range
    ! separation  parameter α  = 1.2  A^-1 gives  the  excess chemical
    ! potential μ = 6.86 kcal (wrong because positive!).  By including
    ! the long  range correlation into  the hc-term one  obtains -4.19
    ! kcal instead.  Here N = 4096, L  = 80 A, ρ = 0.033427745 A^-3, β
    ! = 1.6889 kcal^-1. (You may  need to use --snes-solver "trial" to
    ! get SPC/E water converge).
    !
    ! Note  that,   contrary  to  what  is  being   suggested  in  the
    ! literature, these  numbers depend strongly on  the LJ parameters
    ! of hydrogens. With σ(H) = 0.8  A and ε(H) = 0.046 kcal [Leipzig]
    ! one gets  μ = -4.97  kcal and  with σ(H) =  1.1656 A and  ε(H) =
    ! 0.01553  kcal [Amber] one  gets μ  = -3.66  kcal.  At  least one
    ! source  quotes  yet a  different  number  -3.41  kcal [1].   The
    ! corresponding result for modified TIP3P water with σ(H) = 0.4 A,
    ! and ε(H) = 0.046 kcal is -6.35 kcal.
    !
    ! FIXME: what do we do for charged systems?
    !
    ! [1] "Comparative  Study on Solvation Free  Energy Expressions in
    !     Reference Interaction Site  Model Integral Equation Theory",
    !     Kazuto   Sato,  Hiroshi   Chuman,   and  Seiichiro   Ten-no,
    !     J.  Phys.   Chem.  B,   2005,  109  (36),   pp  17290–17295,
    !     http://dx.doi.org/10.1021/jp053259i
    !
    use foreign, only: HNC => CLOSURE_HNC, KH => CLOSURE_KH
    implicit none
    integer, intent (in) :: method        ! HNC, KH, or anything else
    real (rk), intent (in) :: rho
    real (rk), intent (in) :: h(:, :, :)  ! (nrad, n, m)
    real (rk), intent (in) :: cs(:, :, :) ! (nrad, n, m)
    real (rk), intent (in) :: cl(:, :, :) ! (nrad, n, m)
    real (rk) :: mu(size (h, 1))
    ! *** end of interface ***

    integer :: p, i, j
    real (rk) :: muH, muS, muL, thresh

    select case (method)
    case (KH)
       ! The h² term contributes only in the depletion regions (KH):
       thresh = 0.0
    case (HNC)
       ! The h² term contributes unconditionally (HNC):
       thresh = - huge (thresh)
    case default
       ! There is no h² term otherwise (GF):
       thresh = huge (thresh)
    end select

    do p = 1, size (h, 1)       ! nrad
       muH = 0.0
       muS = 0.0
       muL = 0.0
       do j = 1, size (h, 3)    ! m
          do i = 1, size (h, 2) ! n
             ! The h² term contributes conditionally. Eventually, only
             ! depletion regions  (h < 0)  contribute (KH).  Threshold
             ! is  supposed to  be  0.0 for  KH functional  (depletion
             ! regions contribute),  anywhere between 1 and  +∞ for GF
             ! functional  (no such  term) and  -∞ for  HNC functional
             ! (contributes unconditionally):
             if (-h(p, i, j) > thresh) then
                muH = muH + rho * h(p, i, j)**2 / 2
             endif

             muS = muS + rho * (-cs(p, i, j) - h(p, i, j) * cs(p, i, j) / 2)
             muL = muL + rho * (             - h(p, i, j) * cl(p, i, j) / 2)
          enddo
       enddo

       mu(p) = muH + muS + muL
    enddo
  end function chempot_density

  function chempot (method, rho, h, cs, cl) result (mu)
    !
    ! Computes  the chemical  potential, βμ,  by integration  over the
    ! volume:
    !
    !   βμ = 4πρ ∫ [½h²(r) - c(r) - ½h(r)c(r)] r²dr
    !
    ! Here dr == 1, scale the result by dr³ if that is not the case.
    !
    use fft, only: integrate
    implicit none
    integer, intent (in) :: method        ! HNC, KH, or anything else
    real (rk), intent (in) :: rho
    real (rk), intent (in) :: h(:, :, :)  ! (nrad, n, m)
    real (rk), intent (in) :: cs(:, :, :) ! (nrad, n, m)
    real (rk), intent (in) :: cl(:, :, :) ! (nrad, n, m)
    real (rk) :: mu
    ! *** end of interface ***

    real (rk) :: density (size (h, 1))

    ! Chemical potential density to be integrated:
    density = chempot_density (method, rho, h, cs, cl)

    ! Multiply that by dr³ and divide by β to get the real number:
    mu = integrate (density)
  end function chempot


  subroutine show_sites (name, sites)
    use foreign, only: site, pad
    implicit none
    character (len=*), intent (in) :: name
    type (site), intent (in) :: sites(:)
    ! *** end of interface ***

    integer :: i

    print *, "# ", name

    do i = 1, size (sites)

       write (*, 100) i, pad (sites(i) % name), &
            sites(i) % x, &
            sites(i) % sigma, sites(i) % epsilon, sites(i) % charge
    enddo
100 format (' #', I2, 1X, A5, 6F12.4)
  end subroutine show_sites


  function dipole_density (beta, rho, sites) result (y)
    !
    ! Computes dipole density
    !
    !       4 π    2
    !   y = --- βρμ
    !        9
    !
    use foreign, only: site
    implicit none
    real (rk), intent (in) :: beta, rho
    type (site), intent (in) :: sites(:)
    real (rk) :: y
    ! *** end of interface ***

    !
    ! Here ρμ² has  the dimension of electrostatic energy  which is by
    ! convention    converted   to    kcal    by   multiplying    with
    ! EPSILON0INV.  The inverse  temperature is  supposed to  have the
    ! dimension of 1/kcal here:
    !
    y = 4 * pi * beta * EPSILON0INV * rho * norm2 (dipole (sites))**2 / 9
  end function dipole_density


  function epsilon_rism (beta, rho, sites) result (e)
    !
    ! Computes dielectric constant as  should be predicted by the RISM
    ! theory:
    !
    !   ε = 1 + 3y
    !
    ! where y is the usual dipole density
    !
    !       4 π    2
    !   y = --- βρμ
    !        9
    !
    use foreign, only: site
    implicit none
    real (rk), intent (in) :: beta, rho
    type (site), intent (in) :: sites(:)
    real (rk) :: e
    ! *** end of interface ***

    real (rk) :: y

    y = dipole_density (beta, rho, sites)

    e = 1 + 3 * y
  end function epsilon_rism


  function dipole_factor (e, y) result (a)
    !
    ! Cummings and Stell [CS81] argued that the long-range asymptotics
    ! of the site-site correlation function for a rigid dipolar liquid
    ! differs from the "intuitive" form -βqq/r by a scalar factor A:
    !
    !       1 + ε(3y - 1)
    !   A = -------------
    !         3y(ε - 1)
    !
    ! where ε  and y  have their usual  meaning.  Conversely, if  A is
    ! known or otherwise enforced then the corresponding ε is given by
    !
    !          1 + 3Ay
    !   ε = -------------
    !       1 + 3(A - 1)y
    !
    ! Note that assuming  A = 1 as e.g. is implied  in the HNC closure
    ! leads to the "ideal gas" result
    !
    !   ε = 1 + 3y
    !
    ! Also A is  independent of the site pair  and A -> 2/3 as  y -> 0
    ! because in this regime ε ~ 1 + 3y + 3y².
    !
    ! [CS81] Exact asymptotic form of the site-site direct correlation
    !   function for rigid polar molecules, Peter T. Cummings and
    !   G. Stell, Molecular Physics: An International Journal at the
    !   Interface Between Chemistry and Physics, Volume 44, Issue 2,
    !   1981, pages 529-531,
    !   http://dx.doi.org/10.1080/00268978100102621
    !
    implicit none
    real (rk), intent (in) :: e, y
    real (rk) :: a
    ! *** end of interface ***

    a = (1 + e * (3 * y - 1)) / (3 * y * (e - 1))
  end function dipole_factor


  function dipole (sites) result (m)
    use foreign, only: site
    implicit none
    type (site), intent (in) :: sites(:)
    real (rk) :: m(3)
    ! *** end of interface ***

    integer :: i

    do i = 1, 3
       m(i) = sum (sites % charge * sites % x(i))
    enddo
  end function dipole


  function weights (sites) result (w)
    !
    ! Site weights (volumes) normalized to 1.
    !
    use foreign, only: site
    implicit none
    type (site), intent (in) :: sites(:)
    real (rk) :: w(size (sites))
    ! *** end of interface ***

    ! Choice of the weights as site "volumes":
    w = (sites % sigma)**3

    w = w / sum (w)
  end function weights


  function center (sites) result (c)
    !
    ! We define the center of a  species as a weighted sum of the site
    ! coordinates with weights proportional  to the "volume" of a site
    ! (σ³).  A  "center" appears in  the discussion of  the long-range
    ! site-site correlations for charged sites in a dipole species.  A
    ! neutral  dipole species  has no  charge center.   A  mass center
    ! appears  to  be  off-topic   in  the  discussion  of  long-range
    ! correlations (masses  never appear in the  equations). This kind
    ! of "geometric" center is the best we can offer here. Suggestions
    ! are welcome.
    !
    use foreign, only: site
    implicit none
    type (site), intent (in) :: sites(:)
    real (rk) :: c(3)
    ! *** end of interface ***

    integer :: i
    real (rk) :: w(size (sites))

    ! Choice of the weights:
    w = weights (sites)

    do i = 1, 3
       c(i) = sum (w * sites % x(i))
    enddo
  end function center


  function axes (sites) result (u)
    use foreign, only: site
    use linalg, only: eigv
    use iso_fortran_env, only: error_unit
    implicit none
    type (site), intent (in) :: sites(:)
    real (rk) :: u(3, 3)
    ! *** end of interface ***

    integer :: s, i, j
    real (rk) :: w(size (sites))
    real (rk) :: d(3), c(3), t(3, 3), e(3)

    ! Geometric center of the species:
    c = center (sites)

    ! The same weights as used in center():
    w = weights (sites)

    t = 0.0
    do s = 1, size (sites)
       ! Site coordinates relative to the center:
       d = sites(s) % x - c

       do j = 1, 3
          do i = 1, 3
             t(i, j) = t(i, j) + w(s) * d(i) * d(j)
          enddo
       enddo
    enddo

    ! u(:, :) will be the eigenvectors, and e(:) will be the
    ! eigenvalues:
    call eigv (t, e, u)

    ! FIXME: do  something about degenerate  shapes. The axes  are not
    ! well defined by  the current procedure in such  cases.  At least
    ! issue a warning ...
    if (minval (abs (e(2:) - e(:2))) < 1.0d-7) then
       write (error_unit, *) " WARNING: degenerate shape, e =", e
    endif
  end function axes


  function dipole_axes (sites) result (u)
    use foreign, only: site
    use linalg, only: eigv
    use iso_fortran_env, only: error_unit
    implicit none
    type (site), intent (in) :: sites(:)
    real (rk) :: u(3, 3)
    ! *** end of interface ***

    real (rk) :: m(3)

    ! Z-axis will be collinear with the dipole moment:
    m = dipole (sites)

    ! FIXME: should be just return axes() in this case?
    if (norm2 (m) < 1.0d-7) then
       write (error_unit, *) " WARNING: zero dipole, m =", m
    endif

    ! m is a unit vector from now on:
    m = m / norm2 (m)

    !
    ! The other two axes should be orthogonal to the dipole vector. It
    ! is  unsettling   that  the  two  remaining   axes  are  somewhat
    ! arbitrary. This is just one  of the many possible ways to choose
    ! them:
    !
    block
       integer :: i, j, k, loc(1)
       real (rk) :: v(3, 3), e(3)

       ! This proposes "shape tensor" eigenvectors as local axes:
       v = axes (sites)

       ! We are going to replace the axis closest to the dipole vector
       ! with the dipole  vector. To find such a  vector, first remove
       ! the collinear component from all three:
       do j = 1, 3
          v(:, j) = v(:, j) - m * dot_product (m, v(:, j))
          e(j) = dot_product (v(:, j), v(:, j))
       enddo

       ! The one  that has the  smallest magnitude was the  closest to
       ! the m-axis:
       loc = minloc (e)
       i = loc(1)

       ! Reorder vectors and make m the last one:
       k = 0
       do j = 1, 3
          if (j == i) cycle
          k = k + 1
          u(:, k) = v(:, j)
       enddo
       u(:, 3) = m

       ! Make u(:, 1)  and u(:, 2) orthogonal to  each other (both are
       ! already orthogonal to m):
       u(:, 2) = u(:, 2) - u(:, 1) &
            * dot_product (u(:, 1), u(:, 2)) / dot_product (u(:, 1), u(:, 1))

       ! Make them unit vectors:
       do j = 1, 3
          u(:, j) = u(:, j) / norm2 (u(:, j))
       enddo
    end block

    ! This is another choice of the two remaining vectors:
    block
       integer :: s, i, j
       real (rk) :: w(size (sites))
       real (rk) :: c(3), d(3), dxy(2), t(2, 2), v(2, 2), e(2)

       ! Geometric center of the species:
       c = center (sites)

       ! Choice of the weights:
       w = weights (sites)

       t = 0.0
       do s = 1, size (sites)
          ! Site coordinates relative to the center:
          d = sites(s) % x - c

          ! x- and y-coordinates given the two axes orthogonal to m:
          do j = 1, 2
             dxy(j) = dot_product (d, u(:, j))
          enddo

          !
          ! Increment the  2x2 "shape thensor". This is  the same kind
          ! of shape tensor we used to propose the default (not dipole
          ! related) axes, albeit restricted to two dimensions:
          !
          do j = 1, 2
             do i = 1, 2
                t(i, j) = t(i, j) + w(s) * dxy(i) * dxy(j)
             enddo
          enddo
       enddo

       !
       ! Here the 2x2 eigenvalue problem is being solved: v(:, :) will
       ! be  the   eigenvectors,  and   e(:)  will  be   the  (unused)
       ! eigenvalues:
       !
       call eigv (t, e, v)

       !
       ! Alternative x- and y-axes are just another linear combination
       ! of the old ones (did I mention they are somewhat arbitrary?):
       !
       ! u   := Σ  u   v
       !  ij     k  ik  kj
       !
       u(1:3, 1:2) = matmul (u(1:3, 1:2), v(1:2, 1:2))

       ! FIXME: Still, if say in water, the hydrogens have zero radius
       ! (weights), they do not  contribute to the "shape tensor" (the
       ! species is  effectively a cylinder  or rather a ball  in this
       ! case) and thus  x- and y- axes are  still arbitrary as coming
       ! out of  the eigensolver.  So far  we just issue  a warning in
       ! such cases:
       if (abs (e(2) - e(1)) < 1.0d-7) then
          write (error_unit, *) " WARNING: degenerate shape, e =", e
       endif
    end block
  end function dipole_axes


  function local_coords (sites, u) result (x)
    !
    ! Returns  orientation   independent  coordinates  of   the  sites
    ! relative to the "center" of the species.
    !
    ! If supplied  the axes u(:,  :) as prepared by  dipole_axes() the
    ! Z-axis is collinear with the dipole moment by construction.
    !
    use foreign, only: site
    implicit none
    type (site), intent (in) :: sites(:)
    real (rk), intent (in) :: u(3, 3) ! orthogonal matrix
    real (rk) :: x(3, size (sites))
    ! *** end of interface ***

    integer :: i, j
    real (rk) :: c(3)

    c = center (sites)

    do j = 1, size (sites)
       do i = 1, 3
          x(i, j) = dot_product (u(:, i), sites(j) % x - c)
       enddo
    enddo
  end function local_coords


  function dipole_correction (beta, rho, eps, sites, k) result (xk)
    !
    ! The   site-site  correction  for   dipole  species   in  Fourier
    ! representation (see e.g. Ref. [KH00a]):
    !
    !                                     2
    !   χ  (k) = f (k) h (k) f (k)  ~  o(k )
    !    ab       a     c     b
    !
    ! with a short k-range
    !
    !                                         2
    !   h (k) = (ε - 1 - 3y) / (ρy) exp[- (ks) / 4]  ~  o(1)
    !    c
    !
    ! and the site-specific functions
    !
    !   f (k) = j (kx ) j (ky ) j (kz )  ~  kz / 3
    !    a       0   a   0   a   1   a        a
    !
    ! where the site parameters
    !
    !   (x , y , z )
    !     a   a   a
    !
    ! are the site coordinates bound to a local coordiante system with
    ! z-axis collinear to the dipole  vector. The center of the system
    ! and  the two  other  axes appear  to  be left  arbitrary in  the
    ! literature (I would like a counterexample).
    !
    ! Note that by construction the charge-weighted sum
    !
    !   Σ  q  f (k)  ~  (k / 3) Σ q  z   =  kμ / 3
    !    a  a  a                 a a  a
    !
    ! at small k. And thus, the double sum
    !
    !                                      2
    !  Σ   q  ρ χ  (k) ρ q   ~  h (0) (kρμ) / 9
    !   ab  a    ab       b      c
    !
    !                            2
    !                        =  k  [ε - 1 - 3y] / 4πβ
    !
    !
    ! [KH00a] Potentials of mean force of simple ions in ambient
    !   aqueous solution. I. Three-dimensional reference interaction
    !   site model approach, Andriy Kovalenko and Fumio Hirata,
    !   J. Chem. Phys. 112, 10391 (2000);
    !   http://dx.doi.org/10.1063/1.481676
    !
    use foreign, only: site
    use bessel, only: j0, j1
    implicit none
    real (rk), intent (in) :: beta, rho
    real (rk), intent (in) :: eps        ! target epsilon
    type (site), intent (in) :: sites(:) ! (m)
    real (rk), intent (in) :: k(:)       ! (nrad)
    real (rk) :: xk(size (k), size (sites), size (sites)) ! (nrad, m, m)
    ! *** end of inteface ***

    real (rk), parameter :: s = 0.5 * ANGSTROM ! Ref. [KH00a]
    integer :: i, j
    real (rk) :: x(3, size (sites))
    real (rk) :: f(size (k), size (sites))
    real (rk) :: h(size (k))
    real (rk) :: y              ! dipole density

    ! Make a choice for the center of the species and local axes, then
    ! compute  the local  coordinates.  The (third)  z-axis should  be
    ! collinear with the dipole vector:
    x = local_coords (sites, dipole_axes (sites))

    ! Precompute j0 * j0 * j1 for each site:
    do i = 1, size (sites)
       associate (d => x(:, i))
         f(:, i) = j0 (k * d(1)) * j0 (k * d(2)) * j1 (k * d(3))
       end associate
    enddo

    y = dipole_density (beta, rho, sites)

    ! h(k), common for all sites. Note the the term is proportional to
    ! the difference  of target and rism dielectric  constants --- the
    ! unmodified  RISM theory  gives  e =  1  + 3y  as the  dielectric
    ! constant:
    h = (eps - (1 + 3 * y)) / (rho * y) * exp (- (k * s)**2 / 4)

    do j = 1, size (sites)
       do i = 1, size (sites)
          xk(:, i, j) = f (:, i) * h(:) * f(:, j)
       enddo
    enddo
  end function dipole_correction


  subroutine gnuplot1 (r, g)
    implicit none
    real (rk), intent (in) :: r(:) ! (nrad)
    real (rk), intent (in) :: g(:) ! (nrad)
    ! *** end of interface ***

    integer :: p

    ! This prints a lot of data on tty!
    print *, "# r then g"
    do p = 1, size (g)
       write (*, *) r(p), g(p)
    enddo
  end subroutine gnuplot1


  subroutine gnuplot3 (r, g)
    implicit none
    real (rk), intent (in) :: r(:)       ! (nrad)
    real (rk), intent (in) :: g(:, :, :) ! (nrad, n, m)
    ! *** end of interface ***

    integer :: i, j, p

    ! This prints a lot of data on tty!
    print *, "# r then g(i, j) for",  size (g, 2), "x", size (g, 3), "pairs"
    do p = 1, size (g, 1)
       write (*, *) r(p), ((g(p, i, j), i = 1, size (g, 2)), j = 1, size (g, 3))
    enddo
  end subroutine gnuplot3

end module rism

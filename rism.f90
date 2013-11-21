module rism
  use kinds, only: rk
  implicit none
  private

  ! These are  all bind(c)  and are not  used in Fortran  sources, but
  ! from C-side:
  public :: rism_solvent
  public :: rism_solute
  public :: rism_nrad
  public :: rism_rmax
  public :: rism_upscale
  ! *** END OF INTERFACE ***

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


  function rism_self_energy (n, sites, spec) result (e) bind (c)
    !
    ! Calculate the energy summation between each pair of residues for
    ! a given site.
    !
    ! Needs to be consistent with ./rism.h
    !
    use iso_c_binding, only: c_int, c_double
    use foreign, only: site
    implicit none
    integer (c_int), intent (in), value :: n
    type (site), intent (in) :: sites(n)
    integer (c_int), intent (in) :: spec(n)
    real (c_double) :: e
    ! *** end of interface ***

    e = self_energy (sites, spec)
  end function rism_self_energy


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

    ! FIXME:  The  code operates  with  arrays  of  the [m,  m,  nrad]
    ! shape. However,  for historical reasons  the C-interface assumes
    ! the  transposed shape.   Conversion between  the two  layouts is
    ! perfomed by main().
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

    ! FIXME:  The  code operates  with  arrays  of  the [m,  m,  nrad]
    ! shape. However,  for historical reasons  the C-interface assumes
    ! the  transposed shape.   Conversion between  the two  layouts is
    ! perfomed by main().
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
    use drism, only: epsilon_rism
    implicit none
    type (problem_data), intent (in) :: pd
    type (site), intent (in) :: solvent(:)
    type (site), optional, intent (in) :: solute(:)
    real (rk), optional, intent (out) :: t(:, :, :) ! (nrad, m, m), sic!
    real (rk), optional, intent (inout) :: x(:, :, :) ! (nrad, m, m), sic!
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
       ! Indirect correlation t = h - c aka γ (hence the name):
       real (rk) :: gam(size (solvent), size (solvent), nrad)

       ! Solvent susceptibility χ = ω + ρh:
       real (rk) :: chi(size (solvent), size (solvent), nrad)

       ! Procedures  down  the  call  chain  will set  this  to  valid
       ! assiciation lists:
       type (obj) :: vdict, udict

       if (vv) then
          call rism_vv (pd % closure, nrad, rmax, pd % beta, pd % rho, &
               solvent, gam, chi, vdict)

          ! Output, if requested:
          if (present (t)) then
             t = clayout (gam)
          endif

          if (present (x)) then
             x = clayout (chi)
          endif
       else
          if (.not. present (x)) error stop "no way to get chi!"
          chi = flayout (x)
          vdict = nil
       endif

       if (uv) then
          call rism_uv (pd % closure, nrad, rmax, pd % beta, pd % rho, &
               solvent, chi, solute, udict)
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
    use fft, only: fourier_rows, FT_FW, FT_BW
    use snes, only: snes_default
    use foreign, only: site, comm_rank
    use options, only: getopt
    use lisp, only: obj
    use units, only: pi
    use drism, only: dipole_density, dipole_factor, dipole_correction, &
         epsilon_rism
    implicit none
    integer, intent (in) :: method          ! HNC, KH, or PY
    integer, intent (in) :: nrad            ! grid size
    real (rk), intent (in) :: rmax          ! cell size
    real (rk), intent (in) :: beta          ! inverse temp
    real (rk), intent (in) :: rho
    type (site), intent (in) :: sites(:)    ! (m)
    real (rk), intent (out) :: gam(:, :, :) ! (m, m, nrad)
    real (rk), intent (out) :: chi(:, :, :) ! (m, m, nrad)
    type (obj), intent (out) :: dict
    ! *** end of interface ***


    ! Pair quantities. FIXME: they are symmetric, one should use that:
    real (rk), dimension (size (sites), size (sites), nrad) :: &
         vr, vk, wk, xk, t, c

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
    call force_field (sites, sites, r, k, vr, vk)

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
    call snes_default (t, iterate_t)

    ! Do not assume c has a meaningfull value, it was overwritten with
    ! c(k):
    c = closure (method, beta, vr, t)

    !
    ! Return the  indirect correlation t(r)  aka γ(r) in real  rep and
    ! solvent susceptibility χ(k) in Fourier rep:
    !
    gam = t

    block
       real (rk) :: h(m, m, nrad)

       ! h(k) = c(k) + t(k)
       h = fourier_rows (c + t) * (dr**3 / FT_FW)

       ! χ = ω + ρh
       chi = wk + rho * h
    end block

    ! Done with it, print results. Here solute == solvent:
    call post_process (method, beta, rho, sites, sites, dr, dk, vr, t, &
         A=A, eps=eps, dict=dict, rbc=.false.)

  contains

    function iterate_t (t) result (dt)
      !
      ! Closure over  host variables: r, k,  dr, dk, v,  c, beta, rho,
      ! ... Implements procedure(func1).
      !
      implicit none
      real (rk), intent (in) :: t(:, :, :) ! (m, m, nrad)
      real (rk) :: dt(size (t, 1), size (t, 2), size (t, 3))
      ! *** end of interface ***

      c = closure (method, beta, vr, t)

      ! Forward FT via DST:
      c = fourier_rows (c) * (dr**3 / FT_FW)

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
      dt = fourier_rows (dt) * (dk**3 / FT_BW)

      ! Return the increment that vanishes at convergence:
      dt = dt - t
    end function iterate_t
  end subroutine rism_vv


  subroutine rism_uv (method, nrad, rmax, beta, rho, solvent, chi, solute, dict)
    use snes, only: snes_default
    use foreign, only: site, verbosity
    use lisp, only: obj, cons, sym, num
    use options, only: getopt
    use units, only: pi, angstrom
    implicit none
    integer, intent (in) :: method         ! HNC, KH, or PY
    integer, intent (in) :: nrad           ! grid size
    real (rk), intent (in) :: rmax         ! cell size
    real (rk), intent (in) :: beta         ! inverse temp
    real (rk), intent (in) :: rho
    type (site), intent (in) :: solvent(:) ! (m)
    real (rk), intent(in) :: chi(:, :, :)  ! (m, m, nrad)
    type (site), intent (in) :: solute(:)  ! (n)
    type (obj), intent (out) :: dict
    ! *** end of interface ***

    ! Solute-solvent pair quantities:
    real (rk), dimension (size (solute), size (solvent), nrad) :: &
         v_uvr, v_uvk, t_uvx, c_uvx, expB ! (n, m, nrad)

    ! Solute-solute pair quantities:
    real (rk), dimension (size (solute), size (solute), nrad) :: w_uuk ! (n, n, nrad)

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
    call force_field (solute, solvent, r, k, v_uvr, v_uvk)

    ! Rigid-bond solute-solute correlations on the k-grid:
    w_uuk = omega_fourier (solute, k)

    rbc = getopt ("rbc")

    if (rbc) then
       call bridge (solute, solvent, beta, r, dr, k, dk, expB)
    else
       ! Not used in this case:
       expB = 1.0
    endif


    ! Intitial guess:
    t_uvx = 0.0

    ! Find t such that iterate_t  (t) == 0. FIXME: passing an internal
    ! function as a callback is an F08 feature. GFortran 4.3 on Debian
    ! Lenny does not support that:
    call snes_default (t_uvx, iterate_t)

    ! Done with it, print results:
    call post_process (method, beta, rho, solvent, solute, dr, dk, v_uvr, t_uvx, &
         A=1.0d0, eps=0.0d0, dict=dict, rbc=rbc)

    !
    ! As  a part  of  post-processing, compute  the bridge  correction
    ! using TPT.
    !
    ! As  claimed  in  [KH00c],  in  general  the  solvation  chemical
    ! potential is no  longer of an analytic form  if a nonzero bridge
    ! correction  is  included  into  the  closure.  Direct  numerical
    ! evaluation  of   chemical  potential  can  be   carried  out  by
    ! integration over the LJ  diameter parameter. See expression (20)
    ! in [KH00c].
    !
    ! So with  our current implementation,  calculations involving RBC
    ! could be split into two parts:
    !
    ! 1. Don't apply bridge correction when solving t, only add energy
    !    contribution to excess chemical potential in the spirit of
    !    thermodynamic perturbation theory (TPT).  In this case the
    !    value of solvaiton chemical potential is corrected but
    !    correlation function is not perturbed, for some solutes, e.g.
    !    OH-, one can observe an unphsically high peak for
    !    O(OH-)-H(H2O) pair.
    !
    ! 2. In order to get a reasonable (water) solvent structure,
    !    especially for those solutes which have an atomic site of
    !    large negative charge, we need to "switch on" RBC in closure
    !    when solving OZ equation iteratively. However, in this case
    !    evaluation of chemical potential is not implemented (FIXME).
    !
    ! [KH00c] Hydration free energy of hydrophobic solutes studied by
    !   a reference interaction site model with a repulsive bridge
    !   correction and a thermodynamic perturbation method, Andriy
    !   Kovalenko and Fumio Hirata , J. Chem. Phys. 113, 2793 (2000);
    !   http://dx.doi.org/10.1063/1.1305885
    !
    block
      real (rk) :: e, h(n, m, nrad), expB(n, m, nrad)

      call bridge (solute, solvent, beta, r, dr, k, dk, expB)

      ! I. Call closure (without RBC TPT correction): h = c + t
      h = closure (method, beta, v_uvr, t_uvx) + t_uvx

      e = chempot_bridge (beta, rho, h, expB, r, dr)

      if (verbosity > -1) then
         print *, "# TPT bridge correction =", e, "rbc =", rbc
         ! Add a warning when distribution is perturbed by RBC
         if (rbc) then
           print *, "# WARNING: distribution is perturbed, RBC/TPT is approximate"
         endif
      endif

      ! Cons  a dictionary  entry with  RBC TPT  correction  to result
      ! collection. Note that the correction is "wrong" when RDFs were
      ! obtained in an SCF procedure with rbc == .true.
      dict = cons (cons (sym ("RBC-TPT"), num (e)), dict)

      ! II. Call closure with RBC TPT correction: h = c + t
      h = closure_rbc (method, beta, v_uvr, t_uvx, expB) + t_uvx

      e = chempot_bridge (beta, rho, h, expB, r, dr)

      if (verbosity > -1) then
         print *, "# SCF bridge correction =", e, "rbc =", rbc
         ! Add a warning when iterations were performed without RBC
         if (.not. rbc) then
           print *, "# WARNING: distribution is perturbed, RBC/SCF is approximate"
         endif
      endif

      ! Cons  a dictionary  entry with  RBC SCF  correction  to result
      ! collection.  Note  that this number  may be still  "wrong" for
      ! two reasons:  (a) when rbc  = .false.  the RDFs  are perturbed
      ! only  in  a  post-SCF  fashion  and (b)  literature  seems  to
      ! advertise  the   use  of  RBC   in  thermodynamic  integration
      ! procedure (RBC-TI) instead.
      dict = cons (cons (sym ("RBC-SCF"), num (e)), dict)
    end block

    block
      real (rk) :: se
      integer :: spec(n), i

      ! Get  species ID.   Sites belong  to  the same  species if  the
      ! distance  is below  some fixed  value. FIXME:  here  a literal
      ! constant:
      spec = identify_species (solute, 2.0 * angstrom)

      if (verbosity > 1) then
         do i = 1, n
            print *, solute(i) % name, spec(i)
         enddo
      endif

      se = self_energy (solute, spec)

      print *, "# Self energy = ", se, " kcal/mol"

      dict = cons (cons (sym ("SELF-ENERGY"), num (se)), dict)
    end block

  contains

    function iterate_t (t) result (dt)
      !
      ! Closure over  host variables: r, k,  dr, dk, v,  c, beta, rho,
      ! ... Implements procedure (func1).
      !
      use fft, only: fourier_rows, FT_FW, FT_BW
      implicit none
      real (rk), intent (in) :: t(:, :, :) ! (n, m, nrad)
      real (rk) :: dt(size (t, 1), size (t, 2), size (t, 3))
      ! *** end of interface ***

      if (rbc) then
         c_uvx = closure_rbc (method, beta, v_uvr, t, expB)
      else
         c_uvx = closure (method, beta, v_uvr, t)
      endif

      ! Forward FT via DST:
      c_uvx = fourier_rows (c_uvx) * (dr**3 / FT_FW)

      !
      ! The  real-space representation  encodes  only the  short-range
      ! part  of  the  direct   corrlation.  The  (fixed)  long  range
      ! contribution is added here:
      !
      !   C := C  - βV
      !         S     L
      !
      c_uvx = c_uvx - beta * v_uvk

      !
      ! OZ  equation,  involves   "convolutions",  take  care  of  the
      ! normalization here.   As a  functional of c  this is  a linear
      ! relation  t =  C'(c).  FIXME:  Because of  this  linearity one
      ! could have handled the long range term added to c above and to
      ! resulting t below differently:
      !
      !   t := C'(c + x) + x = C'(c) + [C'(x) + x]
      !
      ! with x  = -βv.   Here the second  term in the  square brackets
      ! derived  from  the  fixed  long-range assymptotics  of  direct
      ! correlation  is constant.   This  would be  just one  addition
      ! after computing C'(c) instead of one before and one after.
      !
      dt = oz_uv_equation_c_t (c_uvx, w_uuk, chi)

      !
      ! Since we plugged  in the Fourier transform of  the full direct
      ! correlation including the long range part into the OZ equation
      ! what we get out is the full indirect correlation including the
      ! long-range part.   Note however  that assymptotically C  and T
      ! differ  by just a  sign ---  the menmonic  is C  + T  is short
      ! range.  Take the assymptote of T out making it short range:
      !
      !   T  := T - βV
      !    S          L
      !
      dt = dt - beta * v_uvk

      ! Inverse FT via DST:
      dt = fourier_rows (dt) * (dk**3 / FT_BW)

      ! Return the increment that vanishes at convergence:
      dt = dt - t
    end function iterate_t

    function jacobian_t (t, dt) result (ddt)
      !
      ! Closure over  host variables: r, k,  dr, dk, v,  c, beta, rho,
      ! ... Implements procedure (func2).
      !
      use fft, only: fourier_rows, FT_FW, FT_BW
      implicit none
      real (rk), intent (in) :: t(:, :, :)  ! (n, m, nrad)
      real (rk), intent (in) :: dt(:, :, :) ! (n, m, nrad)
      real (rk) :: ddt(size (t, 1), size (t, 2), size (t, 3))
      ! *** end of interface ***

      ! In this  procedure the variable c_uvx has  a different meaning
      ! than in iterate_t(). Here  it is a *differential* increment dc
      ! of c due to t -> t + dt.
      if (rbc) then
         ! c_uvx = closure_rbc1 (method, beta, v_uvr, t, dt, expB)
         stop "not implemented"
      else
         c_uvx = closure1 (method, beta, v_uvr, t, dt)
      endif

      ! Forward FT via DST:
      c_uvx = fourier_rows (c_uvx) * (dr**3 / FT_FW)

      !
      ! OZ  equation,  involves   "convolutions",  take  care  of  the
      ! normalization here.   As a  functional of c  this is  a linear
      ! relation. See the comments in interate_t() on why the constant
      ! long-range term does not contribute to the Jacobian.
      !
      ddt = oz_uv_equation_c_t (c_uvx, w_uuk, chi)

      ! Inverse FT via DST:
      ddt = fourier_rows (ddt) * (dk**3 / FT_BW)

      ! Return the increment that vanishes at convergence:
      ddt = ddt - dt
    end function jacobian_t
  end subroutine rism_uv


  subroutine bridge (solute, solvent, beta, r, dr, k, dk, expB)
    !
    ! Calculate repulsive bridge correction exp(B):
    !
    !
    !   exp[B  (r)]  =  Π   exp[-4βε  (σ  / r)¹²] * ω  (r)
    !        ij        l≠j          il  il           lj
    !
    ! Here i is a  solute site index and j and l  are the solvent site
    ! indices.
    !
    ! FIXME: there is some room for improvement with the assymptotics.
    ! If one  assumes that the  convolution with ω(r) does  not change
    ! assymptotics, the bridge function should behave as
    !
    !  B  (r) ~ - β  Σ  v  (r)  ~  o(1 / r¹²)
    !   ij          l≠j  il
    !
    ! towards infinity.  However, in the chem. pot. density expression
    ! as it  appears in TPT there  is some "ringing" on  that level if
    ! you look at r¹² * g * [exp(B) - 1] and it is not due to g (which
    ! may be verified by setting  g = 1). Fortunately enough the shape
    ! of 4πr² * g * [exp(B) - 1] is satisfactorily smooth and decaying
    ! fast enough.
    !
    use fft, only: fourier_rows, FT_BW, FT_FW
    use foreign, only: site
    implicit none
    type (site), intent (in) :: solute(:)     ! (n)
    type (site), intent (in) :: solvent(:)    ! (m)
    real (rk), intent (in) :: beta            ! inverse temperature
    real (rk), intent (in) :: r(:)            ! (nrad)
    real (rk), intent (in) :: k(:)            ! (nrad)
    real (rk), intent (in) :: dr, dk
    real (rk), intent (out) :: expB(:, :, :) ! (n, m, nrad)
    ! *** end of interface ***

    integer :: m, n

    m = size (solvent)
    n = size (solute)

    block
      ! Solvent-solvent pair quantities:
      real (rk), dimension (m, m, size (r)) :: w

      ! Storage for solute-solvent pair quantities (dont take the name
      ! h literally):
      real (rk), dimension (n, m, size (r)) :: f, h

      ! Solvent  rigid-bond   correlations  ω  on   the  k-grid.   The
      ! self-correlation  is  a  δ-function   (or  1  in  the  Fourier
      ! representation).
      w = omega_fourier (solvent, k)

      ! Repulsive part of LJ potential goes to f(i, j, :). First index
      ! is that of a solute site  i, second index is of a solvent site
      ! j:
      call lj_repulsive (solute, solvent, r, f)

      !
      ! A  Mayer's function  f =  exp(-βv) -  1 due  to  the repulsive
      ! branch of the potential.  We will assume that convolution of a
      ! constant with intramolecular correlation ω(r) is unchanged:
      !
      !   ω * (1 + f) = 1 + ω * f
      !
      ! Such a convolution corresponds  to an averaging over sphere at
      ! every point.
      !
      f = expm1 (-beta * f)

      ! Fourier transform of  the Mayer's factor exp(-βv) -  1.  It is
      ! this  factor   (up  to  a  constant)  which   appears  in  the
      ! convolution with the intra-molecular solvent-solvent site-site
      ! correlation ω:
      f = fourier_rows (f) * (dr**3 / FT_FW)

      ! Compute expB(i,  j, :) as a  product over all  solvent sites l
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
                 h(i, j, :) =  f(i, l, :) * w(l, j, :)
              enddo
           enddo

           ! Transform convolutions to the real space:
           h = fourier_rows (h) * (dk**3 / FT_BW)

           ! Here the  product is accumulated.  FIXME:  Note that even
           ! though the factors  with l == j are  computed above, they
           ! are excluded from the product:
           do j = 1, m          ! solvent sites
              do i = 1, n       ! solute sites
                 if (j == l) cycle
                 expB(i, j, :) = expB(i, j, :) * (1 + h(i, j, :))
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
    real (rk), intent (in) :: h(:, :, :)    ! (n, m, nrad)
    real (rk), intent (in) :: expB(:, :, :) ! (n, m, nrad)
    real (rk), intent (in) :: r(:)          ! (nrad)
    real (rk), intent (in) :: dr
    real (rk) :: mu
    ! *** end of interface ***

    integer :: i, j, p
    real (rk) :: acc, density (size (r))

    ! According   to  thermodynamic   perturbation   theory,  chemical
    ! potential density to be integrated:
    do p = 1, size (h, 3)
       acc = 0.0
       do j = 1, size (h, 2)
          do i = 1, size (h, 1)
             acc = acc + (expB(i, j, p) - 1) * (1 + h(i, j, p))
          enddo
       enddo
       density(p) = acc
    enddo

    ! The integration  procedure assumes a  very specific radial  (i +
    ! 1/2) grid.  Multiply that by dr³ and divide by β to get the real
    ! number:
    mu = integrate (density) * (rho * dr**3 / beta)
  end function chempot_bridge


  subroutine pair (a, b, sigma, epsilon)
    !
    ! Place to implement combination rules.
    !
    use foreign, only: site
    use options, only: getopt
    implicit none
    type (site), intent (in) :: a, b
    real (rk), intent (out) :: sigma, epsilon
    ! *** end of interface ***

    integer :: rule

    if (.not. getopt ("comb-rule", rule)) then
       rule = 0
    endif

    select case (rule)
    case (0)
       ! Arithmetic average for sigma:
       sigma = (a % sigma + b % sigma) / 2
       epsilon = sqrt (a % epsilon * b % epsilon)
    case (1)
       ! Both are geometric averages.  The same rule as with geometric
       ! averages of 6-12 coefficients:
       sigma = sqrt (a % sigma * b % sigma)
       epsilon = sqrt (a % epsilon * b % epsilon)
    case default
       sigma = -1.0
       epsilon = huge (epsilon)
       error stop "no such rule!"
    end select
  end subroutine pair


  subroutine lj_repulsive (asites, bsites, r, vr)
    !
    ! only calculate the repulsive part of LJ potential (r^-12 term)
    !
    use foreign, only: site
    implicit none
    type (site), intent (in) :: asites(:)       ! (n)
    type (site), intent (in) :: bsites(:)       ! (m)
    real (rk), intent (in) :: r(:)              ! (nrad)
    real (rk), intent (out) :: vr(:, :, :)      ! (n, m, nrad)
    ! ** end of interface ***

    real (rk) :: sigma, epsilon
    integer :: i, j

    do j = 1, size (bsites)
       do i = 1, size (asites)
          ! Derive parameters of pair interaction:
          call pair (asites(i), bsites(j), sigma, epsilon)

          if (sigma /= 0.0) then
              vr (i, j, :) = epsilon * lj12 (r / sigma)
          else
              vr (i, j, :) = 0.0
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
       A, eps, dict, rbc)
    !
    ! Prints some results.
    !
    use fft, only: fourier_rows, FT_FW, integral
    use linalg, only: polyfit
    use foreign, only: site, verbosity, comm_rank, &
         HNC => CLOSURE_HNC, KH => CLOSURE_KH, PY => CLOSURE_PY
    use lisp, only: obj, cons, nil, sym, num
    use units, only: pi, ALPHA, EPSILON0INV, KCAL, KJOULE
    use drism, only: dipole, center, dipole_axes, local_coords, dipole_density, &
         epsilon_rism, dipole_factor, dipole_correction
    implicit none
    integer, intent (in) :: method         ! HNC, KH or PY
    real (rk), intent (in) :: beta         ! inverse temperature
    real (rk), intent (in) :: rho
    type (site), intent (in) :: solvent(:) ! (m)
    type (site), intent (in) :: solute(:)  ! (n)
    real (rk), intent (in) :: dr, dk       ! grid steps
    real (rk), intent (in) :: v(:, :, :)   ! (n, m, nrad)
    real (rk), intent (in) :: t(:, :, :)   ! (n, m, nrad)
    real (rk), intent (in) :: A            ! long-range scaling factor
    real (rk), value :: eps     ! requested dielectric constant, or 0
    type (obj), intent (out) :: dict
    logical, value :: rbc
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

    n = size (t, 1)             ! number of solute sites
    m = size (t, 2)             ! number of solvent sites
    nrad = size (t, 3)          ! number of radial point

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
       real (rk) :: c(n, m, nrad)
       real (rk) :: cl(n, m, nrad)
       real (rk) :: h(n, m, nrad)
       real (rk) :: g(n, m, nrad)
       real (rk) :: x(nrad), x0(nrad), x1(nrad), x2(nrad), xx(nrad) ! qhq in k-space
       real (rk) :: xd(m, m, nrad)
       real (rk) :: expB(n, m, nrad)
       real (rk) :: ni(n, m, nrad) ! number integral

       ! Dont like to pass redundant  info, recompute r(:) from dr and
       ! k(:) from dk:
       forall (i = 1:nrad)
          r(i) = (2 * i - 1) * dr / 2
          k(i) = (2 * i - 1) * dk / 2
       end forall

       ! For the same reason recomute c and g:
       if (rbc) then
          call bridge (solute, solvent, beta, r, dr, k, dk, expB)
          c = closure_rbc (method, beta, v, t, expB)
       else
          expB = 1.0
          c = closure (method, beta, v, t)
       endif
       h = c + t
       g = 1 + h

       ! Real-space rep of the  long-range correlation. Note the extra
       ! scaling factor A:
       forall (p = 1:nrad, i = 1:n, j = 1:m)
          cl(i, j, p) = - (beta * A) * solute(i) % charge * solvent(j) % charge &
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
          real (rk) :: qu(n), qv(m)
          real (rk) :: y, fac0, fac1, fac2
          real (rk) :: qhq(nrad)
          real (rk) :: hk(n, m, nrad)

          ! Small-k  behavior   of  qh(k)q  which   is  essential  for
          ! dielectric permittivity:
          hk = fourier_rows (h) * (dr**3 / FT_FW)

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
          qv(:) = rho * solvent(:) % charge
          qu(:) = rho * solute(:) % charge

          ! Note  the  extra  EPSILON0INV   factor,  it  seems  to  be
          ! necessary to get the kcal/A³ units right:
          qhq = 0.0
          x = 0.0
          xx = 0.0
          do i = 1, size (solute)
             do j = 1, size (solvent)
                do p = 1, nrad
                   qhq(p) = qhq(p) + qu(i) * h(i, j, p) * qv(j) * EPSILON0INV
                   x(p) = x(p) + qu(i) * hk(i, j, p) * qv(j) * EPSILON0INV
                enddo
             enddo
          enddo

          xd = dipole_correction (beta, rho, eps, solvent, k)

          xx = 0.0
          do i = 1, size (solvent)
             do j = 1, size (solvent)
                do p = 1, nrad
                   xx(p) = xx(p) + qv(i) * xd(i, j, p) * qv(j) * EPSILON0INV
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


       ! Compute number integral
       forall (i = 1:n, j = 1:m)
          ni(i, j, :) = integral (g(i, j, :)) * (rho * dr**3)
       end forall

       ! This prints a lot of data on tty!
       if (verb > 1) then
          print *, "# r(i) then g(i), ni(i), v(i), t(i), c(i), each for",  n, "x", m, "pairs"
          do p = 1, nrad
             write (*, *) r(p), &
                  &     ((g(i, j, p), i=1,n), j=1,m), &
                  &     ((ni(i, j, p), i=1,n), j=1,m), &
                  &     ((v(i, j, p), i=1,n), j=1,m), &
                  &     ((t(i, j, p), i=1,n), j=1,m), &
                  &     ((c(i, j, p), i=1,n), j=1,m), &
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
      use units, only: pi
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


  function omega_fourier (sites, k) result (w)
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
    type (site), intent (in) :: sites(:)                 ! (m)
    real (rk), intent (in) :: k(:)                       ! (nrad)
    real (rk) :: w(size (sites), size (sites), size (k)) ! (m, m, nrad)
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
          w(i, j, :) = j0 (k * rab)
       enddo
    enddo
  end function omega_fourier


  function omega_fourier1 (sites, k, n, dx) result (dw)
    !
    ! Differential of ω(k) in response do displacment dx of site n.
    !
    use foreign, only: site
    use bessel, only: j1
    implicit none
    type (site), intent (in) :: sites(:) ! (m)
    real (rk), intent (in) :: k(:)       ! (nrad)
    integer, intent (in) :: n            ! moving site
    real (rk), intent (in) :: dx(3)      ! site displacement
    real (rk) :: dw(size (sites), size (sites), size (k)) ! (m, m, nrad)
    ! *** end of inteface ***

    real (rk) :: xa(3), xb(3), xab(3), rab, drab
    integer :: i, m

    m = size (sites)

    ! Only  the the  n-th row  and n-th  column of  the matrix  dω are
    ! non-zero.  The diagonal element at  (n, n) is zero though (since
    ! diagonal elements  of ω are  1, their derivative is  zero).  All
    ! other elements vanish:
    dw(:, :, :) = 0.0

    ! Coordinates of the moving site:
    xa = sites(n) % x

    ! Loop over other sites, i /= n:
    do i = 1, m
       if (i == n) cycle

       xb = sites(i) % x

       xab = xa - xb
       rab = norm2 (xab)
       drab = dot_product (xab, dx) / rab

       ! d/dx j0(x) = -j1(x):
       dw(i, n, :) = -j1 (k * rab) * k * drab

       ! The matrix and its differential are symmetric:
       dw(n, i, :) = dw(i, n, :)
    enddo
  end function omega_fourier1


  subroutine force_field (asites, bsites, r, k, vr, vk)
    !
    ! Force field  is represented by two contributions,  the first one
    ! is (hopefully) of short range and is returned in array vr on the
    ! r-grid. The other term is  then of long range and, naturally, is
    ! represented by array vk on the k-grid.
    !
    use foreign, only: site
    use units, only: EPSILON0INV, ALPHA
    implicit none
    type (site), intent (in) :: asites(:)  ! (n)
    type (site), intent (in) :: bsites(:)  ! (m)
    real (rk), intent (in) :: r(:)         ! (nrad)
    real (rk), intent (in) :: k(:)         ! (nrad)
    real (rk), intent (out) :: vr(:, :, :) ! (n, m, nrad)
    real (rk), intent (out) :: vk(:, :, :) ! (n, m, nrad)
    ! *** end of inteface ***

    real (rk) :: epsilon, sigma, charge
    integer :: i, j
    type (site) :: a, b

    do j = 1, size (bsites)
       b = bsites(j)
       do i = 1, size (asites)
          a = asites(i)
          ! Derive pair interaction parameters:
          call pair (a, b, sigma, epsilon)
          charge = a % charge * b % charge

          ! Short range on the r-grid:
          if (sigma /= 0.0) then
             vr(i, j, :) = epsilon * lj (r / sigma) + &
                  EPSILON0INV * charge * coulomb_short (r, ALPHA)
          else
             vr(i, j, :) = &
                  EPSILON0INV * charge * coulomb_short (r, ALPHA)
          endif

          ! Long range on the k-grid:
          vk(i, j, :) = &
               EPSILON0INV * charge * coulomb_long_fourier (k, ALPHA)
       enddo
    enddo

  end subroutine force_field


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


  function distance_matrix (sites) result (rij)
    !
    ! Calculate distance matrix of each atomic pair for a given site.
    !
    use foreign, only: site
    implicit none
    type (site), intent (in) :: sites(:)         ! (m)
    real (rk) :: rij(size (sites), size (sites)) ! (m, m)
    ! *** end of interface ***

    integer i, j, m

    m = size (sites)

    do j = 1, m
       do i = 1, m
          rij(i, j) = norm2 (sites(i) % x - sites(j) % x)
       enddo
    enddo
  end function distance_matrix


  function identify_species (sites, r0) result (spec)
    !
    ! Given the distance matrix rij(m, m), get the species ID for each
    ! site spec(m). If the distance  between two sites is less than r0
    ! (~ 2.0 A), then they belong to the same species.
    !
    use foreign, only: site
    implicit none
    type (site), intent (in) :: sites(:) ! (m)
    real (rk), intent (in) :: r0
    integer :: spec(size (sites)) ! (m)
    ! *** end of interface ***

    real (rk) :: rij(size (sites), size (sites)) ! (m, m)
    integer i, j, k, species

    ! Get distance matrix
    rij = distance_matrix (sites)

    ! Initialize species ID to an  invalid ID (here 0).  Valid IDs are
    ! positive. There is zero identified species at the moment:
    spec(:) = 0
    species = 0
    do j = 1, size (sites)
       ! First  see  if the  current  site  j  belongs to  an  already
       ! identified  species.   Only  loop over  previousely  assigned
       ! sites and skip the diagonal i == j:
       do i = 1, j - 1
          !
          ! Sites i and j belong to  the same species if r(i, j) < r0.
          ! FIXME:  for  the  current  study (aqueous  solution)  it's
          ! first-come-first-served.    What  if   a  site   is  close
          ! connected with more than  one neighbours while they belong
          ! to different species?
          !
          ! It  may happen  that  upon  addition of  the  site j  that
          ! connects two  species into one the number  of species goes
          ! down.  The sites of  the previousely separate species need
          ! to be re-labeled:
          !
          if (rij(i, j) < r0) then
             if (spec(j) == 0) then
                spec(j) = spec(i)
             else if (spec(j) /= spec(i)) then
                ! All sites of species  i need to be re-labled: FIXME:
                ! verify when encontered and remove the statement:
                error stop "unverified yet!"
                do k = 1, j
                   if (spec(k) /= spec(i)) cycle
                   spec(k) = spec(j)
                enddo
             endif
             ! At this point spec(i) == spec(j) and both are valid.
          else
             ! Otherwise both spec(j) and spec(i) remain unchanged.
          endif
       enddo

       ! If site j was assigned, proceed to the next:
       if (spec(j) /= 0) cycle

       ! Otherwise we  have found a  new species --  increment species
       ! count and assign new ID:
       species = species + 1
       spec(j) = species
    enddo
    if (any (spec == 0)) error stop "cannot happen!"
  end function identify_species


  function self_energy (sites, spec) result (e)
    !
    ! Calculate the energy summation between each pair of residues for
    ! a given site.
    !
    use foreign, only: site
    use units, only: ALPHA, EPSILON0INV
    implicit none
    type (site), intent (in) :: sites(:) ! (m)
    integer, intent (in) :: spec(:)      ! (m)
    real (rk) :: e
    ! *** end of interface ***

    integer i, j, m

    m = size (sites)

    ! Loop over upper triangle with i <= j <= m, not to count the same
    ! interaction twice.
    e = 0.0
    do j = 1, m
       do i = 1, j
          ! Self-interaction  of each  site is  automatically excluded
          ! because of this:
          if (spec(i) == spec(j)) cycle
          e = e + pair_energy (sites(i), sites(j))
       enddo
    enddo

    contains

      function pair_energy (a, b) result (e)
        !
        ! Return the interaction energy between two sites: Pari_energy =
        ! LJ + Coulomb_short + Coulomb_long
        !
        use foreign, only: site
        use units, only: EPSILON0INV, ALPHA
        implicit none
        type (site), intent (in) :: a, b
        real (rk) :: e
        ! *** end of interface ***

        real (rk) :: epsilon, sigma, rab, coul

        ! Combining LJ parameters:
        call pair (a, b, sigma, epsilon)

        rab = norm2 (a % x - b % x)

        ! Coulomb long + short energy happens  to be just 1 / rab, but
        ! in this  way we do not  rely on that equality  should it not
        ! hold for some reason:
        coul = EPSILON0INV * a % charge * b % charge * &
             (coulomb_short (rab, ALPHA) + coulomb_long (rab, ALPHA))

        if (sigma /= 0.0) then
          e = coul + epsilon * lj (rab / sigma)
        else
          e = coul
        endif
      end function pair_energy
  end function self_energy


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
    use units, only: pi
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
    use units, only: pi
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


  elemental function closure1 (method, beta, v, t, dt) result (dc)
    use foreign, only: HNC => CLOSURE_HNC, KH => CLOSURE_KH, PY => CLOSURE_PY
    implicit none
    integer, intent (in) :: method
    real (rk), intent (in) :: beta, v, t, dt
    real (rk) :: dc
    ! *** end of interface ***

    ! FIXME: the other two are not yet implemented:
    select case (method)
    case (HNC)
       dc = closure_hnc1 (beta, v, t, dt)
    case (KH)
       dc = closure_kh1 (beta, v, t, dt)
    ! case (PY)
    !    dc = closure_py1 (beta, v, t, dt)
    case default
       dc = huge (dc)          ! FIXME: cannot abort in pure functions
    end select
  end function closure1


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
  ! correlation t = h - c is denoted by greek "γ" in other sources. We
  ! will avoid greek identifiers utill better times.
  !
  !   c = exp (-βv + t) - 1 - t
  !
  elemental function closure_hnc (beta, v, t) result (c)
    implicit none
    real (rk), intent (in) :: beta, v, t
    real (rk) :: c
    ! *** end of interface ***

    ! c = exp (-beta * v + t) - 1 - t
    c = expm1 (-beta * v + t) - t
  end function closure_hnc


  elemental function closure_hnc1 (beta, v, t, dt) result (dc)
    implicit none
    real (rk), intent (in) :: beta, v, t, dt
    real (rk) :: dc
    ! *** end of interface ***

    ! dc = [exp (-beta * v + t) - 1] dt
    dc = expm1 (-beta * v + t) * dt
  end function closure_hnc1

  !
  ! 2)  Kovalenko-Hirata (KH)  closure.   Same as  HNC in  "depletion"
  ! regions but avoids exponential grows:
  !
  !        / exp (-βv + t) - 1 - t, if -βv + t <= 0
  !   c = <
  !        \ -βv, otherwise
  !
  elemental function closure_kh (beta, v, t) result (c)
    implicit none
    real (rk), intent (in) :: beta, v, t
    real (rk) :: c
    ! *** end of interface ***

    real (rk) :: x

    x = -beta * v + t

    ! For x  <= 0 use  exp(x) - 1,  but do not grow  exponentially for
    ! positive x:
    if (x <= 0.0) then
       c = expm1 (x) - t
    else
       c = -beta * v
    endif
  end function closure_kh


  elemental function closure_kh1 (beta, v, t, dt) result (dc)
    implicit none
    real (rk), intent (in) :: beta, v, t, dt
    real (rk) :: dc
    ! *** end of interface ***

    real (rk) :: x

    x = -beta * v + t

    ! For x  <= 0 use  exp(x) - 1,  but do not grow  exponentially for
    ! positive x:
    if (x <= 0.0) then
       ! dc = [exp (-beta * v + t) - 1] dt
       dc = expm1 (x) * dt
    else
       dc = 0.0
    endif
  end function closure_kh1


  ! 3)  Percus-Yevick  (PY)   closure  relation  between  direct-  and
  ! indirect correlation c and t:
  !
  !   c := exp (-βv) [1 + t] - 1 - t
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
  !    c := exp (-βv + t + B) - 1 - t
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
  ! to compute t =  h - c form c:
  !
  !                -1                -1   2
  !   t =  (1 - ρc)   c - c = (1 - ρc)  ρc
  !
  ! If you scale c by h3 beforehand  or pass rho' = rho * h3 and scale
  ! the result  by h3 in addition,  you will compute  exactly what the
  ! older version of the function did:
  !
  function oz_vv_equation_c_t (rho, C, W) result (T)
    implicit none
    real (rk), intent (in) :: rho
    real (rk), intent (in) :: C(:, :, :) ! (m, m, nrad)
    real (rk), intent (in) :: W(:, :, :) ! (m, m, nrad)
    real (rk) :: T(size (C, 1), size (C, 2), size (C, 3))
    ! *** end of interface ***

    integer :: i

    ! There is  no reason to  handle the 1x1 case  differently, except
    ! clarity.  The MxM branch should be able to handle that case too.
    if (size (C, 1) == 1) then
       ! FIXME: it  is implied here  that W =  1. See comments  on the
       ! value of ω(k) for i == j in omega_fourier().
       do i = 1, size (C, 3)
          T(1, 1, i) = oz_vv_equation_c_t_1x1 (rho, C(1, 1, i))
       enddo
    else
       do i = 1, size (C, 3)
          T(:, :, i) = oz_vv_equation_c_t_MxM (rho, C(:, :, i), W(:, :, i))
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
    call sles (m, T, H)

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


  !
  ! RISM equation, here for h and c:
  !
  !   h   =  ω  * c   * [ω + ρ h  ]
  !    uv     u    uv     v     vv
  !
  ! The term is square brackets  is the solvent property (fixed by the
  ! assumption of the infinite dilution)  and is passed as the solvent
  ! susceptibility
  !
  !   χ   =  ω + ρ h
  !    vv     v     vv
  !
  ! The first of the two functions returns the indirect correlation
  !
  !   t   =  h  -  c
  !    uv     uv    uv
  !
  ! whereas the second returns the direct correlation h = t + c.
  !
  function oz_uv_equation_c_t (c_uvk, w_uuk, x_vvk) result (t_uvk)
    !
    ! Solute-solvent OZ: c -> t
    !
    implicit none
    real (rk), intent (in) :: c_uvk(:, :, :) ! (n, m, nrad)
    real (rk), intent (in) :: w_uuk(:, :, :) ! (n, n, nrad)
    real (rk), intent (in) :: x_vvk(:, :, :) ! (m, m, nrad)
    real (rk) :: t_uvk(size (c_uvk, 1), size (c_uvk, 2), size (c_uvk, 3))
    ! *** end of interface ***

    integer :: k

    !
    ! Many associative  matrix multiplies: NxM =  NxN * NxM  * MxM -
    ! NxM or in u/v terms: uv = uu * uv * vv - uv
    !
    !$omp parallel do
    do k = 1, size (c_uvk, 3)
       associate (uv => c_uvk(:, :, k), uu => w_uuk(:, :, k), vv => x_vvk(:, :, k))
         t_uvk(:, :, k) = matmul (uu, matmul (uv, vv)) - uv
       end associate
    enddo
    !$omp end parallel do
  end function oz_uv_equation_c_t


  function oz_uv_equation_c_h (c_uvk, w_uuk, x_vvk) result (h_uvk)
    !
    ! Solute-solvent OZ: c -> h
    !
    implicit none
    real (rk), intent (in) :: c_uvk(:, :, :) ! (n, m, nrad)
    real (rk), intent (in) :: w_uuk(:, :, :) ! (n, n, nrad)
    real (rk), intent (in) :: x_vvk(:, :, :) ! (m, m, nrad)
    real (rk) :: h_uvk(size (c_uvk, 1), size (c_uvk, 2), size (c_uvk, 3))
    ! *** end of interface ***

    integer :: k

    !$omp parallel do
    do k = 1, size (c_uvk, 3)
       associate (uv => c_uvk(:, :, k), uu => w_uuk(:, :, k), vv => x_vvk(:, :, k))
         h_uvk(:, :, k) = matmul (uu, matmul (uv, vv))
       end associate
    enddo
    !$omp end parallel do
  end function oz_uv_equation_c_h


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
    use units, only: pi
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
    real (rk), intent (in) :: h(:, :, :)  ! (n, m, nrad)
    real (rk), intent (in) :: cs(:, :, :) ! (n, m, nrad)
    real (rk), intent (in) :: cl(:, :, :) ! (n, m, nrad)
    real (rk) :: mu(size (h, 3))
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

    do p = 1, size (h, 3)       ! nrad
       muH = 0.0
       muS = 0.0
       muL = 0.0
       do j = 1, size (h, 2)    ! m
          do i = 1, size (h, 1) ! n
             ! The h² term contributes conditionally. Eventually, only
             ! depletion regions  (h < 0)  contribute (KH).  Threshold
             ! is  supposed to  be  0.0 for  KH functional  (depletion
             ! regions contribute),  anywhere between 1 and  +∞ for GF
             ! functional  (no such  term) and  -∞ for  HNC functional
             ! (contributes unconditionally):
             if (-h(i, j, p) > thresh) then
                muH = muH + rho * h(i, j, p)**2 / 2
             endif

             muS = muS + rho * (-cs(i, j, p) - h(i, j, p) * cs(i, j, p) / 2)
             muL = muL + rho * (             - h(i, j, p) * cl(i, j, p) / 2)
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
    real (rk), intent (in) :: h(:, :, :)  ! (n, m, nrad)
    real (rk), intent (in) :: cs(:, :, :) ! (n, m, nrad)
    real (rk), intent (in) :: cl(:, :, :) ! (n, m, nrad)
    real (rk) :: mu
    ! *** end of interface ***

    real (rk) :: density (size (h, 3))

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
    real (rk), intent (in) :: g(:, :, :) ! (n, m, nrad)
    ! *** end of interface ***

    integer :: i, j, p

    ! This prints a lot of data on tty!
    print *, "# r then g(i, j) for",  size (g, 1), "x", size (g, 2), "pairs"
    do p = 1, size (g, 3)
       write (*, *) r(p), ((g(i, j, p), i = 1, size (g, 1)), j = 1, size (g, 2))
    enddo
  end subroutine gnuplot3


  function clayout (x) result (y)
    implicit none
    real (rk), intent (in) :: x(:, :, :) ! (n, m, nrad)
    real (rk) :: y(size (x, 3), size (x, 1), size (x, 2)) ! (nrad, n, m)
    ! *** end of interface ***

    integer :: na, nb

    na = size (x, 1) * size (x, 2)
    nb = size (x, 3)

    call trans (na, nb, x, y)
  end function clayout


  function flayout (y) result (x)
    implicit none
    real (rk), intent (in) :: y(:, :, :) ! (nrad, n, m)
    real (rk) :: x(size (y, 2), size (y, 3), size (y, 1)) ! (n, m, nrad)
    ! *** end of interface ***

    integer :: na, nb

    na = size (y, 2) * size (y, 3)
    nb = size (y, 1)

    call trans (nb, na, y, x)
  end function flayout


  subroutine trans (na, nb, x, y)
    implicit none
    integer, intent (in) :: na, nb
    real (rk), intent (in) :: x(na, nb)
    real (rk), intent (out) :: y(nb, na)
    ! *** end of interface ***

    y = transpose (x)
  end subroutine trans


end module rism

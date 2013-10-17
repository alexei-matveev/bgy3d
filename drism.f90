module drism
  use kinds, only: rk
  implicit none
  private

  public :: dipole
  public :: center
  public :: dipole_axes
  public :: local_coords
  public :: epsilon_rism
  public :: dipole_density
  public :: dipole_factor
  public :: dipole_correction

contains

  function dipole_density (beta, rho, sites) result (y)
    !
    ! Computes dipole density
    !
    !       4 π    2
    !   y = --- βρμ
    !        9
    !
    use foreign, only: site
    use units, only: pi, EPSILON0INV
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
    use units, only: ANGSTROM
    implicit none
    real (rk), intent (in) :: beta, rho
    real (rk), intent (in) :: eps        ! target epsilon
    type (site), intent (in) :: sites(:) ! (m)
    real (rk), intent (in) :: k(:)       ! (nrad)
    real (rk) :: xk(size (sites), size (sites), size (k)) ! (m, m, nrad)
    ! *** end of inteface ***

    real (rk), parameter :: s = 0.5 * ANGSTROM ! Ref. [KH00a]
    integer :: i, j
    real (rk) :: x(3, size (sites))
    real (rk) :: f(size (k), size (sites)) ! (nrad, m)
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
          xk(i, j, :) = f (:, i) * h(:) * f(:, j)
       enddo
    enddo
  end function dipole_correction

end module drism

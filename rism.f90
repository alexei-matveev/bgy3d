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

  call rism1d (rho=0.5d0, beta=1.0d0, rmax=10.0d0, n=32)

  do i = 1, 10**4
     call test_dst (2**12)
  end do

contains

  subroutine rism1d (rho, beta, rmax, n)
    implicit none
    real(rk), intent(in) :: rho, beta, rmax
    integer, intent(in) :: n
    ! *** end of interface ***

    real(rk), parameter :: pi = 4 * atan (1.0_rk)
    real(rk), parameter :: half = 0.5, alpha = 0.02
    real(rk) :: r(n), k(n), v(n), t(n), c(n), t1(n)
    real(rk) :: diff, factor
    integer :: i, iter
    logical :: converged

    forall (i = 1:n) r(i) = (i - half) * rmax / n
    forall (i = 1:n) k(i) = (i - half) * pi / rmax


    ! LJ potential, sigma=1, epsilon=1:
    v = lj (r)

    ! print *, "r=", r
    ! print *, "v=", v

    iter = 0
    converged = .false.
    t = 0.0
    do while (.not. converged)
       iter = iter + 1

       c = closure_hnc (beta, v, t)

       ! Forward DST, unnormalized:
       c = dst (r * c) / k

       ! OZ equation, convolution integral. FIXME: what is the correct
       ! value of factor here?
       factor = 1.0
       stop "FIX the factor here!"
       t1 = oz_equation_c_t (rho, c * factor)

       ! Inverse DST, renormalized:
       t1 = dst (k * t1) / r / (2 * n)

       diff = maxval (abs (t1 - t))
       converged = (diff < 1.0e-10)

       t = alpha * t1 + (1 - alpha) * t
       print *, "iter=", iter, "diff=", diff
    end do

    ! It was overwritten with c(k):
    c = closure_hnc (beta, v, t)

    do i = 1, n
       print *, t(i), c(i)
    enddo
  end subroutine rism1d

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

    ! t = c / (1.0 - rho * c) - c
    t = rho * (c * c) / (1.0 - rho * c)
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

    f = 4 * sr6 * ( sr6 - 1)
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
       stop "does not match"
    endif
  end subroutine test_dst
end program

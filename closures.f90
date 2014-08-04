module closures
  !
  ! FIXME:  public  functions expect  the  inverse  temperature β  but
  ! internally we operate with potential in units of temperature, βv.
  !
  ! Copyright (c) 2013, 2014 Alexei Matveev
  ! Copyright (c) 2013 Bo Li
  !
  use kinds, only: rk
  implicit none
  private

  ! Elemental functions:
  public :: expm1
  public :: closure, closure1
  public :: closure_rbc

  ! bind (c), used from C only
  public :: rism_closure
  ! *** END OF INTERFACE ***

contains

  subroutine rism_closure (method, beta, n, v, t, c) bind (c)
    !
    ! This one is for using from the C-side. Note that Fortran forbids
    ! aliasing between the output c and the input v and t.
    !
    ! Should be consistent with ./rism.h
    !
    use iso_c_binding, only: c_int, c_double
    implicit none
    integer (c_int), intent (in), value :: method, n
    real (c_double), intent (in), value :: beta
    real (c_double), intent (in) :: v(n), t(n)
    real (c_double), intent (out) :: c(n)
    ! *** end of interface ***

    c = closure (method, beta, v, t)
  end subroutine rism_closure


  elemental function closure (method, beta, v, t) result (c)
    use foreign, only: HNC => CLOSURE_HNC, KH => CLOSURE_KH, PY => CLOSURE_PY, &
         PSE0 => CLOSURE_PSE0, PSE7 => CLOSURE_PSE7
    implicit none
    integer, intent (in) :: method
    real (rk), intent (in) :: beta, v, t
    real (rk) :: c
    ! *** end of interface ***

    integer :: n                ! PSE order

    select case (method)
    case (HNC)
       c = closure_hnc (beta * v, t)
    case (KH)
       c = closure_kh (beta * v, t)
       ! c = closure_pse (1, beta * v, t)
    case (PY)
       c = closure_py (beta * v, t)
    case (PSE0:PSE7)
       n = method - PSE0
       c = closure_pse (n, beta * v, t)
    case default
       c = huge (c)            ! FIXME: cannot abort in pure functions
    end select
  end function closure


  elemental function closure1 (method, beta, v, t, dt) result (dc)
    use foreign, only: HNC => CLOSURE_HNC, KH => CLOSURE_KH, PY => CLOSURE_PY, &
         PSE0 => CLOSURE_PSE0, PSE7 => CLOSURE_PSE7
    implicit none
    integer, intent (in) :: method
    real (rk), intent (in) :: beta, v, t, dt
    real (rk) :: dc
    ! *** end of interface ***

    integer :: n                ! PSE order

    ! FIXME: some are not yet implemented:
    select case (method)
    case (HNC)
       dc = closure_hnc1 (beta * v, t, dt)
    case (KH)
       dc = closure_kh1 (beta * v, t, dt)
    ! case (PY)
    !    dc = closure_py1 (beta * v, t, dt)
    case (PSE0:PSE7)
       n = method - PSE0
       dc = closure_pse1 (n, beta * v, t, dt)
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


  pure function exps1 (n, x) result (y)
    !
    ! Truncated  exponential  series,  an  expansion for  exp(x)  -  1
    ! employed by PSE-n and KH closures:
    !
    !                    k
    !           n       x       x     n      n!   k-1
    !   f(x) = Σ       ----  = ----  Σ      ---- x
    !           k = 1   k!      n!    k = 1  k!
    !
    ! otherwise. The alternative expression  has the advantage that we
    ! can  start with  the coefficient  of the  highest  power without
    ! having to compute it first.  Will silently accept negative n and
    ! for finite input x will return zero just as as for n = 0.
    !
    implicit none
    integer, intent (in) :: n
    real (rk), intent (in) :: x
    real (rk) :: y
    ! *** end of interface ***

    ! FIXME: specialize for  n = 0, 1, 2  where y = 0, x,  x + x**2/2,
    ! maybe?
    real (rk) :: c
    integer :: k

    ! Horner scheme:
    c = 1.0            ! n!/n!, leading coefficient
    y = 0.0
    do k = n, 1, -1    ! iterate n times
       y =  x * y + c
       c = c * k       ! n!/(k-1)!
       !
       ! At this point in each round:
       !
       ! 1) y = 1
       ! 2) y == x + n, i.e. two terms
       ! 3) y == x(x + n) + n(n-1), i.e. three terms
       !
       ! ...
       !         n-1      n-2
       ! n) y = x    +  nx    +  ...  +  n!, all n terms
    enddo
    ! Here c = n!:
    y = (x / c) * y
  end function exps1


  !
  ! 1)  Hypernetted Chain  (HNC)  closure relation  to compute  direct
  ! correlation function c in real space.  The OZ equation (elsewhere)
  ! is  the second relation  between the  two unknowns.   The indirect
  ! correlation t = h - c is denoted by greek "γ" in other sources. We
  ! will avoid greek identifiers utill better times.
  !
  !   c = exp (-βv + t) - 1 - t
  !
  pure function closure_hnc (v, t) result (c)
    implicit none
    real (rk), intent (in) :: v, t
    real (rk) :: c
    ! *** end of interface ***

    ! c = exp (-beta * v + t) - 1 - t
    c = expm1 (-v + t) - t
  end function closure_hnc


  pure function closure_hnc1 (v, t, dt) result (dc)
    implicit none
    real (rk), intent (in) :: v, t, dt
    real (rk) :: dc
    ! *** end of interface ***

    ! dc = [exp (-beta * v + t) - 1] dt
    dc = expm1 (-v + t) * dt
  end function closure_hnc1


  !
  ! 2)  Kovalenko-Hirata (KH)  closure.   Same as  HNC in  "depletion"
  ! regions but avoids exponential grows:
  !
  !        / exp (-βv + t) - 1 - t, if -βv + t <= 0
  !   c = <
  !        \ -βv, otherwise
  !
  pure function closure_kh (v, t) result (c)
    implicit none
    real (rk), intent (in) :: v, t
    real (rk) :: c
    ! *** end of interface ***

    real (rk) :: x

    x = -v + t

    ! For x  <= 0 use  exp(x) - 1,  but do not grow  exponentially for
    ! positive x:
    if (x <= 0.0) then
       c = expm1 (x) - t
    else
       c = -v
    endif
  end function closure_kh


  pure function closure_kh1 (v, t, dt) result (dc)
    implicit none
    real (rk), intent (in) :: v, t, dt
    real (rk) :: dc
    ! *** end of interface ***

    real (rk) :: x

    x = -v + t

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
  pure function closure_py (v, t) result (c)
    implicit none
    real (rk), intent (in) :: v, t
    real (rk) :: c
    ! *** end of interface ***

    c = exp (-v) * (1 + t) - 1 - t
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
       c = closure_hnc_rbc (beta * v, t, expB)
    case default
       c = huge (c)            ! FIXME: cannot abort in pure functions
    end select
  end function closure_rbc

  !
  ! HNC closure with repulsive bridge correction:
  !
  !    c := exp (-βv + t + B) - 1 - t
  !
  pure function closure_hnc_rbc (v, t, expB) result (c)
    implicit none
    real (rk), intent (in) :: v, t, expB
    real (rk) :: c
    ! *** end of interface ***

    c = exp (-v + t) * expB - 1 - t
  end function closure_hnc_rbc


  pure function closure_pse (n, v, t) result (c)
    implicit none
    integer, intent (in) :: n
    real (rk), intent (in) :: v, t ! v in temperature units
    real (rk) :: c
    ! *** end of interface ***

    real (rk) :: x

    x = -v + t

    ! For x < 0 f(x) = exp(x) - 1, but does not grow exponentially for
    ! positive x. E.g. for n = 1, and x >= 0 f(x) = x, so that c = x -
    ! t = -v.  FIXME: should we rather  code c = -v + [f(x) - x] where
    ! the term in brackets is o(x²) for small x?
    c = pse (n, x) - t

  end function closure_pse


  pure function closure_pse1 (n, v, t, dt) result (dc)
    !
    ! Suffix 1 here indicated the differential, not the order.
    !
    implicit none
    integer, intent (in) :: n
    real (rk), intent (in) :: v, t, dt
    real (rk) :: dc
    ! *** end of interface ***

    real (rk) :: x

    x = -v + t

    ! For x <= 0 use [exp(x) -  1] * dt, but do not grow exponentially
    ! for positive x. For n <= 1 and x >= 0 this will return zero:
    dc = pse (n - 1, x) * dt
  end function closure_pse1


  pure function pse (n, x) result (y)
    !
    ! This combination of an  exponential and power series is employed
    ! by PSE-n and KH closures:
    !
    !   f(x) = exp(x) - 1  if x < 0
    !
    ! and
    !                    k
    !           n       x
    !   f(x) = Σ       ----
    !           k = 1   k!
    !
    ! otherwise.  Will silently accept negative n and behaves as for n
    ! = 0.
    !
    implicit none
    integer, intent (in) :: n
    real (rk), intent (in) :: x
    real (rk) :: y
    ! *** end of interface ***

    if (x < 0.0) then
       y = expm1 (x)
    else
       y = exps1 (n, x)
    endif
  end function pse

end module closures

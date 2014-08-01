module closures
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
  ! correlation function c in real space.  The OZ equation (elsewhere)
  ! is  the second relation  between the  two unknowns.   The indirect
  ! correlation t = h - c is denoted by greek "γ" in other sources. We
  ! will avoid greek identifiers utill better times.
  !
  !   c = exp (-βv + t) - 1 - t
  !
  pure function closure_hnc (beta, v, t) result (c)
    implicit none
    real (rk), intent (in) :: beta, v, t
    real (rk) :: c
    ! *** end of interface ***

    ! c = exp (-beta * v + t) - 1 - t
    c = expm1 (-beta * v + t) - t
  end function closure_hnc


  pure function closure_hnc1 (beta, v, t, dt) result (dc)
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
  pure function closure_kh (beta, v, t) result (c)
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


  pure function closure_kh1 (beta, v, t, dt) result (dc)
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
  pure function closure_py (beta, v, t) result (c)
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
  pure function closure_hnc_rbc (beta, v, t, expB) result (c)
    implicit none
    real (rk), intent (in) :: beta, v, t, expB
    real (rk) :: c
    ! *** end of interface ***

    c = exp (-beta * v + t) * expB - 1 - t
  end function closure_hnc_rbc

end module closures

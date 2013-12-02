module snes
  use kinds, only: rk
  implicit none
  private

  public :: func1, func2        ! abstract interfaces
  public :: snes_default, snes_picard
  public :: krylov

  !
  ! *** END OF INTERFACE ***
  !

  !
  ! This defines two  iterator interfaces x -> dx  for use in fixpoint
  ! non-linear problems.  One is  the type of Fortran closures another
  ! is the corresponding type of C closures. Such an iterator function
  ! is supposed to return array of zeros (ok, rather small numbers) at
  ! convergence.
  !
  ! Well, "iterator"  is probably a  bad name.  These  "iterators" are
  ! supoosed  to be pure  functions of  the respective  arguments. The
  ! closure context is supposed to be used read-only (or only to store
  ! intermediate quantities derived from the argument at pre-allocated
  ! positions).  Violation of  these constraints will eventually break
  ! the carefully designed logic of some non-linear solvers.
  !
  abstract interface
     function func1 (x) result (dx)
       import rk
       implicit none
       real (rk), intent (in) :: x(:, :, :)
       real (rk) :: dx(size (x, 1), size (x, 2), size (x, 3))
     end function func1

     function func2 (a, b) result (c)
       import rk
       implicit none
       real (rk), intent (in) :: a(:, :, :), b(:, :, :)
       real (rk) :: c(size (a, 1), size (a, 2), size (a, 3))
     end function func2

     ! b = f(a):
     subroutine arrfunc1 (ctx, n, a, b) bind (c)
       use iso_c_binding, only: c_ptr, c_int, c_double
       implicit none
       type (c_ptr), intent (in), value :: ctx
       integer (c_int), intent (in), value :: n
       real (c_double), intent (in), target :: a(n)
       real (c_double), intent (out), target :: b(n)
     end subroutine arrfunc1

     ! c = f(a, b):
     subroutine arrfunc2 (ctx, n, a, b, c) bind (c)
       use iso_c_binding, only: c_ptr, c_int, c_double
       implicit none
       type (c_ptr), intent (in), value :: ctx
       integer (c_int), intent (in), value :: n
       real (c_double), intent (in), target :: a(n), b(n)
       real (c_double), intent (out), target :: c(n)
     end subroutine arrfunc2
  end interface

  !
  ! This is a concrete function, implemented in C:
  !
  interface
     subroutine rism_snes (ctx, f, df, n, x) bind (c)
       !
       ! void rism_snes (void *ctx, ArrFunc1 f, ArrFunc2 df,
       !                 int n, real x[n])
       !
       ! using procedure  (arrfunc1) as objective and  (when not null)
       ! procedure  (arrfunc2)  as  Jacobian.   In  C-lang  these  are
       ! ArrFunc1 and ArrFunc2.
       !
       ! See ./bgy3d-snes.c
       !
       use iso_c_binding, only: c_ptr, c_int, c_double
       implicit none
       type (c_ptr), intent (in), value :: ctx
       procedure (arrfunc1) :: f
       procedure (arrfunc2) :: df
       integer (c_int), intent (in), value :: n
       real (c_double), intent (inout) :: x(n)
     end subroutine rism_snes

     subroutine rism_krylov (ctx, f, n, b, x) bind (c)
       !
       ! void rism_krylov (void *ctx, ArrFunc1 f,
       !                   int n, real b[n], real x[n])
       !
       ! See ./bgy3d-snes.c
       !
       use iso_c_binding, only: c_ptr, c_int, c_double
       implicit none
       type (c_ptr), intent (in), value :: ctx
       procedure (arrfunc1) :: f
       integer (c_int), intent (in), value :: n
       real (c_double), intent (in) :: b(n)
       real (c_double), intent (inout) :: x(n)
     end subroutine rism_krylov
  end interface

  !
  ! A pointer to to the structure of this type will serve as a closure
  ! context to be passed to the C-world. Upon callback the Fortran sub
  ! implementing a procedure (arrfunc1)  can use the procedure pointer
  ! to perform the actual work.  I wish closures could be made simpler
  ! than that.
  !
  type context
     integer :: shape(3) = -1
     procedure (func1), pointer, nopass :: f => NULL()
     procedure (func2), pointer, nopass :: df => NULL()
  end type context

contains

  subroutine snes_picard (x, f)
    !
    ! Simple Picard iteration
    !
    !   x    = x + α f(x )
    !    n+1    n       n
    !
    implicit none
    real (rk), intent (inout) :: x(:, :, :)
    procedure (func1) :: f  ! (x) -> dx
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


  subroutine snes_default (x, f, df)
    !
    ! Delegates the actual work to Petsc by way of C-func rism_snes().
    !
    use iso_c_binding, only: c_loc
    implicit none
    real (rk), intent (inout) :: x(:, :, :)
    procedure (func1) :: f            ! (x) -> dx
    procedure (func2), optional :: df ! (x, dx) -> df
    ! *** end of interface ***

    type (context), target :: ctx
    procedure (arrfunc2), pointer :: jac

    ctx % shape = shape (x)
    ctx % f => f

    if (present (df)) then
       ctx % df => df
       jac => jacobian
    else
       ctx % df => NULL()
       jac => NULL()
    endif

    call rism_snes (c_loc (ctx), objective, jac, size (x), x)
  end subroutine snes_default


  subroutine objective (ctx, n, x, dx) bind (c)
    !
    ! Implements  procedure  (arrfunc1)  and  will be  passed  to  the
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

    ! We  dont  do any  work  ourselves,  just  extract a  pointer  to
    ! procedure (func1) an let it do the rest:
    call c_f_pointer (ctx, f_ctx)

    block
       real (rk), pointer :: y(:, :, :), dy(:, :, :)

       ! F03 way to create aliases:
       associate (n => f_ctx % shape)
         ! F03 pointer association with contiguous array:
         y(1:n(1), 1:n(2), 1:n(3)) => x
         dy(1:n(1), 1:n(2), 1:n(3)) => dx
       end associate

       ! The warning by GF 4.6 and 4.7 is incorrect:
       ! http://gcc.gnu.org/bugzilla/show_bug.cgi?id=55855
       dy = f_ctx % f (y)
    end block
  end subroutine objective


  subroutine jacobian (ctx, n, x, dx, df) bind (c)
    !
    ! Implements  procedure  (arrfunc2)  and  will be  passed  to  the
    ! rism_snes()   together   with   the   suitable   context.    See
    ! snes_default().
    !
    use iso_c_binding, only: c_f_pointer, c_ptr, c_int
    implicit none
    type (c_ptr), intent (in), value :: ctx
    integer (c_int), intent (in), value :: n
    real (rk), intent (in), target :: x(n), dx(n)
    real (rk), intent (out), target :: df(n)
    ! *** end of interface ***

    type (context), pointer :: f_ctx

    ! We  dont  do any  work  ourselves,  just  extract a  pointer  to
    ! procedure (func2) an let it do the rest:
    call c_f_pointer (ctx, f_ctx)

    block
       real (rk), pointer :: y(:, :, :), dy(:, :, :), dz(:, :, :)

       ! F03 way to create aliases:
       associate (n => f_ctx % shape)
         ! F03 pointer association with contiguous array:
         y(1:n(1), 1:n(2), 1:n(3)) => x
         dy(1:n(1), 1:n(2), 1:n(3)) => dx
         dz(1:n(1), 1:n(2), 1:n(3)) => df
       end associate

       ! The warning by GF 4.6 and 4.7 is incorrect:
       ! http://gcc.gnu.org/bugzilla/show_bug.cgi?id=55855
       dz = f_ctx % df (y, dy)
    end block
  end subroutine jacobian


  function krylov0 (f, b) result (x)
    !
    ! Solves for f(x) = b. May be abused as poor man solver for linear
    ! equations.
    !
    implicit none
    procedure (func1) :: f
    real (rk), intent (in) :: b(:, :, :)
    real (rk) :: x(size (b, 1), size (b, 2), size (b, 3))
    ! *** end of interface ***

    ! FIXME: abusing SNES solver for (presumably) linear problem:
    x = 0.0
    call snes_default (x, fb)
  contains

    function fb (x) result (y)
      implicit none
      real (rk), intent (in) :: x(:, :, :)
      real (rk) :: y(size (x, 1), size (x, 2), size (x, 3))
      ! *** end of interface ***

      y = f(x) - b
    end function fb
  end function krylov0


  function krylov (f, b) result (x)
    !
    ! Solves for f(x) = b assuming linear f(x).
    !
    use iso_c_binding, only: c_loc
    implicit none
    procedure (func1) :: f
    real (rk), intent (in) :: b(:, :, :)
    real (rk) :: x(size (b, 1), size (b, 2), size (b, 3))
    ! *** end of interface ***

    type (context), target :: ctx

    ctx % shape = shape (b)
    ctx % f => f
    ! ctx % df is NULL by default

    x = 0.0
    call rism_krylov (c_loc (ctx), objective, size (b), b, x)
  end function krylov

end module snes

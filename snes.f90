module snes
  use kinds, only: rk
  implicit none
  private

  public :: snes_default, snes_picard

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
     function f_iterator (x) result (dx)
       import rk
       implicit none
       real (rk), intent (in) :: x(:, :, :)
       real (rk) :: dx(size (x, 1), size (x, 2), size (x, 3))
     end function f_iterator

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
  end interface

  !
  ! A pointer to to the structure of this type will serve as a closure
  ! context to be passed to the C-world. Upon callback the Fortran sub
  ! implementing a procedure(arrfunc1) can use the procedure pointer
  ! to perform the actual work.  I wish closures could be made simpler
  ! than that.
  !
  type context
     integer :: shape(3)
     procedure (f_iterator), pointer, nopass :: f
  end type context

contains

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
    use iso_c_binding, only: c_ptr, c_loc, c_null_funptr
    implicit none
    procedure (f_iterator) :: f  ! (x) -> dx
    real (rk), intent (inout) :: x(:, :, :)
    ! *** end of interface ***

    type (context), target :: f_ctx
    type (c_ptr) :: ctx
    procedure (arrfunc2), pointer :: df => NULL()

    f_ctx % f => f
    f_ctx % shape = shape (x)
    ctx = c_loc (f_ctx)

    call rism_snes (ctx, iterator, df, size (x), x)
  end subroutine snes_default


  subroutine iterator (ctx, n, x, dx) bind (c)
    !
    ! Implements  procedure(arrfunc1)  and  will  be passed  to  the
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
    ! procedure(f_iterator) an let it do the rest:
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
  end subroutine iterator
end module snes

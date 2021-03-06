module foreign
  !
  ! Copyright (c) 2013, 2014 Alexei Matveev
  !
  use iso_c_binding, only: c_int, c_double, c_char, c_bool
  use lisp, only: obj
  implicit none
  private

  public :: expm1
  public :: bgy3d_problem_data
  public :: bgy3d_problem_data_print

  public :: bgy3d_getopt_test
  public :: bgy3d_getopt_int
  public :: bgy3d_getopt_real
  public :: bgy3d_getopt_scm
  public :: bgy3d_getopt_string

  public :: comm_rank
  public :: comm_size
  public :: comm_allreduce
  public :: comm_set_parallel_x

  public :: verbosity           ! extern int

  integer (c_int), bind (c) :: verbosity ! see bgy3d.c

  public :: pad                 ! module procedure

  ! Keep this in sync with bgy3d-solutes.h:
  type, public, bind (c) :: site
     character (kind=c_char, len=5) :: name ! atom types. What are they used for?
     real (c_double) :: x(3)                ! coordinates
     real (c_double) :: sigma               ! sigma for LJ
     real (c_double) :: epsilon             ! epsilon for LJ
     real (c_double) :: charge              ! charge
     type (obj) :: obj
  end type site

  ! Keep these in sync with bgy3d.h:
  enum, bind (c)
     enumerator :: CLOSURE_HNC, CLOSURE_KH, CLOSURE_PY
     ! FIXME: more scalable solution?
     enumerator :: CLOSURE_PSE0
     enumerator :: CLOSURE_PSE1
     enumerator :: CLOSURE_PSE2
     enumerator :: CLOSURE_PSE3
     enumerator :: CLOSURE_PSE4
     enumerator :: CLOSURE_PSE5
     enumerator :: CLOSURE_PSE6
     enumerator :: CLOSURE_PSE7
  end enum

  public :: CLOSURE_HNC, CLOSURE_KH, CLOSURE_PY, CLOSURE_PSE0, &
       CLOSURE_PSE1, CLOSURE_PSE2, CLOSURE_PSE3, CLOSURE_PSE4, &
       CLOSURE_PSE5, CLOSURE_PSE6, CLOSURE_PSE7

  type, public, bind (c) :: problem_data
     !
     ! See  bgy3d.h.  The  fields prefixed  with  XX are  not used  in
     ! Fortran code.
     !
     integer (c_int) :: XXN(3)  ! global grid size
     real (c_double) :: XXL(3)  ! box size
     real (c_double) :: XXh(3)  ! mesh width
     real (c_double) :: rmax    ! radial extent (1D)
     integer (c_int) :: nrad    ! number of radial points (1D)
     real (c_double) :: beta    ! inverse temperature, 1/kT
     real (c_double) :: rho     ! solvent density

     ! Other staff  that was retrieved by the  solvers themselves from
     ! the (Petsc) environment:
     real (c_double) :: XXlambda   ! Mixing parameter.
     real (c_double) :: XXdamp     ! Scaling factor.
     integer (c_int) :: XXmax_iter ! Maximal number of iterations.
     real (c_double) :: XXnorm_tol ! Convergence threshold.
     integer (c_int) :: closure    ! enum for HNC, KH, or PY
  end type problem_data

  !
  ! These are concrete functions, implemented in C:
  !
  interface
     ! ProblemData bgy3d_problem_data (void)
     function bgy3d_problem_data () result (pd) bind (c)
       import problem_data
       implicit none
       type (problem_data) :: pd
     end function bgy3d_problem_data

     ! void bgy3d_problem_data_print (const ProblemData *PD)
     subroutine bgy3d_problem_data_print (pd) bind (c)
       import problem_data
       implicit none
       type (problem_data), intent (in) :: pd
     end subroutine bgy3d_problem_data_print

     ! From libm:
     pure function expm1 (x) result (y) bind (c)
       import c_double
       implicit none
       real (c_double), intent (in), value :: x
       real (c_double) :: y
     end function expm1

     !
     ! Keep these consistent with bgy3d-getopt.{h,c}:
     !
     ! bool bgy3d_getopt_test (const char key[]);
     ! bool bgy3d_getopt_int (const char key[], int *val);
     ! bool bgy3d_getopt_real (const char key[], double *val);
     ! bool bgy3d_getopt_string (const char key[], int len, char val[len]);
     !
     function bgy3d_getopt_test (key) result (ok) bind (c)
       import c_bool, c_char
       implicit none
       character (kind=c_char), intent (in) :: key(*) ! 0-terminated
       logical (c_bool) :: ok
     end function bgy3d_getopt_test

     function bgy3d_getopt_int (key, val) result (ok) bind (c)
       import c_bool, c_char, c_int
       implicit none
       character (kind=c_char), intent (in) :: key(*) ! 0-terminated
       integer (c_int), intent (inout) :: val
       logical (c_bool) :: ok
     end function bgy3d_getopt_int

     function bgy3d_getopt_real (key, val) result (ok) bind (c)
       import c_bool, c_char, c_double
       implicit none
       character (kind=c_char), intent (in) :: key(*) ! 0-terminated
       real (c_double), intent (inout) :: val
       logical (c_bool) :: ok
     end function bgy3d_getopt_real

     function bgy3d_getopt_string (key, n, val) result (ok) bind (c)
       import c_bool, c_char, c_int
       implicit none
       character (kind=c_char), intent (in) :: key(*) ! 0-terminated
       integer (c_int), intent (in), value :: n
       character (kind=c_char), intent (inout) :: val(n) ! 0-terminated or unchanged
       logical (c_bool) :: ok
     end function bgy3d_getopt_string

     function bgy3d_getopt_scm (key, val) result (ok) bind (c)
       use lisp, only: obj
       import c_bool, c_char, c_double
       implicit none
       character (kind=c_char), intent (in) :: key(*) ! 0-terminated
       type (obj), intent (inout) :: val
       logical (c_bool) :: ok
     end function bgy3d_getopt_scm

     !
     ! Note the absense of communicator in signatures of communication
     ! primitives.  The actual  communicator  is taken  from a  global
     ! variable.  However,  this variable  is, well, a  *variable* and
     ! can be  flipped between MPI_COMM_WORLD  and MPI_COMM_SELF using
     ! the function  comm_set_parallel_x().  This is  VERY fragile, be
     ! careful.
     !
     ! See ./bgy3d.{h,c}.
     !
     function comm_rank () result (rank) bind (c)
       import c_int
       implicit none
       integer (c_int) :: rank
     end function comm_rank

     function comm_size () result (size) bind (c)
       import c_int
       implicit none
       integer (c_int) :: size
     end function comm_size

     subroutine comm_allreduce (n, x) bind (c)
       import c_int, c_double
       implicit none
       integer (c_int), intent (in), value :: n
       real (c_double), intent (inout) :: x(n)
     end subroutine comm_allreduce

     function comm_set_parallel_x (new) result (old) bind (c)
       import c_bool
       implicit none
       logical (c_bool), intent (in), value :: new
       logical (c_bool) :: old
     end function comm_set_parallel_x
  end interface

contains

  function pad (s) result (t)
    !
    ! Returns a string of the  len (s) with everything after the first
    ! null  char   replaced  by  spaces.  May  be   usefull  to  adapt
    ! null-terminated C-strings.
    !
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

end module foreign

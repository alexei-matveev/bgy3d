module linalg
  use kinds, only: rk
  implicit none
  private

  public :: sles
  public :: eigv
  public :: polyfit

  interface
     subroutine dgesv (n, nrhs, a, lda, ipiv, b, ldb, info)
       use kinds, only: rk
       implicit none
       integer, intent (in) :: n, nrhs, lda, ldb
       real (rk), intent (inout) :: a(lda,*)
       integer, intent (out) :: ipiv(*)
       real (rk), intent (inout) :: b(ldb,*)
       integer, intent (out) :: info
     end subroutine dgesv

     subroutine dsyev (jobz, uplo, n, a, lda, w, work, lwork, info)
       use kinds, only: rk
       implicit none
       character (len=1), intent (in) :: jobz, uplo
       integer, intent (in) :: lda, lwork, n
       integer, intent (out) :: info
       real (rk), intent (inout) :: a(lda, *)
       real (rk), intent (out) :: w(*)
       real (rk), intent (out) :: work(*)
     end subroutine dsyev

     subroutine dgetrf (m, n, a, lda, ipiv, info)
       use kinds, only: rk
       implicit none
       integer, intent (in) :: lda, m, n
       integer, intent (out) :: info
       integer, intent (out) :: ipiv( * )
       real (rk), intent (inout) :: a(lda, *)
     end subroutine dgetrf

     subroutine dgetri (n, a, lda, ipiv, work, lwork, info)
       use kinds, only: rk
       implicit none
       integer, intent (in) :: lda, lwork, n
       integer, intent (out) :: info
       integer, intent (in) :: ipiv(*)
       real (rk), intent (out) :: work(lwork)
       real (rk), intent (inout) :: a(lda, *)
     end subroutine dgetri
  end interface

contains

  subroutine sles (m, a, b)
    !
    ! Solve linear equations  A X = B. As in LAPACK  the matrix A is
    ! destroyed and the result is returned in B.
    !
    implicit none
    integer, intent (in) :: m
    real (rk), intent (inout) :: a(m, m), b(m, m)
    ! *** end of interface ***

    integer :: ipiv(m), info

    ! B will  be overwriten with  the result, A will  be overwritten
    ! with its factorization:
    call dgesv (m, m, a, m, ipiv, b, m, info)

    if (info /= 0) then
       block
          integer :: i
          print *, "a="
          do i = 1, m
             print *, a(i, :)
          enddo
          print *, "b="
          do i = 1, m
             print *, b(i, :)
          enddo
          print *, "info=", info
       end block
       stop "dgesv failed, see tty"
    endif
  end subroutine sles

  subroutine eigv (h, e, v)
    implicit none
    real (rk), intent (in)  :: h(:, :) ! (n, n)
    real (rk), intent (out) :: e(:)    ! (n)
    real (rk), intent (out) :: v(:, :) ! (n, n)
    ! *** end of interface ***

    integer :: i
    real (rk) :: a(size (e), size (e))

    a = h

    call dsyev90 (a, e)

    do i = 1, size (e)
       v(:, i) =  a(:, i)
    enddo
  end subroutine eigv


  subroutine dsyev90 (a, e)
    !
    ! Symmetric eigenvalue problem  solver for real symmetric matrices
    ! (LAPACK DSYEV)
    !
    !   A * V = V * e
    !
    ! A (inout):
    !   (in) A upper triangle
    !   (out) eigenvectors
    !
    ! E (out):
    !    eigenvalues
    !
    implicit none
    real (rk), intent (inout) :: a(:,:)
    real (rk), intent (out) :: e(:)
    ! *** end of interface ***

    integer :: n
    integer :: lwork, info
    real (rk), allocatable :: work(:)

    n = size (e)

    allocate (work(1))

    call dsyev ('V', 'U', n, a, n, e, work, -1, info)
    if (info /= 0) error stop "dsyev: returned nonzero (1)"

    lwork = nint (work(1))

    deallocate (work)

    allocate (work(lwork))

    call dsyev ('V', 'U', n, a, n, e, work, lwork, info)
    if (info /= 0) error stop "dsyev: returned nonzero (2)"

    deallocate (work)
  end subroutine dsyev90


  function polyfit (vx, vy, d) result (a)
    implicit none
    real (rk), intent (in) :: vx(:), vy(:)
    integer, intent (in) :: d
    real (rk) :: a(d + 1)
    ! *** end of interface ***

    integer :: n, lwork

    n = d + 1
    lwork = n

    block
       integer :: i, j
       integer :: info
       real (rk) :: x(size (vx), n)
       real (rk) :: xt(n, size (vx))
       real (rk) :: xtx(n, n)
       real (rk) :: work(lwork)
       integer :: ipiv(n)

       ! Prepare the matrix
       do i = 0, d
          do j = 1, size (vx)
             x(j, i + 1) = vx(j)**i
          end do
       end do

       xt  = transpose (x)
       xtx = matmul (xt, x)

       ! Calls to LAPACK subs DGETRF and DGETRI
       call dgetrf (n, n, xtx, size (xtx, 1), ipiv, info)
       if (info /= 0) then
          error stop "dgetrf failed!"
       end if

       call dgetri (n, xtx, size (xtx, 1), ipiv, work, lwork, info)
       if (info /= 0) then
          error stop "dgetri failed!"
       end if

       a = matmul (matmul (xtx, xt), vy)
    end block
  end function

end module linalg

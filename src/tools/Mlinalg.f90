module Mlinalg
!======================================================================================
   !! This module contains wrapper routines to call BLAS/LAPACK and some useful tools
   !! for linear algebra operations
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef, only: dp, one, zero, iu, large_int
   use Mutils, only: str, stop_error
   implicit none
   include "../formats.h"
!-------------------------------------------------------------------------------------- 
   private
   public :: &
      eig, &
      eigvals, &
      eigvalshe, &
      eigh, &
      inv, &
      solve, &
      eye, &
      det, &
      lstsq, &
      diag, &
      trace, &
      svdvals, &
      svd, &
      util_rotate, &
      util_rotate_cc, &
      blas_zgemm, &
      blas_dgemm, &
      util_zgemm, &
      util_matmul, &
      get_large_size
!-------------------------------------------------------------------------------------- 
  interface eig
     !! eigenvalue/-vector problem for general matrices:
     module procedure deig
     module procedure zeig
  end interface eig

  interface eigh
    !! eigenvalue/-vector problem for real symmetric/complex hermitian matrices:
     module procedure deigh_generalized
     module procedure deigh_simple
     module procedure zeigh_generalized
     module procedure zeigh_simple
  end interface eigh

  interface eigvals
    !! eigenvalues for general matrices:
     module procedure deigvals
     module procedure zeigvals
  end interface eigvals

  interface inv
     !! matrix inversion for real/complex matrices:
     module procedure dinv
     module procedure zinv
  end interface inv

  interface solve
     !! solution to linear systems of equation with real/complex coefficients:
     module procedure dsolve
     module procedure zsolve
  end interface solve

  interface det
     !!  determinants of real/complex square matrices:
     module procedure ddet
     module procedure zdet
  end interface det

  interface lstsq
     !! least square solutions the real/complex systems of equations of possibly non-square shape:
     module procedure dlstsq
     module procedure zlstsq
  end interface lstsq

  interface diag
     !! construction of square matrices from the diagonal elements:
     module procedure ddiag
     module procedure zdiag
  end interface diag

  interface trace
     !! trace of real/complex matrices:
     module procedure dtrace
     module procedure ztrace
  end interface trace

  interface svdvals
     !!  singular values of real/complex matrices:
     module procedure dsvdvals
     module procedure zsvdvals
  end interface svdvals

  interface svd
     !! singular value decomposition of real/complex matrices:
     module procedure dsvd
     module procedure zsvd
  end interface svd

  interface util_matmul
    !! matmul interface, using GEMM for larger matrices
    module procedure dutil_matmul, zutil_matmul
  end interface util_matmul


  interface blas_zgemm
     pure subroutine zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
     !! Interfrace to BLAS routine ZGEMM
       import :: dp
       character,intent(in)      :: TRANSA
       character,intent(in)      :: TRANSB
       integer,intent(in)        :: M
       integer,intent(in)        :: N
       integer,intent(in)        :: K
       complex(dp),intent(in)    :: ALPHA
       complex(dp),intent(in)    :: A(LDA,*)
       integer,intent(in)        :: LDA
       complex(dp),intent(in)    :: B(LDB,*)
       integer,intent(in)        :: LDB
       complex(dp),intent(in)    :: BETA
       complex(dp),intent(inout) :: C(LDC,*)
       integer,intent(in)        :: LDC
     end subroutine zgemm
  end interface blas_zgemm

  interface blas_dgemm
     pure subroutine dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
     !! Interfrace to BLAS routine DGEMM
       import :: dp
       character,intent(in)      :: TRANSA
       character,intent(in)      :: TRANSB
       integer,intent(in)        :: M
       integer,intent(in)        :: N
       integer,intent(in)        :: K
       real(dp),intent(in)       :: ALPHA
       real(dp),intent(in)       :: A(LDA,*)
       integer,intent(in)        :: LDA
       real(dp),intent(in)       :: B(LDB,*)
       integer,intent(in)        :: LDB
       real(dp),intent(in)       :: BETA
       real(dp),intent(inout)    :: C(LDC,*)
       integer,intent(in)        :: LDC
     end subroutine dgemm
  end interface blas_dgemm

contains
!--------------------------------------------------------------------------------------
   pure logical function get_large_size(n)
   !! returns .true. if matrix size n is greater than a fixed number (defined in Mdef.f90)
   !! This is useful for matrix operations: for small matrices, Fortran intrinsic matmul
   !! is typically faster than BLAS
      integer,intent(in) :: n !! vector / matrix size

      get_large_size = n >= large_int

   end function get_large_size
!--------------------------------------------------------------------------------------
   subroutine util_zgemm(a, b, c, transa_opt, transb_opt)
   !! Simplilfied interface for BLAS ZGEMM
   !! Computes C = op(A) * op(B)
      complex(kind=dp), intent(in)            :: a(:, :) !! matrix A
      complex(kind=dp), intent(in)            :: b(:, :) !! matrix B
      complex(kind=dp), intent(out)           :: c(:, :) !! matrix C
      character(len=1), intent(in), optional  :: transa_opt !! operation tag for matrix A
      character(len=1), intent(in), optional  :: transb_opt !! operation tag for matrix B
      integer          :: m, n, k
      character(len=1) :: transa, transb

      transa = 'N'
      transb = 'N'
      if (present(transa_opt)) transa = transa_opt
      if (present(transb_opt)) transb = transb_opt

      ! m ... number of rows in Op(A) and C
      ! n ... number of columns in Op(B) and C
      ! k ... number of columns in Op(A) resp. rows in Op(B)
      m = size(c, 1)
      n = size(c, 2)

      if (transa /= 'N') then
         k = size(a, 1)
      else
         k = size(a, 2)
      end if

      call zgemm(transa, transb, m, n, k, one, a, size(a, 1), b, size(b, 1), zero, c, m)

   end subroutine util_zgemm
!--------------------------------------------------------------------------------------
   subroutine dutil_matmul(a,b,c,alpha,beta,large_size)
   !! computes C = alpha * A * B + beta * C for matrices A, B, C
   !! for small matrices, matmul is used, otherwise DGEMM
      real(kind=dp), intent(in)            :: a(:, :) !! matrix A
      real(kind=dp), intent(in)            :: b(:, :) !! matrix B
      real(kind=dp), intent(inout)         :: c(:, :) !! matrix C
      real(kind=dp), intent(in),optional   :: alpha, beta !! constants alpha,beta
      logical,intent(in),optional          :: large_size !! flag for using BLAS for large matrices
      logical :: large_size_,bt_term
      real(dp) :: al,bt
      integer :: m, n, k

      m = size(c, 1)
      n = size(c, 2)
      k = size(a, 2)

      large_size_ = .false.
      if(present(large_size)) large_size_ = large_size

      al = 1.0_dp
      bt = 0.0_dp
      bt_term = present(beta)
      if(present(alpha)) al = alpha
      if(bt_term) bt = beta

      if(large_size_) then
         call blas_dgemm('n', 'n', m, n, k, al, a, size(a, 1), b, size(b, 1), bt, c, m)
      else
         if(bt_term) then
            c = bt * c + al * matmul(a, b)
         else
            c = al * matmul(a, b)
         end if
      end if

   end subroutine dutil_matmul
!--------------------------------------------------------------------------------------
   subroutine zutil_matmul(a,b,c,alpha,beta,large_size)
   !! computes C = alpha * A * B + beta * C for matrices A, B, C
   !! for small matrices, matmul is used, otherwise DGEMM
      complex(kind=dp), intent(in)            :: a(:, :) !! matrix A
      complex(kind=dp), intent(in)            :: b(:, :) !! matrix B
      complex(kind=dp), intent(inout)         :: c(:, :) !! matrix C
      complex(kind=dp), intent(in),optional   :: alpha, beta !! constants alpha,beta
      logical,intent(in),optional             :: large_size !! flag for using BLAS for large matrices
      logical :: large_size_,bt_term
      complex(dp) :: al,bt
      integer :: m, n, k

      m = size(c, 1)
      n = size(c, 2)
      k = size(a, 2)

      large_size_ = .false.
      if(present(large_size)) large_size_ = large_size

      al = one
      bt = zero
      bt_term = present(beta)
      if(present(alpha)) al = alpha
      if(bt_term) bt = beta

      if(large_size_) then
         call blas_zgemm('n', 'n', m, n, k, al, a, size(a, 1), b, size(b, 1), bt, c, m)
      else
         if(bt_term) then
            c = bt * c + al * matmul(a, b)
         else
            c = al * matmul(a, b)
         end if
      end if

   end subroutine zutil_matmul
!--------------------------------------------------------------------------------------
   function util_rotate(nst,rot,A,large_size) result(A_rot)
   !! rotates a matrix A: C = Q^H A Q with a unitary square matrix Q
      integer,intent(in)     :: nst !! matrix size
      complex(dp),intent(in) :: rot(:,:) !! unitary matrix Q
      complex(dp),intent(in) :: A(:,:) !! matrix to be rotated
      logical,intent(in),optional :: large_size !! flag for using BLAS for large matrices
      logical :: use_zgemm
      complex(dp) :: A_rot(nst,nst)
      complex(dp) :: rot_cc(nst,nst)

      use_zgemm = .false.
      if(present(large_size)) use_zgemm = large_size

      if(use_zgemm) then
        call blas_zgemm('c','n',nst,nst,nst,one,rot,nst,A,nst,zero,rot_cc,nst)
        call blas_zgemm('n','n',nst,nst,nst,one,rot_cc,nst,rot,nst,zero,A_rot,nst)
      else
        rot_cc = conjg(transpose(rot))
        A_rot = matmul(rot_cc(1:nst,1:nst), A(1:nst,1:nst))
        A_rot = matmul(A_rot(1:nst,1:nst), rot(1:nst,1:nst))
      end if

   end function util_rotate
!--------------------------------------------------------------------------------------
   function util_rotate_cc(nst,rot,A,large_size) result(A_rot)
   !! rotates a matrix A: C = Q A Q^H with a unitary square matrix Q
      integer,intent(in)     :: nst !! matrix size
      complex(dp),intent(in) :: rot(:,:) !! unitary matrix Q
      complex(dp),intent(in) :: A(:,:) !! matrix to be rotated
      logical,intent(in),optional :: large_size !! flag for using BLAS for large matrices
      logical :: use_zgemm
      complex(dp) :: A_rot(nst,nst)
      complex(dp) :: rot_cc(nst,nst)

      use_zgemm = .false.
      if(present(large_size)) use_zgemm = large_size

      if(use_zgemm) then
        call blas_zgemm('n','n',nst,nst,nst,one,rot,nst,A,nst,zero,rot_cc,nst)
        call blas_zgemm('n','c',nst,nst,nst,one,rot_cc,nst,rot,nst,zero,A_rot,nst)
      else
        rot_cc = conjg(transpose(rot))
        A_rot = matmul(rot(1:nst,1:nst), A(1:nst,1:nst))
        A_rot = matmul(A_rot(1:nst,1:nst), rot_cc(1:nst,1:nst))
      end if

   end function util_rotate_cc
!--------------------------------------------------------------------------------------
  ! TODO: add optional switch for left or right eigenvectors in deig() and zeig()?
  subroutine deig(A, lam, c)
   !! solves real nonsymmetric eigenvalue problem
    real(dp), intent(in) :: A(:, :)  !! matrix for eigenvalue compuation
    complex(dp), intent(out) :: lam(:)  !! eigenvalues: A c = lam c
    complex(dp), intent(out) :: c(:, :)  !! eigenvectors: A c = lam c; c(i,j) = ith component of jth vec.
    ! LAPACK variables for DGEEV:
    real(dp), allocatable ::  At(:,:), vl(:,: ), vr(:,:), wi(:), work(:), wr(:)
    integer :: info, lda, ldvl, ldvr, lwork, n, i

    lda = size(A(:,1))
    n = size(A(1,:))
    call assert_shape(A, [n, n], "solve", "A")
    call assert_shape(c, [n, n], "solve", "c")
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    allocate(At(lda,n), wr(n), wi(n), vl(ldvl,n), vr(ldvr,n), work(lwork))
    At = A

    call dgeev('N', 'V', n, At, lda, wr, wi, vl, ldvl, vr, ldvr, &
         work, lwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "dgeev returned info = ", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "the QR algorithm failed to compute all the"
          write(error_unit,fmt900) "eigenvalues, and no eigenvectors have been computed;"
          write(error_unit,fmt900) "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          write(error_unit,fmt900) "have converged."
       end if
       call stop_error('eig: dgeev error')
    end if

    lam = wr + iu*wi
    ! as DGEEV has a rather complicated way of returning the eigenvectors,
    ! it is necessary to build the complex array of eigenvectors from
    ! two real arrays:
    do i = 1,n
       if(wi(i) > 0.0) then  ! first of two conjugate eigenvalues
          c(:, i) = vr(:, i) + iu*vr(:, i+1)
       elseif(wi(i) < 0.0_dp) then  ! second of two conjugate eigenvalues
          c(:, i) = vr(:, i-1) - iu*vr(:, i)
       else
          c(:, i) = vr(:, i)
       end if
    end do
  end subroutine deig

  subroutine zeig(A, lam, c)
    complex(dp), intent(in) :: A(:, :)  ! matrix to solve eigenproblem for
    complex(dp), intent(out) :: lam(:)  ! eigenvalues: A c = lam c
    complex(dp), intent(out) :: c(:,:)  ! eigenvectors: A c = lam c; c(i,j) = ith component of jth vec.
    ! LAPACK variables:
    integer :: info, lda, ldvl, ldvr, lwork, n, lrwork
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: vl(:,:), vr(:,:), work(:)

    lda = size(A(:,1))
    n = size(A(1,:))
    call assert_shape(A, [n, n], "solve", "A")
    call assert_shape(c, [n, n], "solve", "c")
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    lrwork = 2*n
    allocate(vl(ldvl,n), vr(ldvr,n), work(lwork), rwork(lrwork))
    c = A
    call zgeev('N', 'V', n, c, lda, lam, vl, ldvl, vr, ldvr, work, &
         lwork, rwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "zgeev returned info = ", info
       if(info < 0) then
          write(error_unit,fmt900) "the ",-info, "-th argument had an illegal value."
       else
          write(error_unit,fmt900) "the QR algorithm failed to compute all the"
          write(error_unit,fmt900) "eigenvalues, and no eigenvectors have been computed;"
          write(error_unit,fmt900) "elements and ", info+1, ":", n, " of W contain eigenvalues which have"
          write(error_unit,fmt900) "converged."
       end if
       call stop_error('eig: zgeev error')
    end if
    c = vr
  end subroutine zeig

  subroutine EigValsHE(A,lam)
    complex(dp), intent(in) :: A(:, :)  ! matrix to solve eigenproblem for
    real(dp), intent(out) :: lam(:)  ! eigenvalues: A c = lam c
    ! LAPACK variables:
    integer :: info, lda, lwork, n, lrwork
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:), c(:,:)

    lda = size(A(:,1))
    n = size(A(1,:))
    call assert_shape(A, [n, n], "solve", "A")
    lwork = 2*n  ! TODO: can this size be optimized? query first?
    lrwork = 3*n
    allocate(work(lwork), rwork(lrwork), c(n,n))
    c = A
    call zheev('N', 'U', n, c, lda, lam, work, lwork, rwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "zheev returned info = ", info
       if(info < 0) then
          write(error_unit,fmt900) "the ",-info, "-th argument had an illegal value."
       else
          write(error_unit,fmt900) "the QR algorithm failed to compute all the"
          write(error_unit,fmt900) "eigenvalues, and no eigenvectors have been computed;"
          write(error_unit,fmt900) "elements and ", info+1, ":", n, " of W contain eigenvalues which have"
          write(error_unit,fmt900) "converged."
       end if
       call stop_error('eigvals: zheev error')
    end if

    deallocate(work, rwork, c)
    
  end subroutine EigValsHE

  function deigvals(A) result(lam)
    real(dp), intent(in) :: A(:, :)  ! matrix for eigenvalue compuation
    complex(dp), allocatable :: lam(:)  ! eigenvalues: A c = lam c
    ! LAPACK variables for DGEEV:
    real(dp), allocatable ::  At(:,:), vl(:,: ), vr(:,:), wi(:), work(:), wr(:)
    integer :: info, lda, ldvl, ldvr, lwork, n

    lda = size(A(:,1))
    n = size(A(1,:))
    call assert_shape(A, [n, n], "solve", "A")
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    allocate(At(lda,n), wr(n), wi(n), vl(ldvl,n), vr(ldvr,n), work(lwork), lam(n))
    At = A

    call dgeev('N', 'N', n, At, lda, wr, wi, vl, ldvl, vr, ldvr, &
         work, lwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "dgeev returned info = ", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "the QR algorithm failed to compute all the"
          write(error_unit,fmt900) "eigenvalues, and no eigenvectors have been computed;"
          write(error_unit,fmt900) "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          write(error_unit,fmt900) "have converged."
       end if
       call stop_error('eigvals: dgeev error')
    end if

    lam = wr + iu*wi
  end function deigvals

  function zeigvals(A) result(lam)
    complex(dp), intent(in) :: A(:, :)  ! matrix to solve eigenproblem for
    complex(dp), allocatable :: lam(:)  ! eigenvalues: A c = lam c
    ! LAPACK variables:
    integer :: info, lda, ldvl, ldvr, lwork, n, lrwork
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: At(:,:), vl(:,:), vr(:,:), work(:)

    lda = size(A(:,1))
    n = size(A(1,:))
    call assert_shape(A, [n, n], "solve", "A")
    ldvl = n
    ldvr = n
    lwork = 8*n  ! TODO: can this size be optimized? query first?
    lrwork = 2*n
    allocate(At(lda,n), vl(ldvl,n), vr(ldvr,n), work(lwork), rwork(lrwork), lam(n))
    At = A
    call zgeev('N', 'N', n, At, lda, lam, vl, ldvl, vr, ldvr, work, &
         lwork, rwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "zgeev returned info = ", info
       if(info < 0) then
          write(error_unit,fmt900) "the ",-info, "-th argument had an illegal value."
       else
          write(error_unit,fmt900) "the QR algorithm failed to compute all the"
          write(error_unit,fmt900) "eigenvalues, and no eigenvectors have been computed;"
          write(error_unit,fmt900) "elements and ", info+1, ":", n, " of W contain eigenvalues which have"
          write(error_unit,fmt900) "converged."
       end if
       call stop_error('eig: zgeev error')
    end if
  end function zeigvals

  subroutine deigh_generalized(Am, Bm, lam, c)
    ! solves generalized eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric, Bm symmetric positive definite.
    ! Only the lower triangular part of Am and Bm is used.
    real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    real(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(dp), intent(inout) :: lam(:)   ! eigenvalues: Am c = lam Bm c
    real(dp), intent(inout) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
    integer :: n
    ! lapack variables
    integer :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: Bmt(:,:), work(:)

    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(Bm, [n, n], "eigh", "B")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 1 + 6*n + 2*n**2
    liwork = 3 + 5*n
    allocate(Bmt(n,n), work(lwork), iwork(liwork))
    c = Am; Bmt = Bm  ! Bmt temporaries overwritten by dsygvd
    call dsygvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       write(error_unit,fmt900) "dsygvd returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else if (info <= n) then
          write(error_unit,fmt900) "the algorithm failed to compute an eigenvalue while working"
          write(error_unit,fmt900) "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          write(error_unit,fmt900) "through", mod(info, n+1)
       else
          write(error_unit,fmt900) "The leading minor of order ", info-n, &
               "of B is not positive definite. The factorization of B could ", &
               "not be completed and no eigenvalues or eigenvectors were computed."
       end if
       call stop_error('eigh: dsygvd error')
    end if
  end subroutine deigh_generalized

  subroutine deigh_simple(Am, lam, c)
    ! solves eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric
    ! Only the lower triangular part of Am is used.
    real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam c
    real(dp), intent(inout) :: lam(:)   ! eigenvalues: Am c = lam c
    real(dp), intent(inout) :: c(:,:)   ! eigenvectors: Am c = lam c; c(i,j) = ith component of jth vec.
    integer :: n
    ! lapack variables
    integer :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: work(:)

    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 1 + 6*n + 2*n**2
    liwork = 3 + 5*n
    allocate(work(lwork), iwork(liwork))
    c = Am
    call dsyevd('V','L',n,c,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       write(error_unit,fmt900) "dsyevd returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "the algorithm failed to compute an eigenvalue while working"
          write(error_unit,fmt900) "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          write(error_unit,fmt900) "through", mod(info, n+1)
       end if
       call stop_error('eigh: dsyevd error')
    end if
  end subroutine deigh_simple

  subroutine zeigh_generalized(Am, Bm, lam, c)
    ! solves generalized eigen value problem for all eigenvalues and eigenvectors
    ! Am must by hermitian, Bm hermitian positive definite.
    ! Only the lower triangular part of Am and Bm is used.
    complex(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    complex(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(dp), intent(inout) :: lam(:)      ! eigenvalues: Am c = lam Bm c
    complex(dp), intent(inout) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
    ! lapack variables
    integer :: info, liwork, lrwork, lwork, n
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: Bmt(:,:), work(:)

    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(Bm, [n, n], "eigh", "Bm")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 2*n + n**2
    lrwork = 1 + 5*N + 2*n**2
    liwork = 3 + 5*n
    allocate(Bmt(n,n), work(lwork), rwork(lrwork), iwork(liwork))
    c = Am; Bmt = Bm  ! Bmt temporary overwritten by zhegvd
    call zhegvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,rwork,lrwork,iwork,liwork,info)
    if (info /= 0) then
       write(error_unit,fmt900) "zhegvd returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else if (info <= n) then
          write(error_unit,fmt900) "the algorithm failed to compute an eigenvalue while working"
          write(error_unit,fmt900) "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          write(error_unit,fmt900) "through", mod(info, n+1)
       else
          write(error_unit,fmt900) "The leading minor of order ", info-n, &
               "of B is not positive definite. The factorization of B could ", &
               "not be completed and no eigenvalues or eigenvectors were computed."
       end if
       call stop_error('eigh: zhegvd error')
    end if
  end subroutine zeigh_generalized

  subroutine zeigh_simple(Am, lam, c)
    ! solves eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric
    ! Only the lower triangular part of Am is used.
    complex(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam c
    complex(dp), intent(inout) :: c(:,:)   ! eigenvectors: Am c = lam c; c(i,j) = ith component of jth vec.
    ! LAPACK variables:
    integer :: info, lda, liwork, lrwork, lwork, n
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:)

    ! use LAPACK's zheevd routine
    n = size(Am, 1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(c, [n, n], "eigh", "c")
    lda = max(1, n)
    lwork = 2*n + n**2
    lrwork = 1 + 5*n + 2*n**2
    liwork = 3 + 5*n
    allocate(work(lwork), rwork(lrwork), iwork(liwork))
    c = Am
    call zheevd("V", "L", n, c, lda, lam, work, lwork, rwork, lrwork, &
         iwork, liwork, info)
    if (info /= 0) then
       write(error_unit,fmt900) "zheevd returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "the algorithm failed to compute an eigenvalue while working"
          write(error_unit,fmt900) "through the submatrix lying in rows and columns through"
          write(error_unit,fmt900) info/(n+1), " through ", mod(info, n+1)
       end if
       call stop_error('eigh: zheevd error')
    end if
  end subroutine zeigh_simple

  function dinv(Am) result(Bm)
    real(dp), intent(in) :: Am(:,:)  ! matrix to be inverted
    real(dp) :: Bm(size(Am, 1), size(Am, 2))   ! Bm = inv(Am)
    real(dp), allocatable :: Amt(:,:), work(:)  ! temporary work arrays

    ! LAPACK variables:
    integer ::  info, lda, n, lwork, nb
    integer, allocatable :: ipiv(:)
    integer, external :: ilaenv

    ! use LAPACK's dgetrf and dgetri
    n = size(Am(1, :))
    call assert_shape(Am, [n, n], "inv", "Am")
    lda = n
    nb = ilaenv(1, 'DGETRI', "UN", n, -1, -1, -1)  ! TODO: check UN param
    lwork = n*nb
    if (nb < 1) nb = max(1, n)
    allocate(Amt(n,n), work(lwork), ipiv(n))
    Amt = Am
    call dgetrf(n, n, Amt, lda, ipiv, info)
    if(info /= 0) then
       write(error_unit,fmt900) "dgetrf returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "U(", info, ",", info, ") is exactly zero; The factorization"
          write(error_unit,fmt900) "has been completed, but the factor U is exactly"
          write(error_unit,fmt900) "singular, and division by zero will occur if it is used"
          write(error_unit,fmt900) "to solve a system of equations."
       end if
       call stop_error('inv: dgetrf error')
    end if

    call dgetri(n, Amt, n, ipiv, work, lwork, info)
    if (info /= 0) then
       write(error_unit,fmt900) "dgetri returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "U(", info, ",", info, ") is exactly zero; the matrix is"
          write(error_unit,fmt900) "singular and its inverse could not be computed."
       end if
       call stop_error('inv: dgetri error')
    end if
    Bm = Amt

  end function dinv

  function zinv(Am) result(Bm)
    ! Inverts the general complex matrix Am
    complex(dp), intent(in) :: Am(:,:)   ! Matrix to be inverted
    complex(dp) :: Bm(size(Am, 1), size(Am, 2))   ! Bm = inv(Am)
    integer :: n, nb
    ! lapack variables
    integer :: lwork, info
    complex(dp), allocatable:: Amt(:,:), work(:)
    integer, allocatable:: ipiv(:)
    integer, external :: ilaenv

    n = size(Am, 1)
    call assert_shape(Am, [n, n], "inv", "Am")
    nb = ilaenv(1, 'ZGETRI', "UN", n, -1, -1, -1)  ! TODO: check UN param
    if (nb < 1) nb = max(1, n)
    lwork = n*nb
    allocate(Amt(n,n), ipiv(n), work(lwork))
    Amt = Am
    call zgetrf(n, n, Amt, n, ipiv, info)
    if (info /= 0) then
       write(error_unit,fmt900) "zgetrf returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "U(", info, ",", info, ") is exactly zero; The factorization"
          write(error_unit,fmt900) "has been completed, but the factor U is exactly"
          write(error_unit,fmt900) "singular, and division by zero will occur if it is used"
          write(error_unit,fmt900) "to solve a system of equations."
       end if
       call stop_error('inv: zgetrf error')
    end if
    call zgetri(n, Amt, n, ipiv, work, lwork, info)
    if (info /= 0) then
       write(error_unit,fmt900) "zgetri returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "U(", info, ",", info, ") is exactly zero; the matrix is"
          write(error_unit,fmt900) "singular and its inverse could not be computed."
       end if
       call stop_error('inv: zgetri error')
    end if
    Bm = Amt
  end function zinv

  function dsolve(A, b) result(x)
    ! solves a system of equations A x = b with one right hand side
    real(dp), intent(in) :: A(:,:)  ! coefficient matrix A
    real(dp), intent(in) :: b(:)  ! right-hand-side A x = b
    real(dp), allocatable :: x(:)
    ! LAPACK variables:
    real(dp), allocatable :: At(:,:), bt(:,:)
    integer :: n, info, lda
    integer, allocatable :: ipiv(:)

    n = size(A(1,:))
    lda = size(A(:, 1))  ! TODO: remove lda (which is = n!)
    call assert_shape(A, [n, n], "solve", "A")
    allocate(At(lda,n), bt(n,1), ipiv(n), x(n))
    At = A
    bt(:,1) = b(:)
    call dgesv(n, 1, At, lda, ipiv, bt, n, info)
    if(info /= 0) then
       write(error_unit,fmt900) "dgesv returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "U(", info, ",", info, ") is exactly zero; The factorization"
          write(error_unit,fmt900) "has been completed, but the factor U is exactly"
          write(error_unit,fmt900) "singular, so the solution could not be computed."
       end if
       call stop_error('inv: dgesv error')
    endif
    x = bt(:,1)
  end function dsolve

  function zsolve(A, b) result(x)
    ! solves a system of equations A x = b with one right hand side
    complex(dp), intent(in) :: A(:,:)  ! coefficient matrix A
    complex(dp), intent(in) :: b(:)  ! right-hand-side A x = b
    complex(dp), allocatable :: x(:)
    ! LAPACK variables:
    complex(dp), allocatable :: At(:,:), bt(:,:)
    integer :: n, info, lda
    integer, allocatable :: ipiv(:)

    n = size(A(1,:))
    lda = size(A(:, 1))  ! TODO: remove lda here, too
    call assert_shape(A, [n, n], "solve", "A")
    allocate(At(lda,n), bt(n,1), ipiv(n), x(n))
    At = A
    bt(:,1) = b(:)
    call zgesv(n, 1, At, lda, ipiv, bt, n, info)
    if(info /= 0) then
       write(error_unit,fmt900) "zgesv returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "U(", info, ",", info, ") is exactly zero; The factorization"
          write(error_unit,fmt900) "has been completed, but the factor U is exactly"
          write(error_unit,fmt900) "singular, so the solution could not be computed."
       end if
       call stop_error('inv: zgesv error')
    endif
    x = bt(:,1)
  end function zsolve

  pure function eye(n) result(A)
    ! Returns the identity matrix of size n x n and type real.
    integer, intent(in) :: n
    real(dp) :: A(n, n)
    integer :: i

    A = 0
    do i = 1, n
       A(i, i) = 1
    end do
  end function eye

  function ddet(A) result(x)
    ! compute the determinant of a real matrix using an LU factorization
    real(dp), intent(in) :: A(:, :)
    real(dp) :: x
    integer :: i
    ! LAPACK variables:
    integer :: info, n
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: At(:,:)

    n = size(A(1,:))
    call assert_shape(A, [n, n], "det", "A")
    allocate(At(n,n), ipiv(n))
    At = A
    call dgetrf(n, n, At, n, ipiv, info)
    if(info /= 0) then
       write(error_unit,fmt900) "dgetrf returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "U(", info, ",", info, ") is exactly zero; The factorization"
          write(error_unit,fmt900) "has been completed, but the factor U is exactly"
          write(error_unit,fmt900) "singular, and division by zero will occur if it is used"
          write(error_unit,fmt900) "to solve a system of equations."
       end if
       call stop_error('det: dgetrf error')
    end if

    ! At now contains the LU of the factorization A = PLU
    ! as L has unit diagonal entries, the determinant can be computed
    ! from the product of U's diagonal entries. Additional sign changes
    ! stemming from the permutations P have to be taken into account as well.
    x = 1.0_dp
    do i = 1,n
       if(ipiv(i) /= i) then  ! additional sign change
          x = -x*At(i,i)
       else
          x = x*At(i,i)
       endif
    end do
  end function ddet

  function zdet(A) result(x)
    ! compute the determinant of a real matrix using an LU factorization
    complex(dp), intent(in) :: A(:, :)
    complex(dp) :: x
    integer :: i
    ! LAPACK variables:
    integer :: info, n
    integer, allocatable :: ipiv(:)
    complex(dp), allocatable :: At(:,:)

    n = size(A(1,:))
    call assert_shape(A, [n, n], "det", "A")
    allocate(At(n,n), ipiv(n))
    At = A
    call zgetrf(n, n, At, n, ipiv, info)
    if(info /= 0) then
       write(error_unit,fmt900) "zgetrf returned info =", info
       if (info < 0) then
          write(error_unit,fmt900) "the", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "U(", info, ",", info, ") is exactly zero; The factorization"
          write(error_unit,fmt900) "has been completed, but the factor U is exactly"
          write(error_unit,fmt900) "singular, and division by zero will occur if it is used"
          write(error_unit,fmt900) "to solve a system of equations."
       end if
       call stop_error('inv: zgetrf error')
    end if

    ! for details on the computation, compare the comment in ddet().
    x = 1.0_dp + 0*iu
    do i = 1,n
       if(ipiv(i) /= i) then  ! additional sign change
          x = -x*At(i,i)
       else
          x = x*At(i,i)
       endif
    end do
  end function zdet

  function dlstsq(A, b) result(x)
    ! compute least square solution to A x = b for real A, b
    real(dp), intent(in) :: A(:,:), b(:)
    real(dp), allocatable :: x(:)
    ! LAPACK variables:
    integer :: info, ldb, lwork, m, n, rank
    real(dp) :: rcond
    real(dp), allocatable :: work(:), At(:,:), Bt(:,:)
    integer, allocatable :: jpvt(:)

    m = size(A(:,1)) ! = lda
    n = size(A(1,:))
    ldb = size(b)
    allocate(x(n), At(m,n), Bt(ldb,1), jpvt(n), work(1))
    call dgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         -1, info)  ! query optimal workspace size
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))  ! allocate with ideal size
    rcond = 0.0_dp
    jpvt(:) = 0
    Bt(:,1) = b(:)  ! only one right-hand side
    At(:,:) = A(:,:)
    call dgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         lwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "dgelsy returned info = ", info
       write(error_unit,fmt900) "the ", -info, "-th argument had an illegal value"
       call stop_error('lstsq: dgelsy error')
    endif
    x(:) = Bt(1:n,1)
  end function dlstsq

  function zlstsq(A, b) result(x)
    ! compute least square solution to A x = b for complex A, b
    complex(dp), intent(in) :: A(:,:), b(:)
    complex(dp), allocatable :: x(:)
    ! LAPACK variables:
    integer :: info, ldb, lwork, m, n, rank
    real(dp) :: rcond
    complex(dp), allocatable :: At(:,:), Bt(:,:), work(:)
    real(dp), allocatable :: rwork(:)
    integer, allocatable :: jpvt(:)

    m = size(A(:,1)) ! = lda
    n = size(A(1,:))
    ldb = size(b)
    allocate(x(n), At(m,n), Bt(ldb,1), jpvt(n), work(1), rwork(2*n))
    call zgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         -1, rwork, info)  ! query optimal workspace size
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))  ! allocate with ideal size
    rcond = 0.0_dp
    jpvt(:) = 0
    Bt(:,1) = b(:)  ! only one right-hand side
    At(:,:) = A(:,:)
    call zgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         lwork, rwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "zgelsy returned info = ", info
       write(error_unit,fmt900) "the ", -info, "-th argument had an illegal value"
       call stop_error('lstsq: zgelsy error')
    endif
    x(:) = Bt(1:n,1)
  end function zlstsq

  ! TODO: can assumed types help in Xdiag() and Xtrace()?
  ! TODO: add optional axis parameter in both xdiag() functions
  function ddiag(x) result(A)
    ! construct real matrix from diagonal elements
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: A(:,:)
    integer :: i, n

    n = size(x)
    allocate(A(n,n))
    A(:,:) = 0.0_dp
    forall(i=1:n) A(i,i) = x(i)
  end function ddiag

  function zdiag(x) result(A)
    ! construct complex matrix from diagonal elements
    complex(dp), intent(in) :: x(:)
    complex(dp), allocatable :: A(:,:)
    integer :: i, n

    n = size(x)
    allocate(A(n,n))
    A(:,:) = 0*iu
    forall(i=1:n) A(i,i) = x(i)
  end function zdiag

  ! TODO: add optional axis parameter in both xtrace() functions
  pure function dtrace(A) result(t)
    ! return trace along the main diagonal
    real(dp), intent(in) :: A(:,:)
    real(dp) :: t
    integer :: i

    t = 0.0_dp
    do i = 1,minval(shape(A))
       t = t + A(i,i)
    end do
  end function dtrace

  pure function ztrace(A) result(t)
    ! return trace along the main diagonal
    complex(dp), intent(in) :: A(:,:)
    complex(dp) :: t
    integer :: i

    t = 0*iu
    do i = 1,minval(shape(A))
       t = t + A(i,i)
    end do
  end function ztrace

  function dsvdvals(A) result(s)
    ! compute singular values s_i of a real m x n matrix A
    real(dp), intent(in) :: A(:,:)
    real(dp), allocatable :: s(:)
    ! LAPACK related:
    integer :: info, lwork, m, n
    real(dp), allocatable :: work(:), At(:,:)
    real(dp) :: u(1,1), vt(1,1)  ! not used if only s is to be computed

    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    allocate(At(m,n), s(min(m,n)))
    At(:,:) = A(:, :)  ! A is overwritten in dgesvd

    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call dgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, -1, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    call dgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, lwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "dgesvd returned info = ", info
       if(info < 0) then
          write(error_unit,fmt900) "the ", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "DBDSQR did not converge, there are ", info
          write(error_unit,fmt900) "superdiagonals of an intermediate bidiagonal form B"
          write(error_unit,fmt900) "did not converge to zero. See the description of WORK"
          write(error_unit,fmt900) "in DGESVD's man page for details."
       endif
       call stop_error('svdvals: dgesvd error')
    endif
  end function dsvdvals

  function zsvdvals(A) result(s)
    ! compute singular values s_i of a real m x n matrix A
    complex(dp), intent(in) :: A(:,:)
    real(dp), allocatable :: s(:)
    ! LAPACK related:
    integer :: info, lwork, m, n, lrwork
    complex(dp), allocatable :: work(:), At(:,:)
    real(dp), allocatable :: rwork(:)
    complex(dp) :: u(1,1), vt(1,1)  ! not used if only s is to be computed

    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    lrwork = 5*min(m,n)
    allocate(At(m,n), s(min(m,n)), rwork(lrwork))
    At(:,:) = A(:,:)  ! A is overwritten in zgesvd!

    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call zgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, -1, rwork, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    call zgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, lwork, rwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "zgesvd returned info = ", info
       if(info < 0) then
          write(error_unit,fmt900) "the ", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "ZBDSQR did not converge, there are ", info
          write(error_unit,fmt900) "superdiagonals of an intermediate bidiagonal form B"
          write(error_unit,fmt900) "did not converge to zero. See the description of RWORK"
          write(error_unit,fmt900) "in ZGESVD's man page for details."
       endif
       call stop_error('svdvals: zgesvd error')
    endif
  end function zsvdvals

  subroutine dsvd(A, s, U, Vtransp)
    ! compute the singular value decomposition A = U sigma Vtransp of a
    ! real m x n matrix A
    ! U is m x m
    ! Vtransp is n x n
    ! s has size min(m, n) --> sigma matrix is (n x m) with sigma_ii = s_i
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(out) :: s(:), U(:,:), Vtransp(:,:)
    ! LAPACK related:
    integer :: info, lwork, m, n, ldu
    real(dp), allocatable :: work(:), At(:,:)

    ! TODO: check shapes here and in other routines?

    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    ldu = m
    allocate(At(m,n))
    At(:,:) = A(:,:)  ! use a temporary as dgesvd destroys its input

    call assert_shape(U, [m, m], "svd", "U")
    call assert_shape(Vtransp, [n, n], "svd", "Vtransp")

    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call dgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, -1, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    call dgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, lwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "dgesvd returned info = ", info
       if(info < 0) then
          write(error_unit,fmt900) "the ", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "DBDSQR did not converge, there are ", info
          write(error_unit,fmt900) "superdiagonals of an intermediate bidiagonal form B"
          write(error_unit,fmt900) "did not converge to zero. See the description of WORK"
          write(error_unit,fmt900) "in DGESVD's man page for details."
       endif
       call stop_error('svd: dgesvd error')
    endif
  end subroutine dsvd

  subroutine zsvd(A, s, U, Vtransp)
    ! compute the singular value decomposition A = U sigma V^H of a
    ! complex m x m matrix A
    ! U is m x min(m, n)
    ! Vtransp is n x n
    ! sigma is m x n with with sigma_ii = s_i
    ! note that this routine returns V^H, not V!
    complex(dp), intent(in) :: A(:,:)
    real(dp), intent(out) :: s(:)
    complex(dp), intent(out) :: U(:,:), Vtransp(:,:)
    ! LAPACK related:
    integer :: info, lwork, m, n, ldu, lrwork
    real(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: work(:), At(:,:)

    ! TODO: check shapes here and in other routines?

    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    ldu = m
    lrwork = 5*min(m,n)
    allocate(rwork(lrwork), At(m,n))
    At(:,:) = A(:,:)  ! use a temporary as zgesvd destroys its input

    call assert_shape(U, [m, m], "svd", "U")
    call assert_shape(Vtransp, [n, n], "svd", "Vtransp")

    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call zgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, -1,&
         rwork, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))

    call zgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, &
         lwork, rwork, info)
    if(info /= 0) then
       write(error_unit,fmt900) "zgesvd returned info = ", info
       if(info < 0) then
          write(error_unit,fmt900) "the ", -info, "-th argument had an illegal value"
       else
          write(error_unit,fmt900) "ZBDSQR did not converge, there are ", info
          write(error_unit,fmt900) "superdiagonals of an intermediate bidiagonal form B"
          write(error_unit,fmt900) "did not converge to zero. See the description of WORK"
          write(error_unit,fmt900) "in DGESVD's man page for details."
       endif
       call stop_error('svd: zgesvd error')
    endif
  end subroutine zsvd

end module Mlinalg

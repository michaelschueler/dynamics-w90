module scitools_idrs
!! This is the "Induced Dimension Reduction", IDR(s) (for s=4). IDR(s) is a robust and efficient short recurrence
!! Krylov subspace method for solving large nonsymmetric systems of linear equations. It is described in
!! [Peter Sonneveld and Martin B. van Gijzen, SIAM J. Sci. Comput. 31, 1035 (2008)]. We have adapted the code
!! released by M. B. van Gizjen [http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html].
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use scitools_debug
   use scitools_def,only: dp,iu,one,zero
   use scitools_utils,only: stop_error
   use scitools_linalg,only: util_axpy, util_copy
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: idrs, didrs, zidrs
!--------------------------------------------------------------------------------------
   interface idrs
      module procedure didrs, zidrs
   end interface idrs
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
function didrs( b, s, &
  preconditioner, matrixvector, &
  ddotprod, &
  tolerance, maximum_iterations, variant, print_eachstep, &
  flag, relres, iterations, &
  x0, U0, omega, resvec, H)

  implicit none

  ! Required input parameters:
  real(dp), intent(in) :: b(:, :)  ! system rhs
  integer, intent(in) :: s       ! s parameter
  ! Solution:
  real(dp) :: didrs (size(b,1), size(b,2))
  ! Optional input parameters:
  real(dp), optional, intent(in) :: tolerance
  integer, optional, intent(in) :: maximum_iterations
  character(len=8), optional, intent(in) :: variant
  logical, optional, intent(in) :: print_eachstep
  ! Optional output parameters:
  integer, optional, intent(out) :: flag
  real(dp), optional, intent(out) :: relres
  integer, optional, intent(out) :: iterations
  ! Optional input arrays:
  real(dp), optional, intent(in) :: x0(:, :)
  real(dp), optional, intent(in) :: U0(:, :, :)
  real(dp), optional, intent(in) :: omega(:)
  ! Optional output arrays:
  real(dp), optional, intent(out) :: resvec(:)
  real(dp), optional, intent(out) :: H(:, :)

  interface
    function preconditioner(v)
      import :: dp
      real(dp), dimension(:,:), intent(in)     :: v
      real(dp), dimension(size(v,1),size(v,2)) :: preconditioner
    end function preconditioner
    function matrixvector(v)
      import :: dp
      real(dp), intent(in)       :: v(:,:)
      real(dp)                   :: matrixvector(size(v,1),size(v,2))
    end function matrixvector
    real(dp) function ddotprod(a, b)
      import :: dp
      real(dp), intent(in) :: a(:), b(:)
      real(dp), allocatable :: apsi(:, :), bpsi(:, :)
    end function ddotprod
  end interface

  ! Local arrays:
  real(dp), allocatable           :: P(:,:,:)
  real(dp), allocatable          :: R0(:,:)
  real(dp)                       :: x(size(b,1),size(b,2))
  real(dp)                       :: G(size(b,1),size(b,2),s)
  real(dp)                       :: U(size(b,1),size(b,2),s)
  real(dp)                       :: r(size(b,1),size(b,2))
  real(dp)                       :: v(size(b,1),size(b,2))
  real(dp)                       :: t(size(b,1),size(b,2))
  real(dp)                       :: M(s,s), f(s), mu(s)
  real(dp)                       :: alpha(s), beta(s), gamma(s)

  real(dp)                       :: om, tr
  real(dp)                        :: nr, nt, rho, kappa

  ! Declarations:
  integer               :: n                  !< dimension of the system
  integer               :: nrhs               !< Number of RHS-vectors
  integer               :: maxit              !< maximum number of iterations
  integer               :: method             !< which IDR(s) variant?
  real(dp)                 :: tol                !< actual tolerance
  integer               :: info               !< convergence indicator
  logical               :: print_flag         !< print each iteration
  logical               :: out_flag           !< store flag
  logical               :: out_relres         !< store relres
  logical               :: out_iterations     !< store number of iterations
  logical               :: inispace           !< initial search space
  logical               :: user_omega         !< user defined omega present
  integer               :: n_omega            !< number of user defined omegas
  logical               :: out_resvec         !< store residual norms
  logical               :: out_H              !< store iteration parameters in H
  integer               :: nritz              !< Number of wanted ritz values

  integer               :: iter               !< number of iterations
  integer               :: ii                 !< inner iterations index
  integer               :: jj                 !< G-space index
  real(dp)                 :: normb, normr, tolb !< for tolerance check
  integer               :: i,j,k,l            !< loop counters

  ! Problem size:
  n    = size(b,1)
  ! Number of right-hand-side vectors:
  nrhs = size(b,2)

  ! Check optional input parameters:
  if (present(tolerance)) then
    if (tolerance < 0) then
      call stop_error('The tolerance parameter in idrs routine must be non-negative')
    end if
    tol = tolerance
  else
    tol = 1.0e-6_dp
  end if

  maxit = min(2 * n, 1000)
  if (present(maximum_iterations)) maxit = maximum_iterations

  method = 1 ! biortho
  if (present(variant)) then
    if (variant == 'minsync') then
      method = 2
    else if (variant == 'bicgstab') then
      method = 3
    end if
  endif

  print_flag = .false.
  if(present(print_eachstep)) print_flag = print_eachstep

  ! Initialize the output variables
  out_flag       = present(flag)
  if (out_flag)       flag = -1
  out_relres     = present(relres)
  if (out_relres)      relres = 1.0_dp
  out_iterations = present(iterations)
  if (out_iterations) iterations = 0

  ! Check optional input arrays:
  x = 0.0_dp
  if (present(x0)) x = x0

  U = 0.0_dp
  inispace =  present(U0)
  if (inispace) U = U0

  user_omega = present(omega)
  if (user_omega) then
    n_omega = size(omega)
  end if

  ! Check output arrays
  out_resvec     = present(resvec)
  if (out_resvec) then
    if (maxit+1 > size(resvec)) then
      call stop_error('idrs: Length of vector with residual norms too small, should be maxit+1')
    end if
  end if

  out_H = present(H)
  if (out_H) then
    nritz = size(H,1)-1
    if (size(H,2) /= nritz) then
      call stop_error('Second dimension of H incompatible, with first')
    end if
    H = 0.0_dp
  end if

  ! compute initial residual, set absolute tolerance
  normb = dfrob_norm(b)
  ! Originally, the residue was normalized with ||b||, but to be consistent with other methods,
  ! we remove this feature
  !tolb = tol * normb
  tolb = tol
  r = b - matrixvector(x)
  normr = dfrob_norm(r)

  if (out_resvec) resvec(1) = normr

  ! check if the initial solution is not already a solution within the prescribed
  ! tolerance
  if (normr <= tolb) then
    if (out_iterations) iterations = 0
    if (out_flag)       flag  = 0
    ! Originally, the residue was normalized with ||b||, but to be consistent with other methods,
    ! we remove this feature
    !if (out_relres)     relres = normr/normb
    if (out_relres)     relres = normr
    return
  end if

  ! Define P and kappa (depending on the method)
  if (method == 1) then
    allocate( P(n,nrhs,s))
    call RANDOM_SEED
    call RANDOM_NUMBER(P)
    do j = 1,s
      do k = 1,j-1
        P(:,:,j) = P(:,:,j) - dtrace_dot( P(:,:,k),P(:,:,j)) * P(:,:,k)
      end do
      P(:,:,j) = P(:,:,j)/dfrob_norm(P(:,:,j))
    end do
    kappa = 0.7_dp
  else if (method == 2) then
    ! P is piecewise constant, minimum residual for omega
    kappa = 0.0_dp
  else if (method == 3) then
    !if (s /= 1) stop "s=1 is required for variant bicgstab"
    if(s /= 1) call stop_error("s=1 is required for variant bicgstab")
    allocate(R0(n,nrhs))
    R0 = r
    kappa = 0.0_dp
  end if

  ! Initialize local variables:
  M = 0.0_dp
  om = 1.0_dp
  iter = 0
  info = -1
  jj = 0
  ii = 0

  ! This concludes the initialisation phase

  ! Main iteration loop, build G-spaces:

  do while (info < 0)  ! start of iteration loop

!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Generate s vectors in G_j
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! New right-hand side for small system:
    f = dp_dot(P, R0, r, s)

    do k = 1, s

      ! Update inner iteration counter
      ii = ii + 1

      ! Compute new v
      v = r
      if (jj > 0) then

        ! Solve small system (Note: M is lower triangular) and make v orthogonal to P:
        do i = k,s
          gamma(i) = f(i)
          do j = k, i-1
            gamma(i) = gamma(i) - M(i,j)*gamma(j)
          end do
          gamma(i) = gamma(i)/M(i,i)
          v = v - gamma(i)*G(:,:,i)
        end do

        ! Compute new U(:,:,k)
        t = om * preconditioner(v)
        do i = k,s
          t = t + gamma(i)*U(:,:,i)
        end do
        U(:,:,k) = t

        ! Compute Hessenberg matrix?
        if (out_H .and. ii <= nritz) then
          H(ii-s:ii-k,ii)   = -gamma(k:s)/beta(k:s)
        end if

      else if (.not. inispace) then

        ! Updates for the first s iterations (in G_0):
        U(:, :, k) = preconditioner(v)

      end if

      ! Compute new G(:,:,k), G(:,:,k) is in space G_j
      G(:, :, k) = matrixvector(U(:, :, k))

      ! Bi-Orthogonalise the new basis vectors:
      mu = dp_dot(P, R0, G(:,:,k), s)
      do i = 1, k-1
        alpha(i) = mu(i)
        do j = 1, i-1
          alpha(i) = alpha(i) - M(i,j) * alpha(j)
        end do
        alpha(i) = alpha(i) / M(i,i)
        G(:,:,k) = G(:,:,k) - G(:,:,i) * alpha(i)
        U(:,:,k) = U(:,:,k) - U(:,:,i) * alpha(i)
        mu(k:s)  = mu(k:s)  - M(k:s,i) * alpha(i)
      end do
      M(k:s,k) = mu(k:s)

      ! Compute Hessenberg matrix?
      if (out_H .and. ii <= nritz .and. k > 1) then
        H(ii-k+1:ii-1,ii) =  alpha(1:k-1) / beta(1:k-1)
      end if

      ! Break down?
      if (abs(M(k,k)) <= tiny(tol)) then
        info = 3
        exit
      end if

      ! Make r orthogonal to p_i, i = 1..k, update solution and residual
      beta(k) = f(k) / M(k,k)
      r = r - beta(k) * G(:,:,k)
      x = x + beta(k) * U(:,:,k)

      ! New f = P(prime) *r (first k  components are zero)
      if (k < s) then
        f(k+1:s)   = f(k+1:s) - beta(k)*M(k+1:s,k)
      end if

      ! Compute Hessenberg matrix?
      if (out_H .and. ii <= nritz) then
        H(ii,ii) = 1.0_dp/beta(k)
        l = max(1,ii-s)
        H(l+1:ii+1,ii) = (H(l+1:ii+1,ii) - H(l:ii,ii))
        H(l:ii+1,ii)   = H(l:ii+1,ii)/om
      end if

      ! Check for convergence
      normr = dfrob_norm(r)

      if(print_flag) then
         write(output_unit,'("idrs: it = ",i6," res = ",E12.5)') iter, normr
      end if

      iter = iter + 1
      if (out_resvec) resvec(iter + 1) = normr
      if (normr < tolb) then
        info = 0
        exit
      else if (iter == maxit) then
        info = 1
        exit
      end if

    end do ! Now we have computed s+1 vectors in G_j
    if (info >= 0) then
      exit
    end if

!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute first residual in G_j+1
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Update G-space counter
    jj = jj + 1

    ! Compute first residual in G_j+1
    ! Note: r is already perpendicular to P so v = r

    ! Preconditioning:
    v = preconditioner(r)
    t = matrixvector(v)


    ! Computation of a new omega
    if (user_omega) then
      i = mod(jj,n_omega)
      if (i == 0) i = n_omega
      om = omega(i)
    else if (abs(kappa) <= 1.0e-14_dp) then

      ! Minimal residual (same as in Bi-CGSTAB):
      om = dtrace_dot(t,r)/dtrace_dot(t,t)
    else

      ! 'Maintaining the convergence':
      nr = dfrob_norm(r)
      nt = dfrob_norm(t)
      tr = dtrace_dot(t, r)
      rho = abs(tr / (nt * nr))
      om = tr / (nt * nt)
      if (rho < kappa) then
        om = om * kappa / rho
      end if
    end if
    if (abs(om) <= epsilon(tol)) then
      info = 3
      exit
    end if

    ! Update solution and residual
    r = r - om*t
    x = x + om*v

    ! Check for convergence
    normr = dfrob_norm(r)
    iter = iter + 1
    if (out_resvec) resvec(iter + 1) = normr
    if (normr < tolb) then
      info = 0
    else if (iter == maxit) then
      info = 1
    end if

  end do ! end of while loop

  ! Set output parameters
  r = b - matrixvector(x)
  normr = dfrob_norm(r)

  if (info == 0 .and. normr > tolb) info = 2
  if (out_iterations) iterations = iter
  ! Originally, the residue was normalized with ||b||, but to be consistent with other methods,
  ! we remove this feature
  !if (out_relres)     relres=normr/normb
  if (out_relres)     relres=normr
  if (out_flag)       flag = info

  didrs = x

contains

  !> Trace inner product of complex matrices
  function dtrace_dot(v, w)
    real(dp), intent(in)      :: v(:,:), w(:,:)
    real(dp)                  :: dtrace_dot
    integer k
    dtrace_dot = 0.0_dp
    do k = 1, size(v, 2)
      dtrace_dot = dtrace_dot + ddotprod(v(:, k), w(:, k))
    end do
  end function dtrace_dot

  !> P inner product of complex matrices
  function dp_dot(P, R0, w, s)
    real(dp),    allocatable, intent(in) :: P(:,:,:)
    real(dp), allocatable, intent(in) :: R0(:,:)
    real(dp), intent(in)              :: w(:,:)
    integer, intent(in)             :: s
    real(dp)                          :: dp_dot(s)

    real(dp)                          :: v(s)
    integer                         :: j, k, N, low(s), up(s), step, nrhs

    if (allocated(P)) then
      ! Biortho: P has orthogonal random numbers
      do i = 1, s
        v(i) = 0.0_dp
        do k = 1, size(w, 2)
          v(i) = v(i) + ddotprod(P(:, k, i), w(:, k))
        end do
      end do
    else if (allocated(R0)) then
      ! BiCGSTAB: shadow vector equal to initial residual
      v(1) = 0.0_Dp
      do k = 1, size(w, 2)
        v(1) = v(1) + ddotprod(R0(:, k), w(:, k))
      end do
    else
      ! Minsync: P is piecewise constant
      ! WARNING: the integrals are done here in a peculiar way, not consistent with the
      ! definition of the dot product. So this is probably not working.
      N    = size(w,1)
      nrhs = size(w,2)
      step = N / s
      low(1) = 1
      do i = 1, s-1
        low(i+1) = i * step + 1
        up(i) = i * step
      end do
      up(s) = N

      do i = 1, s
        v(i)  = 0.0_dp
        do j = 1, nrhs
          v(i) = v(i) + sum(w(low(i):up(i),j))
        end do
      end do
    end if

    dp_dot = v
  end function dp_dot

  !> Frobenius norm of complex matrix
  function dfrob_norm(v)
    real(dp), intent(in) :: v(:,:)
    real(dp)             :: dfrob_norm
    integer :: k
    dfrob_norm = 0.0_dp
    do k = 1, size(v, 2)
      dfrob_norm = dfrob_norm + ddotprod(v(:, k), v(:, k))
    end do
    dfrob_norm = sqrt(dfrob_norm)
  end function dfrob_norm

end function didrs




function zidrs( b, s, &
  preconditioner, matrixvector, &
  ddotprod, zdotprod, &
  tolerance, maximum_iterations, variant, print_eachstep, &
  flag, relres, iterations,&
  x0, U0, omega, resvec, H)

  implicit none

  ! Required input parameters:
  complex(dp), intent(in) :: b(:, :)  ! system rhs
  integer, intent(in) :: s       ! s parameter
  ! Solution:
  complex(dp) :: zidrs (size(b,1), size(b,2))
  ! Optional input parameters:
  real(dp), optional, intent(in) :: tolerance
  integer, optional, intent(in) :: maximum_iterations
  character(len=8), optional, intent(in) :: variant
  logical, optional, intent(in) :: print_eachstep
  ! Optional output parameters:
  integer, optional, intent(out) :: flag
  real(dp), optional, intent(out) :: relres
  integer, optional, intent(out) :: iterations
  ! Optional input arrays:
  complex(dp), optional, intent(in) :: x0(:, :)
  complex(dp), optional, intent(in) :: U0(:, :, :)
  complex(dp), optional, intent(in) :: omega(:)
  ! Optional output arrays:
  real(dp), optional, intent(out) :: resvec(:)
  complex(dp), optional, intent(out) :: H(:, :)

  interface
    function preconditioner(v)
      import :: dp
      complex(dp), dimension(:,:), intent(in)     :: v
      complex(dp), dimension(size(v,1),size(v,2)) :: preconditioner
    end function preconditioner
    function matrixvector(v)
      import :: dp
      complex(dp), intent(in)       :: v(:,:)
      complex(dp)                   :: matrixvector(size(v,1),size(v,2))
    end function matrixvector
    real(dp) function ddotprod(a, b)
    import :: dp
      real(dp), intent(in) :: a(:), b(:)
      real(dp), allocatable :: apsi(:, :), bpsi(:, :)
    end function ddotprod
    complex(dp) function zdotprod(a, b)
       import :: dp
      complex(dp), intent(in) :: a(:), b(:)
      complex(dp), allocatable :: apsi(:, :), bpsi(:, :)
    end function zdotprod
  end interface

  ! Local arrays:
  real(dp), allocatable           :: P(:,:,:)
  complex(dp), allocatable          :: R0(:,:)
  complex(dp)                       :: x(size(b,1),size(b,2))
  complex(dp)                       :: G(size(b,1),size(b,2),s)
  complex(dp)                       :: U(size(b,1),size(b,2),s)
  complex(dp)                       :: r(size(b,1),size(b,2))
  complex(dp)                       :: v(size(b,1),size(b,2))
  complex(dp)                       :: t(size(b,1),size(b,2))
  complex(dp)                       :: M(s,s), f(s), mu(s)
  complex(dp)                       :: alpha(s), beta(s), gamma(s)

  complex(dp)                       :: om, tr
  real(dp)                        :: nr, nt, rho, kappa

  ! Declarations:
  integer               :: n                  !< dimension of the system
  integer               :: nrhs               !< Number of RHS-vectors
  integer               :: maxit              !< maximum number of iterations
  integer               :: method             !< which IDR(s) variant?
  real(dp)                 :: tol                !< actual tolerance
  integer               :: info               !< convergence indicator
  logical               :: print_flag         !< print each iteration
  logical               :: out_flag           !< store flag
  logical               :: out_relres         !< store relres
  logical               :: out_iterations     !< store number of iterations
  logical               :: inispace           !< initial search space
  logical               :: user_omega         !< user defined omega present
  integer               :: n_omega            !< number of user defined omegas
  logical               :: out_resvec         !< store residual norms
  logical               :: out_H              !< store iteration parameters in H
  integer               :: nritz              !< Number of wanted ritz values

  integer               :: iter               !< number of iterations
  integer               :: ii                 !< inner iterations index
  integer               :: jj                 !< G-space index
  real(dp)                 :: normb, normr, tolb !< for tolerance check
  integer               :: i,j,k,l            !< loop counters

  ! Problem size:
  n    = size(b,1)
  ! Number of right-hand-side vectors:
  nrhs = size(b,2)

  ! Check optional input parameters:
  if (present(tolerance)) then
    if (tolerance < 0) then
      call stop_error('The tolerance parameter in idrs routine must be non-negative')
    end if
    tol = tolerance
  else
    tol = 1.0e-6_dp
  end if

  maxit = min(2 * n, 1000)
  if (present(maximum_iterations)) maxit = maximum_iterations

  method = 1 ! biortho
  if (present(variant)) then
    if (variant == 'minsync') then
      method = 2
    else if (variant == 'bicgstab') then
      method = 3
    end if
  endif

  print_flag = .false.
  if(present(print_eachstep)) print_flag = print_eachstep

  ! Initialize the output variables
  out_flag       = present(flag)
  if (out_flag)       flag = -1
  out_relres     = present(relres)
  if (out_relres)      relres = one
  out_iterations = present(iterations)
  if (out_iterations) iterations = 0

  ! Check optional input arrays:
  x = zero
  if (present(x0)) x = x0

  U = zero
  inispace =  present(U0)
  if (inispace) U = U0

  user_omega = present(omega)
  if (user_omega) then
    n_omega = size(omega)
  end if

  ! Check output arrays
  out_resvec     = present(resvec)
  if (out_resvec) then
    if (maxit+1 > size(resvec)) then
      call stop_error('idrs: Length of vector with residual norms too small, should be maxit+1')
    end if
  end if

  out_H = present(H)
  if (out_H) then
    nritz = size(H,1)-1
    if (size(H,2) /= nritz) then
      call stop_error('Second dimension of H incompatible, with first')
    end if
    H = zero
  end if

  ! compute initial residual, set absolute tolerance
  normb = zfrob_norm(b)
  ! Originally, the residue was normalized with ||b||, but to be consistent with other methods,
  ! we remove this feature
  !tolb = tol * normb
  tolb = tol
  r = b - matrixvector(x)
  normr = zfrob_norm(r)
  if (out_resvec) resvec(1) = normr

  ! check if the initial solution is not already a solution within the prescribed
  ! tolerance
  if (normr <= tolb) then
    if (out_iterations) iterations = 0
    if (out_flag)       flag  = 0
    ! Originally, the residue was normalized with ||b||, but to be consistent with other methods,
    ! we remove this feature
    !if (out_relres)     relres = normr/normb
    if (out_relres)     relres = normr
    return
  end if

  ! Define P and kappa (depending on the method)
  if (method == 1) then
    allocate( P(n,nrhs,s))
    call RANDOM_SEED
    call RANDOM_NUMBER(P)
    do j = 1,s
      do k = 1,j-1
        P(:,:,j) = P(:,:,j) - dtrace_dot( P(:,:,k),P(:,:,j)) * P(:,:,k)
      end do
      P(:,:,j) = P(:,:,j)/dfrob_norm(P(:,:,j))
    end do
    kappa = 0.7_dp
  else if (method == 2) then
    ! P is piecewise constant, minimum residual for omega
    kappa = zero
  else if (method == 3) then
    !if (s /= 1) stop "s=1 is required for variant bicgstab"
    if(s /=1 ) call stop_error("s=1 is required for variant bicgstab")
    allocate(R0(n,nrhs))
    R0 = r
    kappa = zero
  end if

  ! Initialize local variables:
  M = zero
  om = ONE
  iter = 0
  info = -1
  jj = 0
  ii = 0

  ! This concludes the initialisation phase

  ! Main iteration loop, build G-spaces:

  do while (info < 0)  ! start of iteration loop

!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Generate s vectors in G_j
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! New right-hand side for small system:
    f = zp_dot(P, R0, r, s)

    do k = 1, s

      ! Update inner iteration counter
      ii = ii + 1

      ! Compute new v
      v = r
      if (jj > 0) then

        ! Solve small system (Note: M is lower triangular) and make v orthogonal to P:
        do i = k,s
          gamma(i) = f(i)
          do j = k, i-1
            gamma(i) = gamma(i) - M(i,j)*gamma(j)
          end do
          gamma(i) = gamma(i)/M(i,i)
          v = v - gamma(i)*G(:,:,i)
        end do

        ! Compute new U(:,:,k)
        t = om * preconditioner(v)
        do i = k,s
          t = t + gamma(i)*U(:,:,i)
        end do
        U(:,:,k) = t

        ! Compute Hessenberg matrix?
        if (out_H .and. ii <= nritz) then
          H(ii-s:ii-k,ii)   = -gamma(k:s)/beta(k:s)
        end if

      else if (.not. inispace) then

        ! Updates for the first s iterations (in G_0):
        U(:, :, k) = preconditioner(v)

      end if

      ! Compute new G(:,:,k), G(:,:,k) is in space G_j
      G(:, :, k) = matrixvector(U(:, :, k))

      ! Bi-Orthogonalise the new basis vectors:
      mu = zp_dot(P, R0, G(:,:,k), s)
      do i = 1, k-1
        alpha(i) = mu(i)
        do j = 1, i-1
          alpha(i) = alpha(i) - M(i,j) * alpha(j)
        end do
        alpha(i) = alpha(i) / M(i,i)
        G(:,:,k) = G(:,:,k) - G(:,:,i) * alpha(i)
        U(:,:,k) = U(:,:,k) - U(:,:,i) * alpha(i)
        mu(k:s)  = mu(k:s)  - M(k:s,i) * alpha(i)
      end do
      M(k:s,k) = mu(k:s)

      ! Compute Hessenberg matrix?
      if (out_H .and. ii <= nritz .and. k > 1) then
        H(ii-k+1:ii-1,ii) =  alpha(1:k-1) / beta(1:k-1)
      end if

      ! Break down?
      if (abs(M(k,k)) <= tiny(tol)) then
        info = 3
        exit
      end if

      ! Make r orthogonal to p_i, i = 1..k, update solution and residual
      beta(k) = f(k) / M(k,k)
      r = r - beta(k) * G(:,:,k)
      x = x + beta(k) * U(:,:,k)

      ! New f = P(prime) *r (first k  components are zero)
      if (k < s) then
        f(k+1:s)   = f(k+1:s) - beta(k)*M(k+1:s,k)
      end if

      ! Compute Hessenberg matrix?
      if (out_H .and. ii <= nritz) then
        H(ii,ii) = one/beta(k)
        l = max(1,ii-s)
        H(l+1:ii+1,ii) = (H(l+1:ii+1,ii) - H(l:ii,ii))
        H(l:ii+1,ii)   = H(l:ii+1,ii)/om
      end if

      ! Check for convergence
      normr = zfrob_norm(r)
      iter = iter + 1
      if (out_resvec) resvec(iter + 1) = normr
      if (normr < tolb) then
        info = 0
        exit
      else if (iter == maxit) then
        info = 1
        exit
      end if

    end do ! Now we have computed s+1 vectors in G_j
    if (info >= 0) then
      exit
    end if

!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute first residual in G_j+1
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Update G-space counter
    jj = jj + 1

    ! Compute first residual in G_j+1
    ! Note: r is already perpendicular to P so v = r

    ! Preconditioning:
    v = preconditioner(r)
    t = matrixvector(v)


    ! Computation of a new omega
    if (user_omega) then
      i = mod(jj,n_omega)
      if (i == 0) i = n_omega
      om = omega(i)
    else if (abs(kappa) <= 1.0e-14_dp) then

      ! Minimal residual (same as in Bi-CGSTAB):
      om = ztrace_dot(t,r)/ztrace_dot(t,t)
    else

      ! 'Maintaining the convergence':
      nr = zfrob_norm(r)
      nt = zfrob_norm(t)
      tr = ztrace_dot(t, r)
      rho = abs(tr / (nt * nr))
      om = tr / (nt * nt)
      if (rho < kappa) then
        om = om * kappa / rho
      end if
    end if
    if (abs(om) <= epsilon(tol)) then
      info = 3
      exit
    end if

    ! Update solution and residual
    r = r - om*t
    x = x + om*v

    ! Check for convergence
    normr = zfrob_norm(r)
    iter = iter + 1

   if(print_flag) then
      write(output_unit,'("idrs: it = ",i6," res = ",E12.5)') iter, normr
   end if

    if (out_resvec) resvec(iter + 1) = normr
    if (normr < tolb) then
      info = 0
    else if (iter == maxit) then
      info = 1
    end if

  end do ! end of while loop

  ! Set output parameters
  r = b - matrixvector(x)
  normr = zfrob_norm(r)

  if (info == 0 .and. normr > tolb) info = 2
  if (out_iterations) iterations = iter
  ! Originally, the residue was normalized with ||b||, but to be consistent with other methods,
  ! we remove this feature
  !if (out_relres)     relres=normr/normb
  if (out_relres)     relres=normr
  if (out_flag)       flag = info

  zidrs = x

contains

  !> Trace inner product of complex matrices
  function dtrace_dot(v, w)
    real(dp), intent(in)      :: v(:,:), w(:,:)
    real(dp)                  :: dtrace_dot
    integer k
    dtrace_dot = 0.0_dp
    do k = 1, size(v, 2)
      dtrace_dot = dtrace_dot + ddotprod(v(:, k), w(:, k))
    end do
  end function dtrace_dot

  !> Trace inner product of complex matrices
  function ztrace_dot(v, w)
    complex(dp), intent(in)      :: v(:,:), w(:,:)
    complex(dp)                  :: ztrace_dot
    integer :: k
    ztrace_dot = zero
    do k = 1, size(v, 2)
      ztrace_dot = ztrace_dot + zdotprod(v(:, k), w(:, k))
    end do
  end function ztrace_dot

  !> P inner product of complex matrices
  function zp_dot(P, R0, w, s)
    real(dp),    allocatable, intent(in) :: P(:,:,:)
    complex(dp), allocatable, intent(in) :: R0(:,:)
    complex(dp), intent(in)              :: w(:,:)
    integer, intent(in)             :: s
    complex(dp)                          :: zp_dot(s)

    complex(dp)                          :: v(s)
    integer                         :: j, k, N, low(s), up(s), step, nrhs

    if (allocated(P)) then
      ! Biortho: P has orthogonal random numbers
      do i = 1, s
        v(i) = zero
        do k = 1, size(w, 2)
          v(i) = v(i) + cmplx(ddotprod(P(:, k, i), dble(w(:, k))), ddotprod(P(:, k, i), aimag(w(:, k))), kind=dp)
        end do
      end do
    else if (allocated(R0)) then
      ! BiCGSTAB: shadow vector equal to initial residual
      v(1) = zero
      do k = 1, size(w, 2)
        v(1) = v(1) + zdotprod(R0(:, k), w(:, k))
      end do
    else
      ! Minsync: P is piecewise constant
      ! WARNING: the integrals are done here in a peculiar way, not consistent with the
      ! definition of the dot product. So this is probably not working.
      N    = size(w,1)
      nrhs = size(w,2)
      step = N / s
      low(1) = 1
      do i = 1, s-1
        low(i+1) = i * step + 1
        up(i) = i * step
      end do
      up(s) = N

      do i = 1, s
        v(i)  = zero
        do j = 1, nrhs
          v(i) = v(i) + sum(w(low(i):up(i),j))
        end do
      end do
    end if

    zp_dot = v
  end function zp_dot

  !> Frobenius norm of complex matrix
  function dfrob_norm(v)
    real(dp), intent(in) :: v(:,:)
    real(dp)             :: dfrob_norm
    integer :: k
    dfrob_norm = 0.0_dp
    do k = 1, size(v, 2)
      dfrob_norm = dfrob_norm + ddotprod(v(:, k), v(:, k))
    end do
    dfrob_norm = sqrt(dfrob_norm)
  end function dfrob_norm

  !> Frobenius norm of complex matrix
  function zfrob_norm(v)
    complex(dp), intent(in)      :: v(:,:)
    real(dp)                  :: zfrob_norm
    integer :: k
    zfrob_norm = zero
    do k = 1, size(v, 2)
      zfrob_norm = zfrob_norm + real(zdotprod(v(:, k), v(:, k)),dp)
    end do
    zfrob_norm = sqrt(zfrob_norm)
  end function zfrob_norm

end function zidrs




!======================================================================================
end module scitools_idrs
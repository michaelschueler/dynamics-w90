module scitools_cg_solver
!! Contains conjugate-gradient linear solvers
!======================================================================================
   use scitools_debug
   use scitools_def,only: dp,iu,one,zero
   use scitools_linalg,only: util_axpy, util_copy, util_axpy, util_scal
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: dconjugate_gradients, zconjugate_gradients

   interface dconjugate_gradients
      module procedure dsym_conjugate_gradients, dbi_conjugate_gradients
   end interface dconjugate_gradients

   interface zconjugate_gradients
      module procedure zsym_conjugate_gradients, zbi_conjugate_gradients
   end interface zconjugate_gradients
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine dsym_conjugate_gradients(np, x, b, op, dotp, iter, residue, threshold)
      integer, intent(in)      :: np
      real(dp),  intent(inout) :: x(:)
      real(dp),  intent(in)    :: b(:)
      interface
      subroutine op(x, y)
         use scitools_def,only: dp
         implicit none
         real(dp), intent(in)    :: x(:)
         real(dp), intent(inout) :: y(:)
      end subroutine op
      real(dp) function dotp(x, y)
         use scitools_def,only: dp
         implicit none
         real(dp), intent(in)    :: x(:)
         real(dp), intent(in)    :: y(:)
      end function dotp
      end interface
      integer,          intent(inout) :: iter
      real(dp),  optional, intent(in)    :: threshold
      real(dp), optional, intent(out)   :: residue
      real(dp), allocatable :: r(:), ax(:), p(:), ap(:)
      real(dp)              :: alpha, beta, gamma
      real(dp)                 :: threshold_
      integer               :: max_iter, ip

      threshold_ = 1.0e-6_dp
      if(present(threshold)) threshold_ = threshold

      ALLOCATE( r(1:np))
      ALLOCATE(ax(1:np))
      ALLOCATE( p(1:ubound(x, 1)))
      ALLOCATE(ap(1:np))

      ! Initial residue.
      call op(x, ax)
      do ip = 1, np
         r(ip) = b(ip) - ax(ip)
      end do

      ! Initial search direction.
      call util_copy(np, r, p)

      max_iter = iter
      iter = 1
      do while(iter < max_iter)
         gamma = dotp(r, r)
         if (abs(gamma) < threshold_**2) exit
         call op(p, ap)
         alpha   = gamma/dotp(p, ap)
         call util_axpy(np, -alpha, ap, r)
         call util_axpy(np, alpha, p, x)
         beta    = dotp(r, r)/gamma
         do ip = 1, np
            p(ip) = r(ip) + beta*p(ip)
         end do
         iter    = iter + 1
      end do
      if (present(residue)) residue = sqrt(abs(gamma))

      DEALLOCATE(r)
      DEALLOCATE(ax)
      DEALLOCATE(p)
      DEALLOCATE(ap)

   end subroutine dsym_conjugate_gradients
!--------------------------------------------------------------------------------------
   subroutine zsym_conjugate_gradients(np, x, b, op, dotp, iter, residue, threshold)
      integer, intent(in)      :: np
      complex(dp),  intent(inout) :: x(:)
      complex(dp),  intent(in)    :: b(:)
      interface
      subroutine op(x, y)
         use scitools_def,only: dp
         implicit none
         complex(dp), intent(in)    :: x(:)
         complex(dp), intent(inout) :: y(:)
      end subroutine op
      complex(dp) function dotp(x, y)
         use scitools_def,only: dp
         implicit none
         complex(dp), intent(in)    :: x(:)
         complex(dp), intent(in)    :: y(:)
      end function dotp
      end interface
      integer,          intent(inout) :: iter
      real(dp),  optional, intent(in)    :: threshold
      real(dp), optional, intent(out)    :: residue
      complex(dp), allocatable :: r(:), ax(:), p(:), ap(:)
      complex(dp)              :: alpha, beta, gamma
      real(dp)                 :: threshold_
      integer               :: max_iter, ip

      threshold_ = 1.0e-6_dp
      if(present(threshold)) threshold_ = threshold

      ALLOCATE( r(1:np))
      ALLOCATE(ax(1:np))
      ALLOCATE( p(1:ubound(x, 1)))
      ALLOCATE(ap(1:np))

      ! Initial residue.
      call op(x, ax)
      do ip = 1, np
         r(ip) = b(ip) - ax(ip)
      end do

      ! Initial search direction.
      call util_copy(np, r, p)

      max_iter = iter
      iter = 1
      do while(iter < max_iter)
         gamma = dotp(r, r)
         if (abs(gamma) < threshold_**2) exit
         call op(p, ap)
         alpha   = gamma/dotp(p, ap)
         call util_axpy(np, -alpha, ap, r)
         call util_axpy(np, alpha, p, x)
         beta    = dotp(r, r)/gamma
         do ip = 1, np
            p(ip) = r(ip) + beta*p(ip)
         end do
         iter    = iter + 1
      end do
      if (present(residue)) residue = sqrt(abs(gamma))

      DEALLOCATE(r)
      DEALLOCATE(ax)
      DEALLOCATE(p)
      DEALLOCATE(ap)

   end subroutine zsym_conjugate_gradients
!--------------------------------------------------------------------------------------
   subroutine dbi_conjugate_gradients(np, x, b, op, opt, dotp, iter, residue, threshold)
      integer, intent(in)    :: np
      real(dp),  intent(inout) :: x(:)
      real(dp),  intent(in)    :: b(:)
      interface
         subroutine op(x, y)
            use scitools_def,only: dp
            implicit none
            real(dp), intent(in)    :: x(:)
            real(dp), intent(inout) :: y(:)
         end subroutine op
      end interface
      interface
         subroutine opt(x, y)
            use scitools_def,only: dp
            implicit none
            real(dp), intent(in)    :: x(:)
            real(dp), intent(inout) :: y(:)
         end subroutine opt
         real(dp) function dotp(x, y)
            use scitools_def,only: dp
            implicit none
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: y(:)
         end function dotp
      end interface
      integer,          intent(inout) :: iter
      real(dp), optional, intent(out)   :: residue
      real(dp), optional,  intent(in)    :: threshold

      real(dp), allocatable    :: r(:), rr(:), ax(:), p(:), pp(:), ap(:), atp(:)
      real(dp)                 :: alpha, beta, gamma, err, threshold_
      integer                  :: max_iter


      threshold_ = 1.0e-6_dp
      if(present(threshold)) threshold_ = threshold

      ALLOCATE(  r(1:np))
      ALLOCATE( rr(1:np))
      ALLOCATE( ax(1:np))
      ALLOCATE(  p(1:ubound(x, 1)))
      ALLOCATE( pp(1:ubound(x, 1)))
      ALLOCATE( ap(1:np))
      ALLOCATE(atp(1:np))

      ! Initial residue.
      call op(x, ax)
      ! r <- b - ax, rr <- r
      call util_copy(np, b, r)
      call util_axpy(np, -1.0_dp, ax, r)
      call util_copy(np, r, rr)

      ! Initial search direction.
      call util_copy(np, r, p)
      call util_copy(np, p, pp)

      max_iter = iter
      iter     = 1
      do while(iter < max_iter)
         gamma = dble(dotp(rr, r))
         err   = dble(dotp(r, r))
         if (abs(err) < threshold_**2) exit
         call op (p,  ap)
         call opt(pp, atp)
         alpha = gamma/dble(dotp(pp, ap))
         call util_axpy(np, -alpha, ap, r)         ! r  <- r - alpha*ap
         call util_axpy(np, -alpha, atp, rr)       ! rr <- rr - alpha*atp
         call util_axpy(np, alpha, p, x)           ! x  <- x + alpha*p
         beta = dble(dotp(rr, r))/gamma
         call util_scal(np, beta, p)               ! p  <- r + beta*p
         call util_axpy(np, 1.0_dp, r, p)
         call util_scal(np, beta, pp)              ! pp <- rr + beta*pp
         call util_axpy(np, 1.0_dp, rr, pp)
         iter = iter + 1
      end do
      if (present(residue)) residue = sqrt(abs(err))

      DEALLOCATE(r)
      DEALLOCATE(rr)
      DEALLOCATE(ax)
      DEALLOCATE(p)
      DEALLOCATE(pp)
      DEALLOCATE(ap)
      DEALLOCATE(atp)

   end subroutine dbi_conjugate_gradients
!--------------------------------------------------------------------------------------
   subroutine zbi_conjugate_gradients(np, x, b, op, opt, dotp, iter, residue, threshold)
      integer, intent(in)    :: np
      complex(dp),  intent(inout) :: x(:)
      complex(dp),  intent(in)    :: b(:)
      interface
         subroutine op(x, y)
            use scitools_def,only: dp
            implicit none
            complex(dp), intent(in)    :: x(:)
            complex(dp), intent(inout)   :: y(:)
         end subroutine op
      end interface
      interface
         subroutine opt(x, y)
            use scitools_def,only: dp
            implicit none
            complex(dp), intent(in)    :: x(:)
            complex(dp), intent(inout)   :: y(:)
         end subroutine opt
         complex(dp) function dotp(x, y)
            use scitools_def,only: dp
            implicit none
            complex(dp), intent(in) :: x(:)
            complex(dp), intent(in) :: y(:)
         end function dotp
      end interface
      integer,          intent(inout) :: iter
      real(dp), optional, intent(out)    :: residue
      real(dp), optional,  intent(in)    :: threshold

      complex(dp), allocatable :: r(:), rr(:), ax(:), p(:), pp(:), ap(:), atp(:)
      real(dp)                 :: alpha, beta, gamma, err, threshold_
      integer                  :: max_iter


      threshold_ = 1.0e-6_dp
      if(present(threshold)) threshold_ = threshold

      ALLOCATE(  r(1:np))
      ALLOCATE( rr(1:np))
      ALLOCATE( ax(1:np))
      ALLOCATE(  p(1:ubound(x, 1)))
      ALLOCATE( pp(1:ubound(x, 1)))
      ALLOCATE( ap(1:np))
      ALLOCATE(atp(1:np))

      ! Initial residue.
      call op(x, ax)
      ! r <- b - ax, rr <- r
      call util_copy(np, b, r)
      call util_axpy(np, -one, ax, r)
      call util_copy(np, r, rr)

      ! Initial search direction.
      call util_copy(np, r, p)
      call util_copy(np, p, pp)

      max_iter = iter
      iter     = 1
      err      = 1.0e8_dp
      do while(iter < max_iter .and. abs(err) > threshold_**2)
         gamma = dble(dotp(rr, r))
         err   = dble(dotp(r, r))
         ! if (abs(err) < threshold_**2) exit
         call op (p,  ap)
         call opt(pp, atp)
         alpha = gamma/dble(dotp(pp, ap))
         call util_axpy(np, -alpha*one, ap, r)         ! r  <- r - alpha*ap
         call util_axpy(np, -alpha*one, atp, rr)       ! rr <- rr - alpha*atp
         call util_axpy(np, alpha*one, p, x)           ! x  <- x + alpha*p
         beta = dble(dotp(rr, r))/gamma
         call util_scal(np, beta*one, p)               ! p  <- r + beta*p
         call util_axpy(np, one, r, p)
         call util_scal(np, beta*one, pp)              ! pp <- rr + beta*pp
         call util_axpy(np, one, rr, pp)
         iter = iter + 1
      end do
      if (present(residue)) residue = sqrt(abs(err))

      DEALLOCATE(r)
      DEALLOCATE(rr)
      DEALLOCATE(ax)
      DEALLOCATE(p)
      DEALLOCATE(pp)
      DEALLOCATE(ap)
      DEALLOCATE(atp)

   end subroutine zbi_conjugate_gradients
!--------------------------------------------------------------------------------------

!======================================================================================
end module scitools_cg_solver

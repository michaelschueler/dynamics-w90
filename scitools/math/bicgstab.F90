module scitools_bicgstab
!! Contains stabilized biconjugate-gradient linear solvers
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use scitools_debug
   use scitools_def,only: dp,iu,one,zero
   use scitools_linalg,only: util_axpy, util_copy
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: bicgstab
!--------------------------------------------------------------------------------------
   interface bicgstab
      module procedure dbicgstab, zbicgstab
   end interface bicgstab
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
  ! need to solve asymmetric matrix
  ! http://www.jicfus.jp/wiki/index.php?Bi-CGSTAB%20%E6%B3%95
  ! https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
!--------------------------------------------------------------------------------------
   subroutine dbicgstab(ndim, x, b, matvec, psolve, dotp, info, maxiter, tol, res)
      implicit none
      integer,intent(in) :: ndim
      real(dp),dimension(ndim),intent(inout) :: x
      real(dp),dimension(ndim),intent(in)    :: b
      interface
         subroutine matvec(x,y)
            import :: dp
            real(dp),intent(in) :: x(:)
            real(dp),intent(inout) :: y(:)
         end subroutine matvec

         subroutine psolve(x)
            import :: dp
            real(dp),intent(inout) :: x(:)
         end subroutine psolve

         real(dp) function dotp(x,y)
            import :: dp
            real(dp),intent(in) :: x(:),y(:)         
         end function dotp
      end interface
      integer,intent(inout) :: info
      integer,intent(in),optional :: maxiter
      real(dp),intent(in),optional :: tol
      real(dp),intent(out),optional :: res
      integer :: iter_max
      real(dp) :: tol_,res_
      integer :: i, j, iter
      real(dp),allocatable,dimension(:) :: r, r_new, v, v_new, p, p_new
      real(dp),allocatable,dimension(:) :: r0, h, s, t
      real(dp),allocatable,dimension(:) :: rr
      real(dp),allocatable,dimension(:) :: xx, xx_new, yy
      real(dp) :: rho, rho_new, omega, omega_new
      real(dp) :: alpha, beta
      real(dp) :: r0r, t2, r0v, ts
      real(dp),allocatable,dimension(:) :: pinv_b, Ax

      iter_max = 1000
      if(present(maxiter)) iter_max = maxiter

      tol_ = 1.0e-6_dp
      if(present(tol)) tol_ = tol

      allocate(r(ndim))
      allocate(r_new(ndim))
      allocate(v(ndim))
      allocate(v_new(ndim))
      allocate(p(ndim))
      allocate(p_new(ndim))
      allocate(r0(ndim))
      allocate(h(ndim))
      allocate(s(ndim))
      allocate(t(ndim))
      allocate(xx(ndim))
      allocate(xx_new(ndim))
      allocate(yy(ndim))
      allocate(pinv_b(ndim))
      allocate(Ax(ndim))


      ! left precondition
      ! (P^{-1}A)x = P^{-1}b
      ! P^{-1}b = b_i/A_{ii}
      call util_copy(ndim, b, pinv_b)
      call psolve(pinv_b)

      ! P^{-1}A

      ! initial guess
      call util_copy(ndim, x, xx)

      ! residual
      call matvec(xx, yy)
      call psolve(yy)
      call util_copy(ndim, pinv_b, r)
      call util_axpy(ndim, -1.0_dp, yy, r)
      ! r = pinv_b - yy

      call util_copy(ndim, r, r0)
      r0r = dotp(r0, r)

      if (r0r == 0.0d0) then
         info = 2
         return
      end if
    
      rho   = 1.0d0
      alpha = 1.0d0
      omega = 1.0d0

      v = 0.0d0
      p = 0.0d0

      do iter = 1, iter_max
         rho_new = dotp(r0, r)
         beta = (rho_new/rho)*(alpha/omega)
         ! p_new(i) = r(i) + beta*(p(i)-omega*v(i))
         call util_copy(ndim, r, p_new)
         call util_axpy(ndim, beta, p, p_new)
         call util_axpy(ndim, -beta*omega, v, p_new)

         call matvec(p_new, v_new)
         call psolve(v_new)
         r0v = dotp(r0, v_new)

         alpha = rho_new/r0v

         call util_copy(ndim, xx, h)
         call util_axpy(ndim, alpha, p_new, h)
         !  h(i) = xx(i) + alpha*p_new(i)

         call util_copy(ndim, pinv_b, rr)
         call matvec(h, Ax)
         call psolve(Ax)        
         call util_axpy(ndim, -1.0_dp, Ax, rr)

         res_ = dotp(rr, rr)
         res_ = sqrt(res_)
         if (res_ <= tol_) then
            call util_copy(ndim, h, xx_new)
            exit
         end if

         call util_copy(ndim, r, s)
         call util_axpy(ndim, -alpha, v_new, s)
         ! s(i) = r(i) - alpha*v_new(i)

         call matvec(s, t)
         call psolve(t)

         ts = dotp(t, s)
         t2 = dotp(t, t)
         omega_new = ts/t2

         call util_copy(ndim, h, xx_new)
         call util_axpy(ndim, omega_new, s, xx_new)
         ! xx_new(i) = h(i) + omega_new*s(i)

         call util_copy(ndim, pinv_b, rr)
         call matvec(xx_new, Ax)
         call psolve(Ax)   
         call util_axpy(ndim, -1.0_dp, Ax, rr)         

         res_ = sqrt(dotp(rr, rr))

         if (res_ <= tol) exit

         call util_copy(ndim, s, r_new)
         call util_axpy(ndim, -omega_new, t, r_new)

         call util_copy(ndim, xx_new, xx)
         call util_copy(ndim, p_new, p)
         call util_copy(ndim, r_new, r)
         call util_copy(ndim, v_new, v)

         rho   = rho_new
         omega = omega_new
      end do ! iter

      info = 0
      if (iter >= iter_max .and. res_>tol_) then
         info = 1
      end if

      call util_copy(ndim, xx_new, x)

      deallocate(r)
      deallocate(r_new)
      deallocate(v)
      deallocate(v_new)
      deallocate(p)
      deallocate(p_new)
      deallocate(r0)
      deallocate(h)
      deallocate(s)
      deallocate(t)
      deallocate(xx)
      deallocate(xx_new)
      deallocate(yy)
      deallocate(pinv_b)
      deallocate(Ax)

   end subroutine dbicgstab
!--------------------------------------------------------------------------------------
   subroutine zbicgstab(ndim, x, b, matvec, psolve, dotp, iter, info, tol, res)
      implicit none
      integer,intent(in) :: ndim
      complex(dp),dimension(:),intent(inout) :: x
      complex(dp),dimension(:),intent(in)    :: b
      interface
         subroutine matvec(x,y)
            import :: dp
            complex(dp),intent(in) :: x(:)
            complex(dp),intent(inout) :: y(:)
         end subroutine matvec

         subroutine psolve(x)
            import :: dp
            complex(dp),intent(inout) :: x(:)
         end subroutine psolve

         complex(dp) function dotp(x,y)
            import :: dp
            complex(dp),intent(in) :: x(:),y(:)         
         end function dotp
      end interface
      integer,intent(inout) :: iter
      integer,intent(inout) :: info
      real(dp),intent(in),optional :: tol
      real(dp),intent(out),optional :: res
      integer :: iter_max
      real(dp) :: tol_,res_
      integer :: i, j, its
      complex(dp),allocatable,dimension(:) :: r, r_new, v, v_new, p, p_new
      complex(dp),allocatable,dimension(:) :: r0, h, s, t
      complex(dp),allocatable,dimension(:) :: rr
      complex(dp),allocatable,dimension(:) :: xx, xx_new
      real(dp) :: rho, rho_new, omega, omega_new
      real(dp) :: alpha, beta
      real(dp) :: r0r, t2, r0v, ts
      complex(dp),allocatable,dimension(:) :: pinv_b, Ax

      iter_max = iter

      tol_ = 1.0e-6_dp
      if(present(tol)) tol_ = tol

      allocate(r(ndim))
      allocate(r_new(ndim))
      allocate(v(ndim))
      allocate(v_new(ndim))
      allocate(p(ndim))
      allocate(p_new(ndim))
      allocate(r0(ndim))
      allocate(h(ndim))
      allocate(s(ndim))
      allocate(t(ndim))
      allocate(xx(ndim))
      allocate(xx_new(ndim))
      allocate(pinv_b(ndim))
      allocate(Ax(ndim))

      ! left precondition
      ! (P^{-1}A)x = P^{-1}b
      ! P^{-1}b = b_i/A_{ii}
      ! call util_copy(ndim, b, pinv_b)
      ! call ZCOPY(ndim, b, 1, pinv_b, 1)
      pinv_b = b
      call psolve(pinv_b)

      ! P^{-1}A

      ! initial guess
      ! call util_copy(ndim, x, xx)
      xx = x

      ! residual
      call matvec(xx, Ax)
      call psolve(Ax)
      ! call util_copy(ndim, pinv_b, r)
      ! call util_axpy(ndim, -one, Ax, r)
      r = pinv_b - Ax

      ! call util_copy(ndim, r, r0)
      r0 = r
      r0r = dble(dotp(r0, r))

      if (r0r == 0.0d0) then
         info = 2
         return
      end if
    
      rho   = 1.0d0
      alpha = 1.0d0
      omega = 1.0d0

      v = zero
      p = zero

      do its = 1, iter_max
         rho_new = dble(dotp(r0, r))
         beta = (rho_new/rho)*(alpha/omega)
         p_new(:) = r(:) + beta*(p(:)-omega*v(:))
         ! call util_copy(ndim, r, p_new)
         ! call util_axpy(ndim, beta*one, p, p_new)
         ! call util_axpy(ndim, -beta*omega*one, v, p_new)

         call matvec(p_new, v_new)
         call psolve(v_new)
         r0v = dble(dotp(r0, v_new))

         alpha = rho_new/r0v

         ! call util_copy(ndim, xx, h)
         ! call util_axpy(ndim, alpha*one, p_new, h)
         h(:) = xx(:) + alpha*p_new(:)

         ! call util_copy(ndim, pinv_b, rr)
         call matvec(h, Ax)
         call psolve(Ax)        
         ! call util_axpy(ndim, -one, Ax, rr)
         rr = pinv_b - Ax

         res_ = abs(dotp(rr, rr))
         res_ = sqrt(res_)

         if(info == 1) then
            write(output_unit,'("bicgstab: it = ",i6," res = ",E12.5)') its, res_
         end if

         if (res_ <= tol_) then
            ! call util_copy(ndim, h, xx_new)
            xx_new = h
            exit
         end if

         ! call util_copy(ndim, r, s)
         ! call util_axpy(ndim, -alpha*one, v_new, s)
         s(:) = r(:) - alpha*v_new(:)

         call matvec(s, t)
         call psolve(t)

         ts = dble(dotp(t, s))
         t2 = dble(dotp(t, t))
         omega_new = ts/t2

         ! call util_copy(ndim, h, xx_new)
         ! call util_axpy(ndim, omega_new*one, s, xx_new)
         xx_new(:) = h(:) + omega_new*s(:)

         ! call util_copy(ndim, pinv_b, rr)
         call matvec(xx_new, Ax)
         call psolve(Ax)   
         ! call util_axpy(ndim, -one, Ax, rr)
         rr = pinv_b - Ax         

         res_ = sqrt(abs(dotp(rr, rr)))

         if (res_ <= tol_) exit

         ! call util_copy(ndim, s, r_new)
         ! call util_axpy(ndim, -omega_new*one, t, r_new)
         r_new(:) = s(:) - omega_new*t(:)

         xx    = xx_new
         p     = p_new
         r     = r_new
         v     = v_new
 
         ! call util_copy(ndim, xx_new, xx)
         ! call util_copy(ndim, p_new, p)
         ! call util_copy(ndim, r_new, r)
         ! call util_copy(ndim, v_new, v)

         rho   = rho_new
         omega = omega_new
      end do ! iter

      iter = its

      info = 0
      if (its >= iter_max .and. res_>tol_) then
         info = 1
      end if

      ! call util_copy(ndim, xx_new, x)
      x(:) = xx_new(:)
      if(present(res)) res = res_

      deallocate(r)
      deallocate(r_new)
      deallocate(v)
      deallocate(v_new)
      deallocate(p)
      deallocate(p_new)
      deallocate(r0)
      deallocate(h)
      deallocate(s)
      deallocate(t)
      deallocate(xx)
      deallocate(xx_new)
      deallocate(pinv_b)
      deallocate(Ax)

   end subroutine zbicgstab
!--------------------------------------------------------------------------------------

!======================================================================================
end module scitools_bicgstab
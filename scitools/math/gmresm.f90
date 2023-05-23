module scitools_gmresm
!! Contains real and complex versions of the GMRESM algorithm for iteratively 
!! solving linear equations
!======================================================================================
   use scitools_debug
   use scitools_def,only: dp,iu,one,zero
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: gmresm
!--------------------------------------------------------------------------------------
   interface gmresm
      module procedure dgmresm, zgmresm
   end interface gmresm
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine dgmresm(m,n,x,b,matvec,psolve,dotprd,res,del,its,info)
      implicit none
      integer,          intent(in)    :: m
      integer,          intent(in)    :: n
      real(dp), intent(inout)         :: x(:)
      real(dp), intent(in)            :: b(:)
      interface
      subroutine matvec(x,y)
         import :: dp
         real(dp),intent(in)    :: x(:)
         real(dp),intent(inout) :: y(:)
      end subroutine matvec

      subroutine psolve(x)
         import :: dp
         real(dp),intent(inout) :: x(:)
      end subroutine psolve

      function dotprd(x,y) result(d)
         import :: dp
         real(dp),intent(in)    :: x(:),y(:)
         real(dp) :: d
      end function dotprd
      end interface 
      real(dp), intent(inout) :: res
      real(dp), intent(inout) :: del
      integer,          intent(inout) :: its
      integer,          intent(inout) :: info
      real(dp) :: tol,res_,stgn
      real(dp) :: h_(m+1,m), y(m+1), p(m+1), work(4*m+1)
      real(dp),allocatable :: w(:), z(:)
      real(dp),allocatable :: h(:,:),v(:,:)
      integer :: imx, piv(m), rank, i
      real(dp), save :: beta
      integer, save :: j
      logical :: done   

      allocate(h(m+1,m)); h = 0.0d0
      allocate(v(n,m+1)); v = 0.0d0
      allocate(w(n)); w = 0.0d0
      allocate(z(n)); z = 0.0d0

      if(info==2) then
         call hookstep(j,h,m,beta,del, y)
         ! z = matmul(v(:,1:j),y(1:j))
         call DGEMV('N',n,j,1.0_dp,v,n,y,1,0.0_dp,z,1)
         call psolve(z)
         x = z
         info = 0
         return
      end if  

      tol = res
      imx = its
      its = 0
      v   = 0d0

      1 continue
      res_ = 1d99
      stgn = 1d0 - 1d-14

      beta = sqrt(dotprd(x,x)) 
      if(beta==0d0)  w = 0d0
      if(beta/=0d0)  call matvec(x, w)
      w = b - w
      beta = sqrt(dotprd(w,w)) 
      v(:,1) = w / beta

      h = 0d0
      do j = 1, m
         its = its + 1
         z = v(:,j)      
         call psolve(z)
         call matvec(z, w)
         do i = 1, j
            h(i,j) = dotprd(w,v(:,i))
            w = w - h(i,j)*v(:,i)
         end do
         h(j+1,j) = sqrt(dotprd(w,w))
         v(:,j+1) = w / h(j+1,j)

         p(1) = beta
         p(2:j+1) = 0d0
         h_(1:j+1,1:j) = h(1:j+1,1:j)
         call dgelsy(j+1,j,1,h_,m+1,p,m+1,piv,m,rank,work,4*m+1,i)
         if(i/=0) stop 'gmresm: dgelsy'
         y = p

         ! p(1:j+1) = - matmul(h(1:j+1,1:j),y(1:j))
         call DGEMV('N',j+1,j,-1.0_dp,h,m+1,y,1,0.0_dp,p,1)
         p(1) = p(1) + beta
         res = sqrt(dot_product(p(1:j+1),p(1:j+1)))
         if(info==1) print*, 'gmresm: it=', its,' res=', real(res,dp)

         done = (res<=tol .or. its==imx .or. res>res_)
         if(done .or. j==m) then
            if(del>0d0)  call hookstep(j,h,m,beta,del, y)
            ! z = matmul(v(:,1:j),y(1:j))
            call DGEMV('N',n,j,1.0_dp,v,n,y,1,0.0_dp,z,1)
            call psolve(z)
            x = x + z
            if(its==imx) info = 2
            if(res>res_) info = 1
            if(res<=tol) info = 0
            if(done)     return
            if(del>0d0)  print*, 'gmres: warning! restart affects hookstep'
            goto 1       ! (j==m) restart
         end if
         res_ = res*stgn

      end do   

      deallocate(h,v)

   end subroutine dgmresm
!--------------------------------------------------------------------------------------
   subroutine zgmresm(m,n,x,b,matvec,psolve,dotprd,res,its,info)
      implicit none
      complex(dp),parameter :: mone=-one
      integer,     intent(in)    :: m
      integer,     intent(in)    :: n
      complex(dp), intent(inout) :: x(:)
      complex(dp), intent(in)    :: b(:)
      interface
      subroutine matvec(x,y)
         import :: dp
         complex(dp),intent(in)    :: x(:)
         complex(dp),intent(inout) :: y(:)
      end subroutine matvec

      subroutine psolve(x)
         import :: dp
         complex(dp),intent(inout) :: x(:)
      end subroutine psolve

      function dotprd(x,y) result(d)
         import :: dp
         complex(dp),intent(in)    :: x(:),y(:)
         complex(dp) :: d
      end function dotprd
      end interface  
      real(dp),         intent(inout) :: res
      integer,          intent(inout) :: its
      integer,          intent(inout) :: info
      real(dp) :: tol,res_,stgn
      complex(dp) :: h_(m+1,m), y(m+1), p(m+1)
      real(dp),allocatable :: rwork(:)
      complex(dp),allocatable :: w(:), z(:)
      complex(dp),allocatable :: h(:,:),v(:,:),work(:)
      integer :: imx, piv(m), rank, i, lwork
      real(dp), save :: beta
      integer, save :: j
      logical :: done   

      allocate(h(m+1,m)); h = (0.0d0,0.0d0)
      allocate(v(n,m+1)); v = (0.0d0,0.0d0)
      allocate(w(n)); w = (0.0d0,0.0d0)
      allocate(z(n)); z = (0.0d0,0.0d0)

      lwork = -1
      allocate(work(4*m+1),rwork(2*(m+1)))
      call zgelsy(m+1,m,1,h_,m+1,p,m+1,piv,m,rank,work,lwork,rwork,i)

      lwork = int(dble(work(1)))

      deallocate(work)
      allocate(work(lwork))

      tol = res
      imx = its
      its = 0

      1 continue
      res_ = 1d99
      stgn = 1d0 - 1d-14

      beta = sqrt((dotprd(x,x)))
      if(beta==0d0)  w = (0d0,0d0)
      if(beta/=0d0)  call matvec(x, w)
      w = b - w
      beta = sqrt(abs(dotprd(w,w)) )
      v(:,1) = w / beta

      do j = 1, m
         its = its + 1
         z = v(:,j)      
         call psolve(z)
         call matvec(z, w)
         do i = 1, j
            h(i,j) = dotprd(w,v(:,i))
            w = w - h(i,j)*v(:,i)
         end do
         h(j+1,j) = sqrt(abs(dotprd(w,w)))
         v(:,j+1) = w / h(j+1,j)

         p(1) = beta
         p(2:j+1) = (0d0,0d0)
         h_(1:j+1,1:j) = h(1:j+1,1:j)
         call zgelsy(j+1,j,1,h_,m+1,p,m+1,piv,m,rank,work,lwork,rwork,i)
         if(i/=0) stop 'gmresm: zgelsy'
         y = p

         ! p(1:j+1) = - matmul(h(1:j+1,1:j),y(1:j))
         call ZGEMV('N',j+1,j,mone,h,m+1,y,1,zero,p,1)
         p(1) = p(1) + beta
         res = sqrt(abs(dot_product(p(1:j+1),p(1:j+1))))
         if(info==1) print*, 'gmresm: it=', its,' res=', real(res,dp)

         ! done = (res<=tol .or. its==imx .or. res>res_)
         done = (res<=tol .or. its==imx)
         if(done .or. j==m) then
            ! z = matmul(v(:,1:j),y(1:j))
            call ZGEMV('N',n,j,one,v,n,y,1,zero,z,1)
            call psolve(z)
            x = x + z
            if(its==imx) info = 2
            if(res>res_) info = 1
            if(res<=tol) info = 0
            if(done)     return
            goto 1       ! (j==m) restart
         end if

         res_ = res*stgn

      end do   

      deallocate(h,v,w,z)

      deallocate(work,rwork)

   end subroutine zgmresm
!--------------------------------------------------------------------------------------

!-----------------------------------------------------------------
! replace y with a vector that generates a hookstep
! c.f. Viswanath (2008) arXiv:0809.1498
!-----------------------------------------------------------------
 subroutine hookstep(j,h,m,beta,del, y)
   implicit none
   integer,          intent(in)    :: j, m
   double precision, intent(in)    :: h(:,:)
   double precision, intent(in)    :: beta
   double precision, intent(inout) :: del
   double precision, intent(out)   :: y(j)
   double precision :: a(j+1,j), s(j), u(j+1,j+1), vt(j,j), work(5*(j+1))
   double precision :: p(j+1), q(j), mu, qn
   integer :: info
   
   a = h(1:j+1,1:j)
   
   call dgesvd('A','A',j+1,j,a,j+1,s,u,j+1,vt,j,work,5*(j+1),info)
   if(info/=0) stop 'hookstep: dgesvd'
   
   p(1:j) = beta * u(1,1:j)   

   mu = max(s(j)*s(j)*1d-6,1d-99)
   qn = 1d99
   do while(qn>del)
      mu = mu * 1.1d0
      q = p(1:j)*s/(mu+s*s)
      qn = sqrt(dot_product(q,q))
   end do

   y = matmul(q,vt)

   p = - matmul(h(1:j+1,1:j),y(1:j))
   p(1) = p(1) + beta
   del = sqrt(dot_product(p,p))
 
 end subroutine hookstep
!--------------------------------------------------------------------------------------
end module scitools_gmresm

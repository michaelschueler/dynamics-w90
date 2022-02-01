module Mintegration
!======================================================================================
   use Mdebug
   use Mdef,only: dp,iu,one,zero
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: GetGregW_1D,GetGregW_2D_Trg,GregoryIntegral,GregoryDot,&
      GregoryIntegrate,GregoryCumIntegrate,GetGregoryWeights

   interface GregoryIntegral
      module procedure DGregoryIntegral,ZGregoryIntegral
   end interface GregoryIntegral

   interface GregoryDot
      module procedure GregoryDot_re_re,GregoryDot_c_re,GregoryDot_re_c,GregoryDot_c_c
   end interface GregoryDot
!--------------------------------------------------------------------------------------
   integer,parameter :: p=5
   real(dp),parameter :: wam(0:p-1)=[-(19D0/720D0), 53D0/360D0, -(11D0/30D0), 323D0/360D0, 251D0/720D0]
   real(dp),parameter :: omega(0:p-1) = [251D0/720D0, 299D0/240D0, 211D0/240D0, 739D0/720D0, 1D0]
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   pure function GetStartW() result(w)
      real(dp) :: w(0:p-1,0:p-1)

      w = 0.0_dp
      w(1,0:4) = [251D0/720D0,323D0/360D0,-11D0/30D0,53D0/360D0,-19D0/720D0]
      w(2,0:4) = [29D0/90D0,62D0/45D0,4D0/15D0,2D0/45D0,-1D0/90D0]
      w(3,0:4) = [27D0/80D0, 51D0/40D0, 9D0/10D0, 21D0/40D0, -(3D0/80D0)]
      w(4,0:4) = [14D0/45D0, 64D0/45D0, 8D0/15D0, 64D0/45D0, 14D0/45D0]

   end function GetStartW
!--------------------------------------------------------------------------------------
   subroutine GetGregoryWeights(n,w)
      integer,intent(in)     :: n
      real(dp),intent(inout) :: w(0:)
      integer :: j,k
      real(dp) :: wstart(0:2*p-2,0:p-1)

      w = 0.0_dp

      wstart(0:p-1,0:p-1) = GetStartW()
      do k=0,p-2
         do j=0,k
            wstart(p+k,j) = wstart(p+k-1,j)
         end do
         do j=k+1,p-1
            wstart(p+k,j) = wstart(p+k-1,j) + wam(j-k-1)
         end do
      end do

      if(n <= 2*p-2) then
         do j=0,p-1
            w(j) = wstart(n,j)
         end do
      else
         do j=0,p-1
            w(j) = wstart(2*p-2,j)
         end do
      end if
      if(n >= p) then
         do j=max(n-p+1,p),n
            w(j) = omega(n-j)
         end do
         do j=p,n-p
            w(j) = 1.0_dp
         end do
      end if

   end subroutine GetGregoryWeights
!--------------------------------------------------------------------------------------
   function GregoryIntegrate(n,h,f) result(fint)
      integer,parameter :: p=5
      integer,intent(in)   :: n
      real(dp),intent(in)  :: h
      real(dp),intent(in)  :: f(0:)
      real(dp) :: fint
      integer :: j,k
      real(dp) :: w
      real(dp) :: wstart(0:2*p-2,0:p-1)

      fint = 0.0_dp
      wstart(0:p-1,0:p-1) = GetStartW()
      do k=0,p-2
         do j=0,k
            wstart(p+k,j) = wstart(p+k-1,j)
         end do
         do j=k+1,p-1
            wstart(p+k,j) = wstart(p+k-1,j) + wam(j-k-1)
         end do
      end do

      if(n <= 2*p-2) then
         do j=0,p-1
            w = wstart(n,j)
            fint = fint + w * f(j)
         end do
      else
         do j=0,p-1
            w = wstart(2*p-2,j)
            fint = fint + w * f(j)
         end do
      end if
      if(n >= p) then
         do j=max(n-p+1,p),n
            w = omega(n-j)
            fint = fint + w * f(j)
         end do
         do j=p,n-p
            fint = fint + f(j)
         end do
      end if

      fint = h * fint

   end function GregoryIntegrate
!--------------------------------------------------------------------------------------
   function GregoryCumIntegrate(nmax,h,f) result(fint)
      integer,parameter :: p=5
      integer,intent(in)   :: nmax
      real(dp),intent(in)  :: h
      real(dp),intent(in)  :: f(0:)
      real(dp) :: fint(0:nmax)
      integer :: n,j,k
      real(dp) :: w
      real(dp) :: wstart(0:2*p-2,0:p-1)

      fint = 0.0_dp
      wstart(0:p-1,0:p-1) = GetStartW()
      do k=0,p-2
         do j=0,k
            wstart(p+k,j) = wstart(p+k-1,j)
         end do
         do j=k+1,p-1
            wstart(p+k,j) = wstart(p+k-1,j) + wam(j-k-1)
         end do
      end do

      do n=0,nmax
         if(n <= 2*p-2) then
            do j=0,p-1
               w = wstart(n,j)
               fint(n) = fint(n) + w * f(j)
            end do
         else
            do j=0,p-1
               w = wstart(2*p-2,j)
               fint(n) = fint(n) + w * f(j)
            end do
         end if
         if(n >= p) then
            do j=max(n-p+1,p),n
               w = omega(n-j)
               fint(n) = fint(n) + w * f(j)
            end do
            do j=p,n-p
               fint(n) = fint(n) + f(j)
            end do
         end if
      end do

      fint = h * fint

   end function GregoryCumIntegrate
!--------------------------------------------------------------------------------------
   subroutine GetGregW_1D(nmax,w,transp)
      integer,parameter :: p=5
      integer,intent(in)     :: nmax
      real(dp),allocatable,intent(inout) :: w(:,:) 
      logical,intent(in),optional :: transp
      logical :: transp_ = .false.
      integer :: n,k,j
      real(dp) :: wstart(0:2*p-2,0:p-1)

      if(present(transp)) transp_ = transp

      if(.not.allocated(w)) allocate(w(0:nmax,0:nmax))

      w = 0.0_dp

      wstart(0:p-1,0:p-1) = GetStartW()
      do k=0,p-2
         do j=0,k
            wstart(p+k,j) = wstart(p+k-1,j)
         end do
         do j=k+1,p-1
            wstart(p+k,j) = wstart(p+k-1,j) + wam(j-k-1)
         end do
      end do

      if(transp_) then
         do n=0,nmax
            if(n <= 2*p-2) then
               do j=0,p-1
                  w(j,n) = wstart(n,j)
               end do
            else
               do j=0,p-1
                  w(j,n) = wstart(2*p-2,j)
               end do
            end if
            if(n >= p) then
               do j=max(n-p+1,p),n
                  w(j,n) = omega(n-j)
               end do

               do j=p,n-p
                  w(j,n) = 1.0_dp
               end do
            end if
         end do
      else
         do n=0,nmax
            if(n <= 2*p-2) then
               do j=0,p-1
                  w(n,j) = wstart(n,j)
               end do
            else
               do j=0,p-1
                  w(n,j) = wstart(2*p-2,j)
               end do
            end if
            if(n >= p) then
               do j=max(n-p+1,p),n
                  w(n,j) = omega(n-j)
               end do

               do j=p,n-p
                  w(n,j) = 1.0_dp
               end do
            end if
         end do
      end if

   end subroutine GetGregW_1D 
!--------------------------------------------------------------------------------------
   subroutine GetGregW_2D_Trg(nmax,w2d)
      integer,intent(in)     :: nmax
      real(dp),allocatable,intent(inout) :: w2d(:,:) 
      integer :: n,j
      real(dp),allocatable :: w1d(:,:)
      
      allocate(w1d(0:nmax,0:nmax))
      call GetGregW_1D(nmax,w1d)

      if(.not.allocated(w2d)) allocate(w2d(0:nmax,0:nmax))

      w2d = 0.0_dp
      do n=0,nmax
         do j=0,nmax
            w2d(n,j) = w1d(nmax,n)*w1d(n,j)
         end do
      end do

      deallocate(w1d)

   end subroutine GetGregW_2D_Trg
!--------------------------------------------------------------------------------------
 pure real(dp) function DGregoryIntegral(n,h,f)
    integer,parameter :: np=5
    real(dp),dimension(6),parameter :: wg=[&
         0.31559193121693124_dp,&
         1.3921792328042324_dp,&
         0.62397486772486765_dp,&
         1.2440806878306878_dp,&
         0.90990410052910042_dp,&
         1.0142691798941799_dp]
    integer,intent(in)  :: n
    real(dp),intent(in) :: h
    real(dp),intent(in) :: f(0:)
    integer :: i

    DGregoryIntegral = sum(wg(1:np+1)*f(0:np))
    do i=n-np,n
       DGregoryIntegral = DGregoryIntegral + wg(n-i+1)*f(i)
    end do
    DGregoryIntegral = DGregoryIntegral + sum(f(np+1:n-np-1))
    DGregoryIntegral = h*DGregoryIntegral
    
  end function DGregoryIntegral
!--------------------------------------------------------------------------------------
  pure complex(dp) function ZGregoryIntegral(n,h,f)
    integer,parameter :: np=5
    real(dp),dimension(6),parameter :: wg=[&
         0.31559193121693124_dp,&
         1.3921792328042324_dp,&
         0.62397486772486765_dp,&
         1.2440806878306878_dp,&
         0.90990410052910042_dp,&
         1.0142691798941799_dp]
    integer,intent(in)  :: n
    real(dp),intent(in) :: h
    complex(dp),intent(in) :: f(0:)
    integer :: i

    ZGregoryIntegral = sum(wg(1:np+1)*f(0:np))
    do i=n-np,n
       ZGregoryIntegral = ZGregoryIntegral + wg(n-i+1)*f(i)
    end do
    ZGregoryIntegral = ZGregoryIntegral + sum(f(np+1:n-np-1))
    ZGregoryIntegral = h*ZGregoryIntegral
    
  end function ZGregoryIntegral
!--------------------------------------------------------------------------------------
  pure real(dp) function GregoryDot_re_re(n,h,f,g)
    integer,parameter :: np=5
    real(dp),dimension(6),parameter :: wg=[&
         0.31559193121693124_dp,&
         1.3921792328042324_dp,&
         0.62397486772486765_dp,&
         1.2440806878306878_dp,&
         0.90990410052910042_dp,&
         1.0142691798941799_dp]
    integer,intent(in)  :: n
    real(dp),intent(in) :: h
    real(dp),intent(in) :: f(0:)
    real(dp),intent(in) :: g(0:)
    integer :: i
    
    GregoryDot_Re_Re = sum(wg(1:np+1)*f(0:np)*g(0:np))
    do i=n-np,n
       GregoryDot_Re_Re = GregoryDot_Re_Re + wg(n-i+1)*f(i)*g(i)
    end do
    GregoryDot_Re_Re = GregoryDot_Re_Re + dot_product(f(np+1:n-np-1),g(np+1:n-np-1))
    GregoryDot_Re_Re = h*GregoryDot_Re_Re
    
  end function GregoryDot_Re_Re
!--------------------------------------------------------------------------------------
  pure complex(dp) function GregoryDot_c_re(n,h,f,g)
    integer,parameter :: np=5
    real(dp),dimension(6),parameter :: wg=[&
         0.31559193121693124_dp,&
         1.3921792328042324_dp,&
         0.62397486772486765_dp,&
         1.2440806878306878_dp,&
         0.90990410052910042_dp,&
         1.0142691798941799_dp]
    integer,intent(in)     :: n
    real(dp),intent(in)    :: h
    complex(dp),intent(in) :: f(0:)
    real(dp),intent(in)    :: g(0:)
    integer :: i
    
    GregoryDot_c_Re = sum(wg(1:np+1)*conjg(f(0:np))*g(0:np))
    do i=n-np,n
       GregoryDot_c_Re = GregoryDot_c_Re + wg(n-i+1)*conjg(f(i))*g(i)
    end do
    GregoryDot_c_Re = GregoryDot_c_Re + dot_product(f(np+1:n-np-1),g(np+1:n-np-1))
    GregoryDot_c_Re = h*GregoryDot_c_Re
    
  end function GregoryDot_c_Re
!--------------------------------------------------------------------------------------
  pure complex(dp) function GregoryDot_re_c(n,h,f,g)
    integer,parameter :: np=5
    real(dp),dimension(6),parameter :: wg=[&
         0.31559193121693124_dp,&
         1.3921792328042324_dp,&
         0.62397486772486765_dp,&
         1.2440806878306878_dp,&
         0.90990410052910042_dp,&
         1.0142691798941799_dp]
    integer,intent(in)     :: n
    real(dp),intent(in)    :: h
    real(dp),intent(in) :: f(0:)
    complex(dp),intent(in)    :: g(0:)
    integer :: i
    
    GregoryDot_re_c = sum(wg(1:np+1)*f(0:np)*g(0:np))
    do i=n-np,n
       GregoryDot_re_c = GregoryDot_re_c + wg(n-i+1)*f(i)*g(i)
    end do
    GregoryDot_re_c = GregoryDot_re_c + dot_product(f(np+1:n-np-1),g(np+1:n-np-1))
    GregoryDot_re_c = h*GregoryDot_re_c
    
  end function GregoryDot_re_c
!--------------------------------------------------------------------------------------
  pure complex(dp) function GregoryDot_c_c(n,h,f,g)
    integer,parameter :: np=5
    real(dp),dimension(6),parameter :: wg=[&
         0.31559193121693124_dp,&
         1.3921792328042324_dp,&
         0.62397486772486765_dp,&
         1.2440806878306878_dp,&
         0.90990410052910042_dp,&
         1.0142691798941799_dp]
    integer,intent(in)     :: n
    real(dp),intent(in)    :: h
    complex(dp),intent(in) :: f(0:)
    complex(dp),intent(in) :: g(0:)
    integer :: i
    
    GregoryDot_c_c = sum(wg(1:np+1)*conjg(f(0:np))*g(0:np))
    do i=n-np,n
       GregoryDot_c_c = GregoryDot_c_c + wg(n-i+1)*conjg(f(i))*g(i)
    end do
    GregoryDot_c_c = GregoryDot_c_c + dot_product(f(np+1:n-np-1),g(np+1:n-np-1))
    GregoryDot_c_c = h*GregoryDot_c_c
    
  end function GregoryDot_c_c
!--------------------------------------------------------------------------------------


!======================================================================================
end module Mintegration

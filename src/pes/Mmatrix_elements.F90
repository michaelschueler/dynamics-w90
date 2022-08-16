module Mmatrix_elements
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero,one
   use Mquadrature,only: integral_1d
   use Mspecial,only: spherical_bessel_jn
   use Mwignerd,only: ylm_cart
   use Mangcoeff,only: ClebGord,ThreeYlm
   use Mradialwf,only: radialwf_t
   use Mscattwf,only: scattwf_t
   implicit none
   include "../units_inc.f90"
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: MatrixElement_Dipole_Momentum, ScattMatrixElement_Length
   public :: MatrixElement_StructFac
!-------------------------------------------------------------------------------------- 
   real(dp),parameter :: quad_tol=1.0e-12_dp
!--------------------------------------------------------------------------------------   
contains
!--------------------------------------------------------------------------------------
   elemental integer function minusone_n(n)
      integer,intent(in) :: n

      if(mod(n,2) == 0) then
         minusone_n = 1
      else
         minusone_n = -1
      end if

   end function minusone_n
!--------------------------------------------------------------------------------------
   function MatrixElement_StructFac(rwfA,rwfB,l1,m1,l2,m2,qvec) result(Sq)
      real(dp),parameter :: small=1.0e-10_dp
      type(radialwf_t),intent(in)  :: rwfA,rwfB
      integer,intent(in)           :: l1,m1,l2,m2
      real(dp),intent(in)          :: qvec(3)
      complex(dp)                  :: Sq
      integer :: l3,m3,l_min,l_max,bessel_indx
      real(dp) :: qnrm,radint,Rmax
      complex(dp) :: Sqm

      l_min = abs(l1-l2)
      l_max = l1 + l2 
      Rmax = max(rwfA%Rmax, rwfB%Rmax)

      qnrm = norm2(qvec)

      Sq = zero
      do l3=l_min,l_max
         bessel_indx = l3
         call integral_1d(radfunc_bessel,small,Rmax,quad_tol,radint)
         Sqm = zero
         do m3=-l3,l3
            if(m3 + m2 /= m1) cycle
            Sqm = Sqm + conjg(ylm_cart(l3,m3,qvec)) &
               * ThreeYlm(l1,-m1,l3,m3,l2,m2)
         end do
         Sq = Sq + iu**l3 * Sqm * radint
      end do

      if(mod(m1,2) /= 0) Sq = -Sq
      Sq = 2.0_dp * DPI * Sq

      !.................................................
      contains
      !.................................................
      real(dp) function radfunc_bessel(r)
         real(dp),intent(in) :: r
         radfunc_bessel = r**2 * rwfA%Eval_Bessel(bessel_indx,qnrm,r) * rwfB%Eval(r)
      end function radfunc_bessel
      !.................................................

   end function MatrixElement_StructFac
!--------------------------------------------------------------------------------------
   function AngularMatrixElement(l1,m1,l2,m2) result(mel)
      integer,intent(in)  :: l1,m1,l2,m2
      complex(dp) :: mel(3)
      real(dp) :: cff

      mel = 0.0_dp
      if((abs(l1 - l2) > 1) .or. (l1 == l2) .or. (abs(m1 - m2) > 1)) then
         return
      end if

      if(l1 == 0) then
         if(m2 == -1) then
            mel(1) = 1.0_dp/sqrt(6.0_dp)
            mel(2) = -iu/sqrt(6.0_dp)
         elseif(m2 == 1) then
            mel(1) = -1.0_dp/sqrt(6.0_dp)
            mel(2) = -iu/sqrt(6.0_dp)
         else 
            mel(3) = 1.0_dp/sqrt(3.0_dp)
         end if
         mel = - iu * mel 
         return
      end if

      if(l2 == 0) then
         if(m1 == -1) then
            mel(1) = 1.0_dp/sqrt(6.0_dp)
            mel(2) = iu/sqrt(6.0_dp)
         elseif(m1 == 1) then
            mel(1) = -1.0_dp/sqrt(6.0_dp)
            mel(2) = iu/sqrt(6.0_dp)
         else 
            mel(3) = 1.0_dp/sqrt(3.0_dp)
         end if

         mel = -iu * mel
         return
      end if

      if(m2-m1 == 1) then
         cff = ThreeYlm(l1,-m1,1,-1,l2,m2) / sqrt(2.0_dp)
         mel(1) = cff
         mel(2) = iu * cff
      elseif(m2-m1 == -1) then
         cff = ThreeYlm(l1,-m1,1,+1,l2,m2) / sqrt(2.0_dp)
         mel(1) = -cff
         mel(2) = iu * cff
      elseif(m2 == m1) then
         mel(3) = ThreeYlm(l1,-m1,1,0,l2,m2)
      end if

      mel = -iu * sqrt(2.0_dp * Dpi/3.0_dp) * mel
      if(mod(m1,2) /= 0) mel = -mel

   end function AngularMatrixElement
!--------------------------------------------------------------------------------------
   function ScattMatrixElement_Momentum(swf,rwf,l0,m0,kvec) result(Md)
      type(scattwf_t),intent(in)   :: swf
      type(radialwf_t),intent(in)  :: rwf
      integer,intent(in)           :: l0,m0
      real(dp),intent(in)          :: kvec(3)
      complex(dp),dimension(3)     :: Md
      integer :: l,m,l_min,l_max
      integer :: l_indx
      real(dp) :: knrm
      real(dp) :: radint1,radint2
      complex(dp) :: mel_ang(3),exphi
      real(dp),allocatable :: radint(:)

      Md = zero
      l_min = max(l0 - 1, 0)
      l_max = l0 + 1
      knrm = norm2(kvec)

      allocate(radint(lmin:lmax))
      do l=l_min,l_max
         l_indx = l      
         call integral_1d(radfunc1,0.0_dp,rwf%Rmax,quad_tol,radint1)
         call integral_1d(radfunc2,0.0_dp,rwf%Rmax,quad_tol,radint2)         
         radint(l) = radint1 + (0.5_dp* (l*(l+1) - l0*(l0+1)) + 1.0_dp) * radint2
      end do

      do l=l_min,l_max
         exphi = conjg(swf%Phase(l,knrm))
         do m=-l,l
            mel_ang = AngularMatrixElement(l,m,l0,m0)
            Md(1:3) = Md(1:3) + exphi * mel_ang(1:3) * radint(l) * Ylm_cart(l,m,kvec)
         end do
      end do

      Md = QPI * Md
      !.................................................
      contains
      !.................................................
      real(dp) function radfunc1(r)
         real(dp),intent(in) :: r

         radfunc1 = r**2 * swf%Eval(l_indx,knrm,r) * rwf%Eval(r,idx=1)

      end function radfunc1
      !.................................................
      real(dp) function radfunc2(r)
         real(dp),intent(in) :: r

         radfunc2 = r * swf%Eval(l_indx,knrm,r) * rwf%Eval(r)

      end function radfunc2
      !.................................................

   end function ScattMatrixElement_Momentum
!--------------------------------------------------------------------------------------
   function ScattMatrixElement_Length(swf,rwf,l0,m0,kvec) result(Mk)
      real(dp),parameter :: small=1.0e-10_dp
      type(scattwf_t),intent(in)   :: swf
      type(radialwf_t),intent(in)  :: rwf
      integer,intent(in)           :: l0,m0
      real(dp),intent(in)          :: kvec(3)
      complex(dp)                  :: Mk(3)
      integer :: l,l_min,l_max
      integer :: l_indx
      real(dp) :: knrm,rint
      real(dp) :: gnt(3)
      integer :: neval,ier
      real(dp) :: epsabs,epsrel,abserr
      complex(dp) :: exphi

      Mk = zero
      l_min = max(l0 - 1, 0)
      l_max = l0 + 1
      knrm = norm2(kvec)

      ! evaluate 
      epsabs = 1.0e-10_dp; epsrel = 1.0e-10_dp
      do l=l_min,l_max
         l_indx = l
         call integral_1d(radfunc,0.0_dp,rwf%Rmax,quad_tol,rint)         
         gnt(1) = ThreeYlm(l,-m0+1,1,-1,l0,m0)
         gnt(2) = ThreeYlm(l,-m0-1,1,1,l0,m0)
         gnt(3) = ThreeYlm(l,-m0,1,0,l0,m0)
         exphi = conjg(swf%Phase(l,knrm))

         Mk(1) = Mk(1) + exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
            - gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint / sqrt(2.0d0)
         Mk(2) = Mk(2) + iu*exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
            + gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint / sqrt(2.0d0)
         Mk(3) = Mk(3) + exphi * gnt(3) * rint * conjg(Ylm_cart(l,-m0,kvec))
      end do

      Mk = QPI * sqrt(QPI/3.0d0) * Mk

      !.................................................
      contains
      !.................................................
      real(dp) function radfunc(r)
         real(dp),intent(in) :: r

         radfunc = r**3 * rwf%Eval(r) * swf%Eval(l_indx, knrm, r)

      end function radfunc
      !.................................................

   end function ScattMatrixElement_Length
!-------------------------------------------------------------------------------------- 

!--------------------------------------------------------------------------------------
   function MatrixElement_Dipole_finiteQ(rwfA,rwfB,l1,m1,l2,m2,qvec) result(Md)
      real(dp),parameter :: small=1.0e-10_dp
      type(radialwf_t),intent(in)  :: rwfA,rwfB
      integer,intent(in)           :: l1,m1,l2,m2
      real(dp),intent(in)          :: qvec(3)
      complex(dp),dimension(3)     :: Md
      integer :: l3,m3,lam_min,lam_max,l_min,l_max,lam,mu
      integer :: bessel_indx
      real(dp) :: qnrm,rint,Rmax
      real(dp),allocatable :: rint1(:),rint2(:)
      complex(dp) :: yfac,Mdq(3),mel_ang(3)
      real(dp) :: tol

      Md = zero
      qnrm = norm2(qvec)

      if(qnrm < small) then
         Md = MatrixElement_Dipole(rwfA,rwfB,l1,m1,l2,m2)
         return 
      end if

      if(l2 == 0) then
         l_min = 0
      else
         l_min = l2 - 1
      end if
      l_max = l2 + 1

      Rmax = max(rwfA%Rmax, rwfB%Rmax)
      
      ! precompute radial integrals
      lam_max = l1 + l2 + 2
      allocate(rint1(0:lam_max),rint2(0:lam_max))
      do lam=0,lam_max
         bessel_indx = lam
         call integral_1d(radfunc1_bessel,0.0_dp,Rmax,quad_tol,rint1(lam))
         call integral_1d(radfunc2_bessel,0.0_dp,Rmax,quad_tol,rint2(lam))
      end do

      ! evaluate 
      do l3=l_min,l_max
         lam_min = abs(l3-l1); lam_max = l3+l1
         do lam=lam_min,lam_max
            rint = rint1(lam) + (0.5_dp* (l2*(l2+1) - l3*(l3+1)) + 1.0_dp) * rint2(lam)
            do m3=-l3,l3
               mel_ang = AngularMatrixElement(l3,m3,l2,m2)
               do mu=-lam,lam
                  yfac = iu**lam * conjg(ylm_cart(lam,mu,qvec))
                   Md(1:3) = Md(1:3) + yfac * minusone_n(mu + m3) &
                     * ThreeYlm(l3,-m3,lam,-mu,l1,m1) * rint * mel_ang
               end do
            end do
         end do
      end do

      Md = 2.0_dp * DPI * Md

      deallocate(rint1,rint2)

      !.................................................
      contains
      !.................................................
      real(dp) function radfunc1_bessel(r)
         real(dp),intent(in) :: r

         radfunc1_bessel = r**2 * rwfA%Eval_Bessel(bessel_indx,qnrm,r) * rwfB%Eval(r,idx=1)

      end function radfunc1_bessel
      !.................................................
      real(dp) function radfunc2_bessel(r)
         real(dp),intent(in) :: r

         radfunc2_bessel = r * rwfA%Eval_Bessel(bessel_indx,qnrm,r) * rwfB%Eval(r)

      end function radfunc2_bessel
      !.................................................

   end function MatrixElement_Dipole_finiteQ
!--------------------------------------------------------------------------------------
   function MatrixElement_PWOVLP(rwf,l0,m0,kvec) result(Sk)
      real(dp),parameter :: small=1.0e-10_dp
      type(radialwf_t),intent(in)  :: rwf
      integer,intent(in)           :: l0,m0
      real(dp),intent(in)          :: kvec(3)
      complex(dp)                  :: Sk
      integer :: bessel_indx
      real(dp) :: knrm,rint

      knrm = norm2(kvec)
      bessel_indx = l0

      call integral_1d(radfunc_bessel,0.0_dp,rwf%Rmax,quad_tol,rint)
      Sk = 2.0_dp * DPI * (-iu)**l0 * ylm_cart(l0,m0,kvec) * rint

      !.................................................
      contains
      !.................................................
      real(dp) function radfunc_bessel(r)
         real(dp),intent(in) :: r

         radfunc_bessel = r**2 * rwf%Eval_Bessel(bessel_indx,knrm,r)

      end function radfunc_bessel
      !.................................................

   end function MatrixElement_PWOVLP
!--------------------------------------------------------------------------------------
   function MatrixElement_PWDIP(rwf,l0,m0,kvec) result(Mk)
      real(dp),parameter :: small=1.0e-10_dp
      type(radialwf_t),intent(in)  :: rwf
      integer,intent(in)           :: l0,m0
      real(dp),intent(in)          :: kvec(3)
      complex(dp)                  :: Mk(3)
      integer :: l,l_min,l_max
      integer :: bessel_indx
      real(dp) :: knrm,rint
      real(dp) :: gnt(3)
      integer :: neval,ier
      real(dp) :: epsabs,epsrel,abserr

      Mk = zero
      l_min = max(l0 - 1, 0)
      l_max = l0 + 1
      knrm = norm2(kvec)

      ! evaluate 
      epsabs = 1.0e-10_dp; epsrel = 1.0e-10_dp
      do l=l_min,l_max
         bessel_indx = l
         call integral_1d(radfunc_bessel,0.0_dp,rwf%Rmax,quad_tol,rint)
         gnt(1) = ThreeYlm(l,-m0+1,1,-1,l0,m0)
         gnt(2) = ThreeYlm(l,-m0-1,1,1,l0,m0)
         gnt(3) = ThreeYlm(l,-m0,1,0,l0,m0)
         Mk(1) = Mk(1) + (-iu)**l * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
            - gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint / sqrt(2.0d0)
         Mk(2) = Mk(2) + iu*(-iu)**l * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
            + gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint / sqrt(2.0d0)
         Mk(3) = Mk(3) + (-iu)**l * gnt(3) * rint * conjg(Ylm_cart(l,-m0,kvec))
      end do

      Mk = QPI * sqrt(QPI/3.0d0) * Mk

      !.................................................
      contains
      !.................................................
      real(dp) function radfunc_bessel(r)
         real(dp),intent(in) :: r

         radfunc_bessel = r**3 * rwf%Eval_Bessel(bessel_indx,knrm,r)

      end function radfunc_bessel
      !.................................................

   end function MatrixElement_PWDIP
!--------------------------------------------------------------------------------------
   function MatrixElement_Dipole_cart(Rmax,wfA,wfB) result(Md)
      use Mbsplines,only: spline3d_t
      use Mutils,only: linspace
      use Mintegration,only: GetGregoryWeights
      integer,parameter :: ksp=4
      real(dp),intent(in)      :: Rmax
      complex(dp),intent(in)   :: wfA(:,:,:)
      complex(dp),intent(in)   :: wfB(:,:,:)
      complex(dp),dimension(3) :: Md
      integer :: iflag,inbvx(2),inbvy(2),inbvz(2),iloy(2),iloz(2)
      integer  :: nx,ny,nz,ix,iy,iz
      real(dp) :: dx,dy,dz
      real(dp),allocatable :: xs(:),ys(:),zs(:)
      real(dp),allocatable :: wx(:),wy(:),wz(:)
      complex(dp) :: dwfBd
      type(spline3d_t) :: wfB_spl_real,wfB_spl_imag

      nx = size(wfA,1)
      ny = size(wfA,2)
      nz = size(wfA,3)

      xs = linspace(-Rmax, Rmax, nx)
      ys = linspace(-Rmax, Rmax, ny)
      zs = linspace(-Rmax, Rmax, nz)
      dx = xs(2) - xs(1)
      dy = ys(2) - ys(1)
      dz = zs(2) - zs(1)

      allocate(wx(0:nx-1),wy(0:ny-1),wz(0:nz-1))
      call GetGregoryWeights(nx-1,wx)
      call GetGregoryWeights(ny-1,wy)
      call GetGregoryWeights(nz-1,wz)

      iflag = 0
      call wfB_spl_real%Init(xs,ys,zs,dble(wfB),ksp,ksp,ksp,iflag)
      call wfB_spl_imag%Init(xs,ys,zs,aimag(wfB),ksp,ksp,ksp,iflag)

      Md = zero

      !$OMP PARALLEL PRIVATE(iflag,inbvx,inbvy,inbvz,iloy,iloz,dwfBd)

      inbvx=1; inbvy=1; inbvz=1; iloy=1; iloz=1

      !$OMP DO REDUCTION(+:Md) COLLAPSE(3)
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               dwfBd = wfB_spl_real%Eval(xs(ix),ys(iy),zs(iz),1,0,0,&
                  iflag,inbvx(1),inbvy(1),inbvz(1),iloy(1),iloz(1))
               dwfBd = dwfBd + iu * wfB_spl_imag%Eval(xs(ix),ys(iy),zs(iz),1,0,0,&
                  iflag,inbvx(2),inbvy(2),inbvz(2),iloy(2),iloz(2))               
               Md(1) = Md(1) - iu * wx(ix-1)*wy(iy-1)*wz(iz-1) * conjg(wfA(ix,iy,iz)) * dwfBd
            end do
         end do
      end do
      !$OMP END DO

      inbvx=1; inbvy=1; inbvz=1; iloy=1; iloz=1

      !$OMP DO REDUCTION(+:Md) COLLAPSE(3)
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               dwfBd = wfB_spl_real%Eval(xs(ix),ys(iy),zs(iz),0,1,0,&
                  iflag,inbvx(1),inbvy(1),inbvz(1),iloy(1),iloz(1))
               dwfBd = dwfBd + iu * wfB_spl_imag%Eval(xs(ix),ys(iy),zs(iz),0,1,0,&
                  iflag,inbvx(2),inbvy(2),inbvz(2),iloy(2),iloz(2))               
               Md(2) = Md(2) - iu * wx(ix-1)*wy(iy-1)*wz(iz-1)* conjg(wfA(ix,iy,iz)) * dwfBd
            end do
         end do
      end do
      !$OMP END DO

      inbvx=1; inbvy=1; inbvz=1; iloy=1; iloz=1

      !$OMP DO REDUCTION(+:Md) COLLAPSE(3)
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               dwfBd = wfB_spl_real%Eval(xs(ix),ys(iy),zs(iz),0,0,1,&
                  iflag,inbvx(1),inbvy(1),inbvz(1),iloy(1),iloz(1))
               dwfBd = dwfBd + iu * wfB_spl_imag%Eval(xs(ix),ys(iy),zs(iz),0,0,1,&
                  iflag,inbvx(2),inbvy(2),inbvz(2),iloy(2),iloz(2))               
               Md(3) = Md(3) - iu * wx(ix-1)*wy(iy-1)*wz(iz-1)*conjg(wfA(ix,iy,iz)) * dwfBd
            end do
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Md = dx * dy * dz * Md

      deallocate(xs, ys, zs, wx, wy, wz)
      call wfB_spl_real%Clean()
      call wfB_spl_imag%Clean()

   end function MatrixElement_Dipole_cart
!--------------------------------------------------------------------------------------
   function MatrixElement_Dipole_finiteQ_cart(Rmax,wfA,wfB,qvec) result(Md)
      use Mbsplines,only: spline3d_t
      use Mutils,only: linspace
      use Mintegration,only: GetGregoryWeights
      integer,parameter :: ksp=4
      real(dp),intent(in)      :: Rmax
      complex(dp),intent(in)   :: wfA(:,:,:)
      complex(dp),intent(in)   :: wfB(:,:,:)
      real(dp),intent(in)      :: qvec(:)
      complex(dp),dimension(3) :: Md
      integer :: iflag,inbvx(2),inbvy(2),inbvz(2),iloy(2),iloz(2)
      integer  :: nx,ny,nz,ix,iy,iz
      real(dp) :: dx,dy,dz
      real(dp),allocatable :: xs(:),ys(:),zs(:)
      real(dp),allocatable :: wx(:),wy(:),wz(:)
      complex(dp) :: dwfBd
      type(spline3d_t) :: wfB_spl_real,wfB_spl_imag

      nx = size(wfA,1)
      ny = size(wfA,2)
      nz = size(wfA,3)

      xs = linspace(-Rmax, Rmax, nx)
      ys = linspace(-Rmax, Rmax, ny)
      zs = linspace(-Rmax, Rmax, nz)
      dx = xs(2) - xs(1)
      dy = ys(2) - ys(1)
      dz = zs(2) - zs(1)

      allocate(wx(0:nx-1),wy(0:ny-1),wz(0:nz-1))
      call GetGregoryWeights(nx-1,wx)
      call GetGregoryWeights(ny-1,wy)
      call GetGregoryWeights(nz-1,wz)

      iflag = 0
      call wfB_spl_real%Init(xs,ys,zs,dble(wfB),ksp,ksp,ksp,iflag)
      call wfB_spl_imag%Init(xs,ys,zs,aimag(wfB),ksp,ksp,ksp,iflag)

      Md = zero

      !$OMP PARALLEL PRIVATE(iflag,inbvx,inbvy,inbvz,iloy,iloz,dwfBd)

      inbvx=1; inbvy=1; inbvz=1; iloy=1; iloz=1

      !$OMP DO REDUCTION(+:Md) COLLAPSE(3)
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               dwfBd = wfB_spl_real%Eval(xs(ix),ys(iy),zs(iz),1,0,0,&
                  iflag,inbvx(1),inbvy(1),inbvz(1),iloy(1),iloz(1))
               dwfBd = dwfBd + iu * wfB_spl_imag%Eval(xs(ix),ys(iy),zs(iz),1,0,0,&
                  iflag,inbvx(2),inbvy(2),inbvz(2),iloy(2),iloz(2))               
               Md(1) = Md(1) - iu * wx(ix-1)*wy(iy-1)*wz(iz-1) * conjg(wfA(ix,iy,iz)) &
                  * dwfBd * exp(iu*(qvec(1)*xs(ix) + qvec(2)*ys(iy) + qvec(3)*zs(iz) ))
            end do
         end do
      end do
      !$OMP END DO

      inbvx=1; inbvy=1; inbvz=1; iloy=1; iloz=1

      !$OMP DO REDUCTION(+:Md) COLLAPSE(3)
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               dwfBd = wfB_spl_real%Eval(xs(ix),ys(iy),zs(iz),0,1,0,&
                  iflag,inbvx(1),inbvy(1),inbvz(1),iloy(1),iloz(1))
               dwfBd = dwfBd + iu * wfB_spl_imag%Eval(xs(ix),ys(iy),zs(iz),0,1,0,&
                  iflag,inbvx(2),inbvy(2),inbvz(2),iloy(2),iloz(2))               
               Md(2) = Md(2) - iu * wx(ix-1)*wy(iy-1)*wz(iz-1)* conjg(wfA(ix,iy,iz)) &
                  * dwfBd * exp(iu*(qvec(1)*xs(ix) + qvec(2)*ys(iy) + qvec(3)*zs(iz) ))
            end do
         end do
      end do
      !$OMP END DO

      inbvx=1; inbvy=1; inbvz=1; iloy=1; iloz=1

      !$OMP DO REDUCTION(+:Md) COLLAPSE(3)
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               dwfBd = wfB_spl_real%Eval(xs(ix),ys(iy),zs(iz),0,0,1,&
                  iflag,inbvx(1),inbvy(1),inbvz(1),iloy(1),iloz(1))
               dwfBd = dwfBd + iu * wfB_spl_imag%Eval(xs(ix),ys(iy),zs(iz),0,0,1,&
                  iflag,inbvx(2),inbvy(2),inbvz(2),iloy(2),iloz(2))               
               Md(3) = Md(3) - iu * wx(ix-1)*wy(iy-1)*wz(iz-1)*conjg(wfA(ix,iy,iz)) &
                  * dwfBd * exp(iu*(qvec(1)*xs(ix) + qvec(2)*ys(iy) + qvec(3)*zs(iz) ))
            end do
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Md = dx * dy * dz * Md

      deallocate(xs, ys, zs, wx, wy, wz)
      call wfB_spl_real%Clean()
      call wfB_spl_imag%Clean()

   end function MatrixElement_Dipole_finiteQ_cart
!--------------------------------------------------------------------------------------
   function MatrixElement_StructFac_cart(Rmax,wfA,wfB,qvec) result(Sq)
      use Mutils,only: linspace
      use Mintegration,only: GetGregoryWeights
      real(dp),intent(in)      :: Rmax
      complex(dp),intent(in)   :: wfA(:,:,:)
      complex(dp),intent(in)   :: wfB(:,:,:)
      real(dp),intent(in)      :: qvec(:)
      complex(dp)              :: Sq
      integer  :: nx,ny,nz,ix,iy,iz
      real(dp) :: dx,dy,dz,qdotr
      real(dp),allocatable :: xs(:),ys(:),zs(:)
      real(dp),allocatable :: wx(:),wy(:),wz(:)

      nx = size(wfA,1)
      ny = size(wfA,2)
      nz = size(wfA,3)

      call assert(size(wfB,1) == nx)
      call assert(size(wfB,2) == ny)
      call assert(size(wfB,3) == nz)

      xs = linspace(-Rmax, Rmax, nx)
      ys = linspace(-Rmax, Rmax, ny)
      zs = linspace(-Rmax, Rmax, nz)
      dx = xs(2) - xs(1)
      dy = ys(2) - ys(1)
      dz = zs(2) - zs(1)

      allocate(wx(0:nx-1),wy(0:ny-1),wz(0:nz-1))
      call GetGregoryWeights(nx-1,wx)
      call GetGregoryWeights(ny-1,wy)
      call GetGregoryWeights(nz-1,wz)

      Sq = zero

      !$OMP PARALLEL PRIVATE(qdotr)
      !$OMP DO REDUCTION(+:Sq) COLLAPSE(3)
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               qdotr = qvec(1)*xs(ix) + qvec(2)*ys(iy) + qvec(3)*zs(iz) 
               Sq = Sq + wx(ix-1)*wy(iy-1)*wz(iz-1) * conjg(wfA(ix,iy,iz)) &
                  * wfB(ix,iy,iz) * exp(iu*qdotr)
            end do
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Sq = dx * dy * dz * Sq

      deallocate(xs, ys, zs, wx, wy, wz)

   end function MatrixElement_StructFac_cart
!--------------------------------------------------------------------------------------
   function MatrixElement_PWOVLP_cart(Rmax,wf,kvec) result(Sk)
      use Mutils,only: linspace
      use Mintegration,only: GetGregoryWeights
      real(dp),intent(in)      :: Rmax
      complex(dp),intent(in)   :: wf(:,:,:)
      real(dp),intent(in)      :: kvec(:)
      complex(dp)              :: Sk
      integer  :: nx,ny,nz,ix,iy,iz
      real(dp) :: dx,dy,dz,kdotr
      real(dp),allocatable :: xs(:),ys(:),zs(:)
      real(dp),allocatable :: wx(:),wy(:),wz(:)

      nx = size(wf,1)
      ny = size(wf,2)
      nz = size(wf,3)

      xs = linspace(-Rmax, Rmax, nx)
      ys = linspace(-Rmax, Rmax, ny)
      zs = linspace(-Rmax, Rmax, nz)
      dx = xs(2) - xs(1)
      dy = ys(2) - ys(1)
      dz = zs(2) - zs(1)

      allocate(wx(0:nx-1),wy(0:ny-1),wz(0:nz-1))
      call GetGregoryWeights(nx-1,wx)
      call GetGregoryWeights(ny-1,wy)
      call GetGregoryWeights(nz-1,wz)

      Sk = zero

      !$OMP PARALLEL PRIVATE(kdotr)
      !$OMP DO REDUCTION(+:Sk) COLLAPSE(3)
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               kdotr = kvec(1)*xs(ix) + kvec(2)*ys(iy) + kvec(3)*zs(iz) 
               Sk = Sk + wx(ix-1)*wy(iy-1)*wz(iz-1) * exp(-iu*kdotr) * wf(ix,iy,iz)
            end do
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Sk = dx * dy * dz * Sk

      deallocate(xs, ys, zs, wx, wy, wz)

   end function MatrixElement_PWOVLP_cart
!--------------------------------------------------------------------------------------
   function MatrixElement_PWDIP_cart(Rmax,wf,kvec) result(Mk)
      use Mutils,only: linspace
      use Mintegration,only: GetGregoryWeights
      real(dp),intent(in)      :: Rmax
      complex(dp),intent(in)   :: wf(:,:,:)
      real(dp),intent(in)      :: kvec(:)
      complex(dp)              :: Mk(3)
      integer  :: nx,ny,nz,ix,iy,iz
      real(dp) :: dx,dy,dz,kdotr
      real(dp),allocatable :: xs(:),ys(:),zs(:)
      real(dp),allocatable :: wx(:),wy(:),wz(:)

      nx = size(wf,1)
      ny = size(wf,2)
      nz = size(wf,3)

      xs = linspace(-Rmax, Rmax, nx)
      ys = linspace(-Rmax, Rmax, ny)
      zs = linspace(-Rmax, Rmax, nz)
      dx = xs(2) - xs(1)
      dy = ys(2) - ys(1)
      dz = zs(2) - zs(1)

      allocate(wx(0:nx-1),wy(0:ny-1),wz(0:nz-1))
      call GetGregoryWeights(nx-1,wx)
      call GetGregoryWeights(ny-1,wy)
      call GetGregoryWeights(nz-1,wz)

      Mk = zero

      !$OMP PARALLEL PRIVATE(kdotr)
      !$OMP DO REDUCTION(+:Mk) COLLAPSE(3)
      do iz=1,nz
         do iy=1,ny
            do ix=1,nx
               kdotr = kvec(1)*xs(ix) + kvec(2)*ys(iy) + kvec(3)*zs(iz) 
               Mk(1) = Mk(1) + wx(ix-1)*wy(iy-1)*wz(iz-1) * exp(-iu*kdotr) &
                  * wf(ix,iy,iz) * xs(ix)
               Mk(2) = Mk(2) + wx(ix-1)*wy(iy-1)*wz(iz-1) * exp(-iu*kdotr) &
                  * wf(ix,iy,iz) * ys(iy)
               Mk(3) = Mk(3) + wx(ix-1)*wy(iy-1)*wz(iz-1) * exp(-iu*kdotr) &
                  * wf(ix,iy,iz) * zs(iz)
            end do
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Mk = dx * dy * dz * Mk

      deallocate(xs, ys, zs, wx, wy, wz)

   end function MatrixElement_PWDIP_cart
!--------------------------------------------------------------------------------------
  recursive function lacz_gamma(a) result(g)
    double complex, intent(in) :: a
    double complex :: g
 
    double precision, parameter :: sq2p=sqrt(0.5D0*Qpi),pi=0.25D0*QPi
    integer, parameter :: cg = 7
 
    ! these precomputed values are taken by the sample code in Wikipedia,
    ! and the sample itself takes them from the GNU Scientific Library
    double precision, dimension(0:8), parameter :: p = &
         (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
         771.32342877765313, -176.61502916214059, 12.507343278686905, &
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)
 
    double complex :: t, w, x
    integer :: i
 
    x = a
 
    if ( dreal(x) < 0.5D0 ) then
       g = pi/( sin(pi*x) * lacz_gamma(1D0-x) )
    else
       x = x - 1D0
       t = p(0)
       do i=1, cg+1
          t = t + p(i)/(x+real(i,8))
       end do
       w = x + real(cg,8) + 0.5D0
       g = sq2p* w**(x+0.5D0) * exp(-w) * t
    end if
  end function lacz_gamma
!--------------------------------------------------------------------------------------


!======================================================================================
end module Mmatrix_elements
module Mmatrix_elements
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero,one
   use Mutils,only: linspace
   use Mquadrature,only: integral_1d
   use Mspecial,only: spherical_bessel_jn
   use Mlebedev_quad,only: Lebedev_integral
   use Mwignerd,only: ylm,ylm_cart
   use Mangcoeff,only: ClebGord,ThreeYlm
   use Mradialwf,only: radialwf_t
   use Mvector_bsplines,only: cplx_vector_spline_t, cplx_matrix_spline_t
   use Mscattwf,only: scattwf_t
   use Mradialintegral,only: radialinteg_t
   implicit none
   include "../units_inc.f90"
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: ScattMatrixElement_Momentum, ScattMatrixElement_Length
   public :: ApplyDipoleOp, ApplyMomentumOp, ApplyLengthOp
   public :: dipole_lambda_projection, dipole_lambda_BesselTransform, ScattMatrixElement_Lambda
#ifdef MPI
   public :: dipole_lambda_BesselTransform_mpi
#endif

   interface ScattMatrixElement_Momentum
      module procedure ScattMatrixElement_Momentum_comp, ScattMatrixElement_Momentum_precomp
   end interface ScattMatrixElement_Momentum

   interface ScattMatrixElement_Length
      module procedure ScattMatrixElement_Length_comp, ScattMatrixElement_Length_precomp
   end interface ScattMatrixElement_Length  
!-------------------------------------------------------------------------------------- 
   real(dp),parameter :: quad_tol=1.0e-8_dp
   integer,parameter :: wf_slater=0, wf_grid=1
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
   function ScattMatrixElement_Momentum_comp(swf,rwf,l0,m0,kvec) result(Md)
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

      allocate(radint(l_min:l_max))
      do l=l_min,l_max
         if(l == l0) cycle
         l_indx = l      
         call integral_1d(radfunc1,0.0_dp,rwf%Rmax,quad_tol,radint1)
         call integral_1d(radfunc2,0.0_dp,rwf%Rmax,quad_tol,radint2)         
         radint(l) = radint1 + (0.5_dp* (l0*(l0+1) - l*(l+1)) + 1.0_dp) * radint2
      end do

      do l=l_min,l_max
         if(l == l0) cycle
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

   end function ScattMatrixElement_Momentum_comp
!--------------------------------------------------------------------------------------
   function ScattMatrixElement_Momentum_precomp(swf,radialinteg,l0,m0,kvec) result(Md)
      type(scattwf_t),intent(in)        :: swf
      type(radialinteg_t),intent(in) :: radialinteg
      integer,intent(in)                :: l0,m0
      real(dp),intent(in)               :: kvec(3)
      complex(dp),dimension(3)          :: Md
      integer :: l,m
      real(dp) :: knrm,rint(2)
      complex(dp) :: mel_ang(3),exphi

      Md = zero
      knrm = norm2(kvec)

      call radialinteg%Eval_mom(knrm,rint)

      l = l0 - 1
      exphi = conjg(swf%Phase(l,knrm))
      if(l >= 0) then
         do m=-l,l
            mel_ang = AngularMatrixElement(l,m,l0,m0)
            Md(1:3) = Md(1:3) + exphi * mel_ang(1:3) * rint(1) * Ylm_cart(l,m,kvec)
         end do         
      end if

      l = l0 + 1
      exphi = conjg(swf%Phase(l,knrm))
      do m=-l,l
         mel_ang = AngularMatrixElement(l,m,l0,m0)
         Md(1:3) = Md(1:3) + exphi * mel_ang(1:3) * rint(2) * Ylm_cart(l,m,kvec)
      end do          

      Md = QPI * Md


   end function ScattMatrixElement_Momentum_precomp
!--------------------------------------------------------------------------------------
   function ScattMatrixElement_Length_comp(swf,rwf,l0,m0,kvec) result(Mk)
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
         if(l == l0) cycle
         l_indx = l
         call integral_1d(radfunc,0.0_dp,rwf%Rmax,quad_tol,rint)       
         gnt(1) = minusone_n(-m0+1) * ThreeYlm(l,-m0+1,1,-1,l0,m0)  
         gnt(2) = minusone_n(-m0-1) * ThreeYlm(l,-m0-1,1,1,l0,m0)
         gnt(3) = minusone_n(-m0) * ThreeYlm(l,-m0,1,0,l0,m0)
         exphi = conjg(swf%Phase(l,knrm))

         Mk(1) = Mk(1) + exphi * rint * (gnt(1) * Ylm_cart(l,-m0+1,kvec) &
            - gnt(2) * Ylm_cart(l,-m0-1,kvec)) / sqrt(2.0d0)

         Mk(2) = Mk(2) + iu * exphi * rint * (gnt(1) * Ylm_cart(l,-m0+1,kvec) &
            + gnt(2) * Ylm_cart(l,-m0-1,kvec)) / sqrt(2.0d0)

         Mk(3) = Mk(3) + exphi * rint * gnt(3) * Ylm_cart(l,-m0,kvec)

         ! Mk(1) = Mk(1) + exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
         !    - gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint / sqrt(2.0d0)
         ! Mk(2) = Mk(2) + iu*exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
         !    + gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint / sqrt(2.0d0)
         ! Mk(3) = Mk(3) + exphi * gnt(3) * rint * conjg(Ylm_cart(l,-m0,kvec))
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

   end function ScattMatrixElement_Length_comp
!-------------------------------------------------------------------------------------- 
   function ScattMatrixElement_Length_precomp(swf,radialinteg,l0,m0,kvec) result(Mk)
      type(scattwf_t),intent(in)        :: swf
      type(radialinteg_t),intent(in)    :: radialinteg
      integer,intent(in)                :: l0,m0
      real(dp),intent(in)               :: kvec(3)
      complex(dp)                       :: Mk(3)
      integer :: l
      real(dp) :: knrm,rint(2)
      real(dp) :: gnt(3)
      complex(dp) :: exphi

      Mk = zero
      knrm = norm2(kvec)

      call radialinteg%Eval_len(knrm,rint)

      l = l0 - 1
      if(l >= 0) then
         ! gnt(1) = ThreeYlm(l,-m0+1,1,-1,l0,m0)
         ! gnt(2) = ThreeYlm(l,-m0-1,1,1,l0,m0)
         ! gnt(3) = ThreeYlm(l,-m0,1,0,l0,m0)
         gnt(1) = minusone_n(-m0+1) * ThreeYlm(l,-m0+1,1,-1,l0,m0)  
         gnt(2) = minusone_n(-m0-1) * ThreeYlm(l,-m0-1,1,1,l0,m0)
         gnt(3) = minusone_n(-m0) * ThreeYlm(l,-m0,1,0,l0,m0)
         exphi = conjg(swf%Phase(l,knrm))

         Mk(1) = Mk(1) + exphi * rint(1) * (gnt(1) * Ylm_cart(l,m0-1,kvec) &
            - gnt(2) * Ylm_cart(l,m0+1,kvec)) / sqrt(2.0d0)

         Mk(2) = Mk(2) + iu * exphi * rint(1) * (gnt(1) * Ylm_cart(l,m0-1,kvec) &
            + gnt(2) * Ylm_cart(l,m0+1,kvec)) / sqrt(2.0d0)

         Mk(3) = Mk(3) + exphi * rint(1) * gnt(3) * Ylm_cart(l,m0,kvec)

         ! Mk(1) = Mk(1) + exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
         !    - gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint(1) / sqrt(2.0d0)
         ! Mk(2) = Mk(2) + iu*exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
         !    + gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint(1) / sqrt(2.0d0)
         ! Mk(3) = Mk(3) + exphi * gnt(3) * rint(1) * conjg(Ylm_cart(l,-m0,kvec))
      end if

      l = l0 + 1
      ! gnt(1) = ThreeYlm(l,-m0+1,1,-1,l0,m0)
      ! gnt(2) = ThreeYlm(l,-m0-1,1,1,l0,m0)
      ! gnt(3) = ThreeYlm(l,-m0,1,0,l0,m0)
      gnt(1) = minusone_n(-m0+1) * ThreeYlm(l,-m0+1,1,-1,l0,m0)  
      gnt(2) = minusone_n(-m0-1) * ThreeYlm(l,-m0-1,1,1,l0,m0)
      gnt(3) = minusone_n(-m0) * ThreeYlm(l,-m0,1,0,l0,m0)
      exphi = conjg(swf%Phase(l,knrm))

      exphi = conjg(swf%Phase(l,knrm))

      ! Mk(1) = Mk(1) + exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
      !    - gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint(2) / sqrt(2.0d0)
      ! Mk(2) = Mk(2) + iu*exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
      !    + gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint(2) / sqrt(2.0d0)
      ! Mk(3) = Mk(3) + exphi * gnt(3) * rint(2) * conjg(Ylm_cart(l,-m0,kvec))              

      Mk(1) = Mk(1) + exphi * rint(2) * (gnt(1) * Ylm_cart(l,m0-1,kvec) &
         - gnt(2) * Ylm_cart(l,m0+1,kvec)) / sqrt(2.0d0)

      Mk(2) = Mk(2) + iu * exphi * rint(2) * (gnt(1) * Ylm_cart(l,m0-1,kvec) &
         + gnt(2) * Ylm_cart(l,m0+1,kvec)) / sqrt(2.0d0)

      Mk(3) = Mk(3) + exphi * rint(2) * gnt(3) * Ylm_cart(l,m0,kvec)

      Mk = QPI * sqrt(QPI/3.0d0) * Mk

   end function ScattMatrixElement_Length_precomp
!-------------------------------------------------------------------------------------- 
   subroutine ApplyDipoleOp(rwf,l0,m0,dipwf_m1,dipwf_p1,gauge)
      integer,parameter :: gauge_len=0, gauge_mom=1
      type(radialwf_t),intent(in)            :: rwf
      integer,intent(in)                     :: l0,m0
      type(cplx_matrix_spline_t),intent(out) :: dipwf_m1,dipwf_p1
      integer,intent(in)                     :: gauge

      select case(gauge)
      case(gauge_len)
         call ApplyLengthOp(rwf,l0,m0,dipwf_m1,dipwf_p1)
      case(gauge_mom)
         call ApplyMomentumOp(rwf,l0,m0,dipwf_m1,dipwf_p1)
      end select

   end subroutine ApplyDipoleOp
!-------------------------------------------------------------------------------------- 
   subroutine ApplyLengthOp(rwf,l0,m0,dipwf_m1,dipwf_p1)
      integer,parameter :: Nr=400
      type(radialwf_t),intent(in)            :: rwf
      integer,intent(in)                     :: l0,m0
      type(cplx_matrix_spline_t),intent(out) :: dipwf_m1,dipwf_p1
      integer :: l,m,i,ir,sig
      real(dp) :: rmax,rv
      real(dp) :: gnt(3)
      real(dp),allocatable :: rs(:)
      complex(dp),allocatable :: Rfun(:,:,:)

      Rmax = rwf%Rmax
      rs = linspace(0.0_dp, Rmax, Nr)

      !.......................................
      !            l = l0 + 1
      !.......................................
      l = l0 + 1

      gnt(1) = ThreeYlm(l,-m0+1,1,-1,l0,m0)
      gnt(2) = ThreeYlm(l,-m0-1,1,1,l0,m0)
      gnt(3) = ThreeYlm(l,-m0,1,0,l0,m0)      

      allocate(Rfun(Nr,2*l+1,3)); Rfun = zero
      do m= -l, l
         i = m + l + 1
         sig = minusone_n(m)
         if(m == m0-1) then
            do ir=1,Nr
               rv = rs(ir) * rwf%Eval(rs(ir))
               Rfun(ir,i,1) = Rfun(ir,i,1) + sig * gnt(1) * rv / sqrt(2.0d0)
               Rfun(ir,i,2) = Rfun(ir,i,2) + sig * iu * gnt(1) * rv / sqrt(2.0d0)
            end do
         elseif(m == m0 + 1) then
            do ir=1,Nr
               rv = rs(ir) * rwf%Eval(rs(ir))
               Rfun(ir,i,1) = Rfun(ir,i,1) - sig * gnt(2) * rv / sqrt(2.0d0)
               Rfun(ir,i,2) = Rfun(ir,i,2) + sig * iu * gnt(2) * rv / sqrt(2.0d0)
            end do            
         elseif(m == m0) then
            do ir=1,Nr
               rv = rs(ir) * rwf%Eval(rs(ir))
               Rfun(ir,i,3) = Rfun(ir,i,3) + sig * gnt(3) * rv
            end do
         end if
      end do
      Rfun = sqrt(QPI/3.0d0) * Rfun 

      call dipwf_p1%Init(rs, Rfun, 2*l+1, 3)

      deallocate(Rfun)

      !.......................................
      !            l = l0 - 1
      !.......................................
      l = l0 - 1

      if(l0 > 0) then
         gnt(1) = ThreeYlm(l,-m0+1,1,-1,l0,m0)
         gnt(2) = ThreeYlm(l,-m0-1,1,1,l0,m0)
         gnt(3) = ThreeYlm(l,-m0,1,0,l0,m0)    

         allocate(Rfun(Nr,2*l+1,3)); Rfun = zero

         do m= -l, l
            i = m + l + 1
            sig = minusone_n(m)
            if(m == m0-1) then
               do ir=1,Nr
                  rv = rs(ir) * rwf%Eval(rs(ir))
                  Rfun(ir,i,1) = Rfun(ir,i,1) + sig * gnt(1) * rv / sqrt(2.0d0)
                  Rfun(ir,i,2) = Rfun(ir,i,2) + sig * iu * gnt(1) * rv / sqrt(2.0d0)
               end do
            elseif(m == m0 + 1) then
               do ir=1,Nr
                  rv = rs(ir) * rwf%Eval(rs(ir))
                  Rfun(ir,i,1) = Rfun(ir,i,1) - sig * gnt(2) * rv / sqrt(2.0d0)
                  Rfun(ir,i,2) = Rfun(ir,i,2) + sig * iu * gnt(2) * rv / sqrt(2.0d0)
               end do            
            elseif(m == m0) then
               do ir=1,Nr
                  rv = rs(ir) * rwf%Eval(rs(ir))
                  Rfun(ir,i,3) = Rfun(ir,i,3) + sig * gnt(3) * rv
               end do
            end if
         end do
         Rfun = sqrt(QPI/3.0d0) * Rfun 

         call dipwf_m1%Init(rs, Rfun, 2*l+1, 3)
      else
         allocate(Rfun(Nr,1,3)); Rfun = zero
         call dipwf_m1%Init(rs, Rfun, 1, 3)
      end if

      deallocate(Rfun)
      deallocate(rs)

   end subroutine ApplyLengthOp
!-------------------------------------------------------------------------------------- 

!-------------------------------------------------------------------------------------- 
   subroutine ApplyMomentumOp(rwf,l0,m0,dipwf_m1,dipwf_p1)
      integer,parameter :: Nr=400
      real(dp),parameter :: rsmall=1.0e-6_dp
      type(radialwf_t),intent(in)            :: rwf
      integer,intent(in)                     :: l0,m0
      type(cplx_matrix_spline_t),intent(out) :: dipwf_m1,dipwf_p1
      integer :: l,m,i,ir,idir
      real(dp) :: rmax,rv1,rv2
      complex(dp) :: mel_ang(3)
      real(dp),allocatable :: rs(:)
      complex(dp),allocatable :: Rfun(:,:,:)

      Rmax = rwf%Rmax
      rs = linspace(0.0_dp, Rmax, Nr)

      !.......................................
      !            l = l0 + 1
      !.......................................
      l = l0 + 1 

      allocate(Rfun(Nr,2*l+1,3)); Rfun = zero
      do m= -l, l
         i = m + l + 1
         mel_ang = AngularMatrixElement(l,m,l0,m0)
         do ir=1,Nr
            rv1 = rwf%Eval(rs(ir),idx=1)
            if(rs(ir) < rsmall) then
               rv2 = rwf%Eval(rs(ir))/(rs(ir) + rsmall)
            else
               rv2 = rwf%Eval(rs(ir))/rs(ir)
            end if

            do idir=1,3
               Rfun(ir,i,idir) = (rv1 + (0.5_dp* (l0*(l0+1) - l*(l+1)) + 1.0_dp) * rv2) &
                  * mel_ang(idir)
            end do
         end do
      end do

      call dipwf_p1%Init(rs, Rfun, 2*l+1, 3)

      deallocate(Rfun)
      !.......................................
      !            l = l0 - 1
      !.......................................
      l = l0 - 1

      if(l0 > 0) then
         allocate(Rfun(Nr,2*l+1,3)); Rfun = zero
         do m= -l, l
            i = m + l + 1
            mel_ang = AngularMatrixElement(l,m,l0,m0)
            do ir=1,Nr
               rv1 = rwf%Eval(rs(ir),idx=1)
               if(rs(ir) < rsmall) then
                  rv2 = rwf%Eval(rs(ir))/(rs(ir) + rsmall)
               else
                  rv2 = rwf%Eval(rs(ir))/rs(ir)
               end if

               do idir=1,3
                  Rfun(ir,i,idir) = (rv1 + (0.5_dp* (l0*(l0+1) - l*(l+1)) + 1.0_dp) * rv2) &
                      * mel_ang(idir)
               end do
            end do
         end do

         call dipwf_m1%Init(rs, Rfun, 2*l+1, 3)
      else
         allocate(Rfun(Nr,1,3)); Rfun = zero
         call dipwf_m1%Init(rs, Rfun, 1, 3)
      end if

      deallocate(Rfun)
      deallocate(rs)

   end subroutine ApplyMomentumOp
!-------------------------------------------------------------------------------------- 
   subroutine dipole_lambda_projection(lmax,l1,dipwf,lambda,rc,lam_dip_Rfun)
      real(dp),parameter :: epsabs=1.0e-8_dp,lam_small=1.0e-6_dp
      integer,intent(in)                     :: lmax
      integer,intent(in)                     :: l1
      type(cplx_matrix_spline_t),intent(in)  :: dipwf
      real(dp),intent(in)                    :: lambda
      real(dp),intent(in)                    :: rc
      complex(dp),intent(inout)              :: lam_dip_Rfun(:,:)
      integer :: nstep,idir,l,m,icmp
      integer :: inbvx(2)
      complex(dp),allocatable :: Rfun(:,:)

      nstep = max(nint(lambda * rc),1)

      lam_dip_Rfun = zero
      do idir=1,3
         icmp = 0
         do l=0,lmax
            do m=-l,l
               icmp = icmp + 1
               if(abs(m) <= l1) then
                  inbvx = 1
                  lam_dip_Rfun(icmp,idir) = QPI * Lebedev_integral(spherical_func,epsabs,.false.)
               end if
            end do
         end do
      end do

      !.................................................
      contains
      !.................................................
         complex(dp) function spherical_func(theta,phi)
            real(dp),intent(in) :: theta,phi
            integer :: i1,istep
            complex(dp) :: dfun
            complex(dp) :: dipr

            i1 = m + l1 + 1
            dipr = dipwf%Eval_component(rc,i1,idir,inbvx=inbvx)
            dfun = ylm(l1,m,phi,theta) * dipr

            if(lambda > lam_small) then
               do istep=1,nstep
                  dfun = dfun * exp(lambda * rc * cos(phi) / nstep)
               end do
            end if

            spherical_func = conjg(ylm(l,m,phi,theta)) * dfun

         end function spherical_func
      !.................................................      

   end subroutine dipole_lambda_projection
!-------------------------------------------------------------------------------------- 
   function ScattMatrixElement_Lambda(lmax,swf,bessel_integ,kvec) result(Mk)
      integer,intent(in)                    :: lmax
      type(scattwf_t),intent(in)            :: swf
      type(cplx_matrix_spline_t),intent(in) :: bessel_integ
      real(dp),intent(in)                   :: kvec(3)
      complex(dp)                           :: Mk(3)
      integer :: l,m,ncmp,icmp,idir
      real(dp) :: knrm
      complex(dp) :: rintz,exphi

      ncmp = (lmax+1)**2

      Mk = zero
      knrm = norm2(kvec)

      icmp = 0
      do l=0,lmax
         exphi = conjg(swf%Phase(l,knrm))
         do m=-l,l
            icmp = icmp + 1      

            do idir=1,3
               rintz = bessel_integ%Eval_component(knrm,icmp,idir)
               Mk(idir) = Mk(idir) + QPI * exphi * rintz * Ylm_cart(l,m,kvec)
            end do

         end do
      end do      
 
   end function ScattMatrixElement_Lambda
!-------------------------------------------------------------------------------------- 
   subroutine dipole_lambda_BesselTransform(lmax,lam_dip,swf,kmin,kmax,nk,lam_dip_bessel)
      real(dp),parameter :: small=1.0e-10_dp
      integer,intent(in)                     :: lmax
      type(cplx_matrix_spline_t),intent(in)  :: lam_dip
      type(scattwf_t),intent(in)             :: swf
      real(dp),intent(in)                    :: kmin,kmax
      integer,intent(in)                     :: nk
      type(cplx_matrix_spline_t),intent(out) :: lam_dip_bessel
      real(dp),allocatable :: ks(:)
      integer :: ik,nlm,ilm,icmp,idir,dir_flag,l,m
      integer :: inbvx(2)
      real(dp) :: Rmax,yr,yi,knrm
      integer,allocatable :: lvals(:)
      complex(dp),allocatable :: besslk(:,:,:)

      nlm = (lmax+1)**2
      Rmax = lam_dip%xlim(2) - small
      allocate(lvals(nlm))

      ilm = 0
      do l=0,lmax
         do m=-l,l
            ilm = ilm + 1
            lvals(ilm) = l
         end do
      end do

      ks = linspace(kmin,kmax,nk)
      allocate(besslk(nk,nlm,3))

      do ilm=1,nlm
         do idir=1,3
            do ik=1,nk
               icmp = ilm
               dir_flag = idir
               knrm = ks(ik)
               inbvx = 1
               ! call integral_1d(radial_func_re,0.0_dp,Rmax,quad_tol,yr)
               ! call integral_1d(radial_func_im,0.0_dp,Rmax,quad_tol,yi) 
               ! besslk(ik,ilm,idir) = cmplx(yr, yi, kind=dp)
               besslk(ik,ilm,idir) = ZGregoryKintegral(radial_func,knrm,0.0_dp,Rmax)
            end do
         end do
      end do

      call lam_dip_bessel%Init(ks, besslk, nlm, 3)

      deallocate(besslk)
      deallocate(ks)
      deallocate(lvals)

      !.................................................
      contains
      !.................................................      
         complex(dp) function radial_func(r)
            real(dp),intent(in) :: r
            complex(dp) :: y

            y = lam_dip%Eval_component(r,icmp,dir_flag,inbvx=inbvx)
            ! y = lam_dip%Eval_component(r,icmp,dir_flag)
            radial_func = r**2 * y * swf%Eval(lvals(icmp), knrm, r)

         end function radial_func
      !.................................................
         real(dp) function radial_func_re(r)
            real(dp),intent(in) :: r
            complex(dp) :: y

            y = lam_dip%Eval_component(r,icmp,dir_flag,inbvx=inbvx)
            ! y = lam_dip%Eval_component(r,icmp,dir_flag)
            radial_func_re = r**2 * dble(y)* swf%Eval(lvals(ilm), knrm, r)

         end function radial_func_re
      !.................................................      
         real(dp) function radial_func_im(r)
            real(dp),intent(in) :: r
            complex(dp) :: y

            y = lam_dip%Eval_component(r,icmp,dir_flag,inbvx=inbvx)
            ! y = lam_dip%Eval_component(r,icmp,dir_flag)
            radial_func_im = r**2 * aimag(y) * swf%Eval(lvals(ilm), knrm, r)

         end function radial_func_im
      !.................................................    

   end subroutine dipole_lambda_BesselTransform

!-------------------------------------------------------------------------------------- 
#ifdef MPI
   subroutine dipole_lambda_BesselTransform_mpi(lmax,lam_dip,swf,kmin,kmax,nk,lam_dip_bessel)
      use mpi
      use Marray1d_dist,only: dist_array1d_t,GetDisplSize1D 
      real(dp),parameter :: small=1.0e-10_dp
      integer,intent(in)                     :: lmax
      type(cplx_matrix_spline_t),intent(in)  :: lam_dip
      type(scattwf_t),intent(in)             :: swf
      real(dp),intent(in)                    :: kmin,kmax
      integer,intent(in)                     :: nk
      type(cplx_matrix_spline_t),intent(out) :: lam_dip_bessel
      real(dp),allocatable :: ks(:)
      integer :: ik,nlm,ilm,idir,l,m,ix,ix_glob
      integer :: num_indx,num_indx_loc
      integer :: inbvx(2)
      real(dp) :: Rmax,yr,yi,knrm
      integer,allocatable :: lvals(:),indx_tab(:,:)
      complex(dp),allocatable :: besslk_loc(:),besslk_glob(:),besslk(:,:,:)
      ! .. parallelization ..
      integer,parameter  :: master=0,from_master=1,from_worker=2
      integer :: ntasks,taskid,ierr
      integer :: nsize
      integer,allocatable :: size_loc(:),displ(:)
      type(dist_array1d_t) :: indx_dist

      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

      nlm = (lmax+1)**2
      Rmax = lam_dip%xlim(2) - small
      allocate(lvals(nlm))

      ilm = 0
      do l=0,lmax
         do m=-l,l
            ilm = ilm + 1
            lvals(ilm) = l
         end do
      end do

      ks = linspace(kmin,kmax,nk)

      num_indx = nlm * 3 * nk
      allocate(indx_tab(num_indx,3))

      ix = 0
      do ilm=1,nlm
         do idir=1,3
            do ik=1,nk
               ix = ix + 1
               indx_tab(ix,1) = ilm
               indx_tab(ix,2) = idir
               indx_tab(ix,3) = ik
            end do
         end do
      end do      

      call indx_dist%Init(ntasks,taskid,num_indx)
      num_indx_loc = indx_dist%N_loc(taskid)

      allocate(besslk_loc(num_indx_loc))

      do ix=1,num_indx_loc
         ix_glob = indx_dist%Indx_Loc2Glob(taskid,ix)
         ilm = indx_tab(ix_glob,1)
         idir = indx_tab(ix_glob,2)
         ik = indx_tab(ix_glob,3)
         knrm = ks(ik)
         inbvx = 1
         besslk_loc(ix) = ZGregoryKintegral(radial_func,knrm,0.0_dp,Rmax)
      end do

      allocate(displ(0:ntasks-1),size_loc(0:ntasks-1))
      call GetDisplSize1D(indx_dist%N_loc,1,displ,nsize,size_loc)

      allocate(besslk_glob(num_indx))
      call MPI_Allgatherv(besslk_loc,nsize,MPI_DOUBLE_COMPLEX,besslk_glob,size_loc,displ,&
            MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)

      allocate(besslk(nk,nlm,3))
      do ix=1,num_indx
         ilm = indx_tab(ix,1)
         idir = indx_tab(ix,2)
         ik = indx_tab(ix,3)
         besslk(ik,ilm,idir) = besslk_glob(ix)         
      end do

      call lam_dip_bessel%Init(ks, besslk, nlm, 3)

      deallocate(indx_tab)
      deallocate(displ,size_loc)
      deallocate(besslk_loc)
      deallocate(besslk_glob)
      deallocate(besslk)
      deallocate(ks)
      deallocate(lvals)
      call indx_dist%Clean()

      !.................................................
      contains
      !.................................................      
         complex(dp) function radial_func(r)
            real(dp),intent(in) :: r
            complex(dp) :: y

            y = lam_dip%Eval_component(r,ilm,idir,inbvx=inbvx)
            radial_func = r**2 * y * swf%Eval(lvals(ilm), knrm, r)

         end function radial_func
      !.................................................
   end subroutine dipole_lambda_BesselTransform_mpi
#endif
!-------------------------------------------------------------------------------------- 
   function ZGregoryKintegral(f,k,a,b,nsample) result(res)
      use Mintegration,only: GregoryIntegral
      real(dp),intent(in) :: k
      real(dp),intent(in) :: a,b
      integer,intent(in),optional :: nsample
      complex(dp) :: res
      interface 
         function f(x) 
            use Mdef,only: dp
            real(dp),intent(in) :: x
            complex(dp) :: f
         end function f
      end interface
      integer :: nsample_
      integer :: n,i
      real(dp) :: dx,xi
      complex(dp),allocatable :: fx(:)

      nsample_ = 30
      if(present(nsample)) nsample_ = nsample

      dx = DPI/k / nsample_
      n = max(nint((b-a) / dx), 80)
      allocate(fx(0:n))
      do i=0,n
         xi = a + (b-a) * i / dble(n)
         fx(i) = f(xi)
      end do

      dx = (b-a) / dble(n)

      res = GregoryIntegral(n,dx,fx)

      deallocate(fx)

   end function ZGregoryKintegral
!-------------------------------------------------------------------------------------- 


!======================================================================================
end module Mmatrix_elements
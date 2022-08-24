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
   use Mradialintegral,only: radialintegral_t
   implicit none
   include "../units_inc.f90"
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: ScattMatrixElement_Momentum, ScattMatrixElement_Length

   interface ScattMatrixElement_Momentum
      module procedure ScattMatrixElement_Momentum_comp, ScattMatrixElement_Momentum_precomp
   end interface ScattMatrixElement_Momentum

   interface ScattMatrixElement_Length
      module procedure ScattMatrixElement_Length_comp, ScattMatrixElement_Length_precomp
   end interface ScattMatrixElement_Length
!-------------------------------------------------------------------------------------- 
   real(dp),parameter :: quad_tol=1.0e-8_dp
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
   function ScattMatrixElement_Momentum_precomp(swf,radialintegral,l0,m0,kvec) result(Md)
      type(scattwf_t),intent(in)        :: swf
      type(radialintegral_t),intent(in) :: radialintegral
      integer,intent(in)                :: l0,m0
      real(dp),intent(in)               :: kvec(3)
      complex(dp),dimension(3)          :: Md
      integer :: l,m
      real(dp) :: knrm,rint(2)
      complex(dp) :: mel_ang(3),exphi

      Md = zero
      knrm = norm2(kvec)

      call radialintegral%Eval_mom(knrm,rint)

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

   end function ScattMatrixElement_Length_comp
!-------------------------------------------------------------------------------------- 
   function ScattMatrixElement_Length_precomp(swf,radialintegral,l0,m0,kvec) result(Mk)
      type(scattwf_t),intent(in)        :: swf
      type(radialintegral_t),intent(in) :: radialintegral
      integer,intent(in)                :: l0,m0
      real(dp),intent(in)               :: kvec(3)
      complex(dp)                       :: Mk(3)
      integer :: l
      real(dp) :: knrm,rint(2)
      real(dp) :: gnt(3)
      complex(dp) :: exphi

      Mk = zero
      knrm = norm2(kvec)

      call radialintegral%Eval_len(knrm,rint)

      l = l0 - 1
      if(l >= 0) then
         gnt(1) = ThreeYlm(l,-m0+1,1,-1,l0,m0)
         gnt(2) = ThreeYlm(l,-m0-1,1,1,l0,m0)
         gnt(3) = ThreeYlm(l,-m0,1,0,l0,m0)
         exphi = conjg(swf%Phase(l,knrm))

         Mk(1) = Mk(1) + exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
            - gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint(1) / sqrt(2.0d0)
         Mk(2) = Mk(2) + iu*exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
            + gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint(1) / sqrt(2.0d0)
         Mk(3) = Mk(3) + exphi * gnt(3) * rint(1) * conjg(Ylm_cart(l,-m0,kvec))
      end if

      l = l0 + 1
      gnt(1) = ThreeYlm(l,-m0+1,1,-1,l0,m0)
      gnt(2) = ThreeYlm(l,-m0-1,1,1,l0,m0)
      gnt(3) = ThreeYlm(l,-m0,1,0,l0,m0)
      exphi = conjg(swf%Phase(l,knrm))

      Mk(1) = Mk(1) + exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
         - gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint(2) / sqrt(2.0d0)
      Mk(2) = Mk(2) + iu*exphi * (gnt(1)*conjg(Ylm_cart(l,-m0+1,kvec)) &
         + gnt(2)*conjg(Ylm_cart(l,-m0-1,kvec))) * rint(2) / sqrt(2.0d0)
      Mk(3) = Mk(3) + exphi * gnt(3) * rint(2) * conjg(Ylm_cart(l,-m0,kvec))              

      Mk = QPI * sqrt(QPI/3.0d0) * Mk

   end function ScattMatrixElement_Length_precomp
!-------------------------------------------------------------------------------------- 



!======================================================================================
end module Mmatrix_elements
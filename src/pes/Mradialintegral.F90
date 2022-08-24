module Mradialintegral
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp
   use Mbsplines,only: spline1d_t
   use Mquadrature,only: integral_1d
   use Mradialwf,only: radialwf_t
   use Mscattwf,only: scattwf_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: radialintegral_t

   integer,parameter :: gauge_len=0, gauge_mom=1
   real(dp),parameter :: quad_tol=1.0e-8_dp
!--------------------------------------------------------------------------------------   
   type radialintegral_t   
      integer  :: l0
      real(dp) :: kmin,kmax
      type(spline1d_t) :: len_spl_m1,len_spl_p1
      type(spline1d_t) :: mom1_spl_m1,mom1_spl_p1,mom2_spl_m1,mom2_spl_p1
   contains
      procedure,public :: Init
      procedure,public :: Eval_len
      procedure,public :: Eval_mom      
   end type radialintegral_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,l0,kmin,kmax,swf,rwf,nk,gauge)
      use Mutils,only: linspace
      integer,parameter   :: kx=4
      class(radialintegral_t)     :: me
      integer,intent(in)          :: l0
      real(dp),intent(in)         :: kmin,kmax
      type(scattwf_t),intent(in)  :: swf
      type(radialwf_t),intent(in) :: rwf
      integer,intent(in),optional :: nk
      integer,intent(in),optional :: gauge      
      integer :: nk_,gauge_
      integer :: l_indx,ik
      integer :: iflag
      real(dp) :: knrm
      real(dp),allocatable :: ks(:),rint(:),rint1(:),rint2(:)

      nk_ = 100
      if(present(nk)) nk_ = nk

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      me%l0 = l0

      me%kmin = kmin; me%kmax = kmax
      ks = linspace(kmin, kmax, nk_)

      select case(gauge_)
      case(gauge_len)
         allocate(rint(nk_))

         l_indx = l0 - 1
         if(l_indx >= 0) then
            do ik=1,nk_
               knrm = ks(ik)
               call integral_1d(radfunc_len,0.0_dp,rwf%Rmax,quad_tol,rint(ik))  
            end do
         else
            rint = 0.0_dp
         end if

         iflag = 0
         call me%len_spl_m1%Init(ks,rint,kx,iflag)

         l_indx = l0 + 1
         do ik=1,nk_
            knrm = ks(ik)
            call integral_1d(radfunc_len,0.0_dp,rwf%Rmax,quad_tol,rint(ik))  
         end do

         iflag = 0
         call me%len_spl_p1%Init(ks,rint,kx,iflag)

         deallocate(rint)

      case(gauge_mom)
         allocate(rint1(nk_),rint2(nk_))
         l_indx = l0 - 1
         if(l_indx >= 0) then
            do ik=1,nk_
               knrm = ks(ik)
               call integral_1d(radfunc_mom1,0.0_dp,rwf%Rmax,quad_tol,rint1(ik))  
               call integral_1d(radfunc_mom2,0.0_dp,rwf%Rmax,quad_tol,rint2(ik))
            end do
         else
            rint1 = 0.0_dp
            rint2 = 0.0_dp
         end if         

         iflag = 0
         call me%mom1_spl_m1%Init(ks,rint1,kx,iflag)        
         iflag = 0 
         call me%mom2_spl_m1%Init(ks,rint2,kx,iflag)

         l_indx = l0 + 1
         do ik=1,nk_
            knrm = ks(ik)
            call integral_1d(radfunc_mom1,0.0_dp,rwf%Rmax,quad_tol,rint1(ik))  
            call integral_1d(radfunc_mom2,0.0_dp,rwf%Rmax,quad_tol,rint2(ik))
         end do         

         iflag = 0
         call me%mom1_spl_p1%Init(ks,rint1,kx,iflag)        
         iflag = 0 
         call me%mom2_spl_p1%Init(ks,rint2,kx,iflag)

         deallocate(rint1,rint2)         
      end select

      deallocate(ks)
      !.................................................
      contains
      !.................................................
      real(dp) function radfunc_len(r)
         real(dp),intent(in) :: r

         radfunc_len = r**3 * rwf%Eval(r) * swf%Eval(l_indx, knrm, r)

      end function radfunc_len
      !.................................................
      real(dp) function radfunc_mom1(r)
         real(dp),intent(in) :: r

         radfunc_mom1 = r**2 * swf%Eval(l_indx,knrm,r) * rwf%Eval(r,idx=1)

      end function radfunc_mom1
      !.................................................
      real(dp) function radfunc_mom2(r)
         real(dp),intent(in) :: r

         radfunc_mom2 = r * swf%Eval(l_indx,knrm,r) * rwf%Eval(r)

      end function radfunc_mom2
      !.................................................

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine Eval_len(me,k,rint)
      class(radialintegral_t)     :: me
      real(dp),intent(in)         :: k
      real(dp),intent(out)        :: rint(2)
      integer :: iflag,inbvx

      ! print*, k, me%kmin, me%kmax
      if(k < me%kmin .or. k > me%kmax) then
         write(output_unit,fmt700) "Eval_len: out of bounds"
         rint = 0.0_dp
         return
      end if

      inbvx = 1
      rint(1) = me%len_spl_m1%Eval(k,0,iflag,inbvx)
      inbvx = 1
      rint(2) = me%len_spl_p1%Eval(k,0,iflag,inbvx)     

   end subroutine Eval_len
!--------------------------------------------------------------------------------------
   subroutine Eval_mom(me,k,rint)
      class(radialintegral_t)     :: me
      real(dp),intent(in)         :: k
      real(dp),intent(out)        :: rint(2)
      integer :: iflag,inbvx
      integer :: l
      real(dp) :: rint1,rint2

      if(k < me%kmin .or. k > me%kmax) then
         write(output_unit,fmt700) "Eval_mom: out of bounds"
         rint = 0.0_dp
         return
      end if

      l = me%l0 - 1

      if(l >= 0) then
         inbvx = 1
         rint1 = me%mom1_spl_m1%Eval(k,0,iflag,inbvx)
         inbvx = 1
         rint2 = me%mom2_spl_m1%Eval(k,0,iflag,inbvx)

         rint(1) = rint1 + (0.5_dp* (me%l0*(me%l0+1) - l*(l+1)) + 1.0_dp) * rint2 
      else
         rint(1) = 0.0_Dp
      end if 

      l = me%l0 + 1
      inbvx = 1
      rint1 = me%mom1_spl_p1%Eval(k,0,iflag,inbvx)
      inbvx = 1
      rint2 = me%mom2_spl_p1%Eval(k,0,iflag,inbvx)
      rint(2) = rint1 + (0.5_dp* (me%l0*(me%l0+1) - l*(l+1)) + 1.0_dp) * rint2 

   end subroutine Eval_mom
!--------------------------------------------------------------------------------------
 
!======================================================================================
end module Mradialintegral
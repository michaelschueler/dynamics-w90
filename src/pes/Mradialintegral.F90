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
      integer  :: lmin,lmax
      real(dp) :: kmin,kmax
      type(spline1d_t),allocatable,dimension(:) :: len_spl,mom1_spl,mom2_spl
   contains
      procedure,public :: Init
      procedure,public :: Eval_len
      procedure,public :: Eval_mom      
   end type radialintegral_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,lmin,lmax,kmin,kmax,swf,rwf,nk,gauge)
      use Mutils,only: linspace
      integer,parameter   :: kx=4
      class(radialintegral_t)     :: me
      integer,intent(in)          :: lmin,lmax
      real(dp),intent(in)         :: kmin,kmax
      type(scattwf_t),intent(in)  :: swf
      type(radialwf_t),intent(in) :: rwf
      integer,intent(in),optional :: nk
      integer,intent(in),optional :: gauge      
      integer :: nk_,gauge_
      integer :: l,l_indx,ik
      integer :: iflag
      real(dp) :: knrm
      real(dp),allocatable :: ks(:),rint(:),rint1(:),rint2(:)

      nk_ = 100
      if(present(nk)) nk_ = nk

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      me%lmin = lmin
      me%lmax = lmax

      me%kmin = kmin; me%kmax = kmax
      ks = linspace(kmin, kmax, nk_)

      select case(gauge_)
      case(gauge_len)
         allocate(me%len_spl(lmin:lmax))
         allocate(rint(nk_))
         do l=lmin,lmax
            l_indx = l
            do ik=1,nk_
               knrm = ks(ik)
               call integral_1d(radfunc_len,0.0_dp,rwf%Rmax,quad_tol,rint(ik))  
            end do

            iflag = 0
            call me%len_spl(l)%Init(ks,rint,kx,iflag)
         end do
         deallocate(rint)
      case(gauge_mom)
         allocate(me%mom1_spl(lmin:lmax),me%mom2_spl(lmin:lmax))
         allocate(rint1(nk_),rint2(nk_))
         do l=lmin,lmax
            l_indx = l
            do ik=1,nk_
               knrm = ks(ik)
               call integral_1d(radfunc_mom1,0.0_dp,rwf%Rmax,quad_tol,rint1(ik))  
               call integral_1d(radfunc_mom2,0.0_dp,rwf%Rmax,quad_tol,rint2(ik))  
            end do

            iflag = 0
            call me%mom1_spl(l)%Init(ks,rint1,kx,iflag)
            iflag = 0
            call me%mom2_spl(l)%Init(ks,rint2,kx,iflag)
         end do
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
   function Eval_len(me,l,k) result(rint)
      class(radialintegral_t)     :: me
      integer,intent(in)          :: l
      real(dp),intent(in)         :: k
      real(dp)                    :: rint
      integer :: iflag,inbvx

      ! print*, k, me%kmin, me%kmax
      if(k < me%kmin .or. k > me%kmax) then
         write(output_unit,fmt700) "Eval_len: out of bounds"
         rint = 0.0_dp
         return
      end if

      if(l < me%lmin .or. l > me%lmax) then
         write(output_unit,fmt700) "Eval_len: l out of bounds"         
         rint = 0.0_dp
         return
      end if

      inbvx = 1
      rint = me%len_spl(l)%Eval(k,0,iflag,inbvx)

   end function Eval_len
!--------------------------------------------------------------------------------------
   function Eval_mom(me,l,l0,k) result(rint)
      class(radialintegral_t)     :: me
      integer,intent(in)          :: l,l0
      real(dp),intent(in)         :: k
      real(dp)                    :: rint
      integer :: iflag,inbvx
      real(dp) :: rint1,rint2

      if(k < me%kmin .or. k > me%kmax) then
         write(output_unit,fmt700) "Eval_mom: out of bounds"
         rint = 0.0_dp
         return
      end if

      if(l < me%lmin .or. l > me%lmax) then
         write(output_unit,fmt700) "Eval_mom: l out of bounds"         
         rint = 0.0_dp
         return
      end if

      inbvx = 1
      rint1 = me%mom1_spl(l)%Eval(k,0,iflag,inbvx)
      inbvx = 1
      rint2 = me%mom2_spl(l)%Eval(k,0,iflag,inbvx)

      rint = rint1 + (0.5_dp* (l0*(l0+1) - l*(l+1)) + 1.0_dp) * rint2      

   end function Eval_mom
!--------------------------------------------------------------------------------------
 
!======================================================================================
end module Mradialintegral
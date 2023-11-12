module pes_radialintegral
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use scitools_def,only: dp
   use scitools_utils,only: stop_error
   use scitools_bsplines,only: spline1d_t
   use scitools_quadrature,only: integral_1d
   use pes_radialwf,only: radialwf_t
   use pes_scattwf,only: scattwf_t
   use pes_scatt_input,only: scatt_input_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: radialinteg_t, radialinteg_gen_t

   integer,parameter :: gauge_len=0, gauge_mom=1
   real(dp),parameter :: quad_tol=1.0e-8_dp
!--------------------------------------------------------------------------------------   
   type radialinteg_t   
      integer  :: l0
      real(dp) :: kmin,kmax
      type(spline1d_t) :: len_spl_m1,len_spl_p1
      type(spline1d_t) :: mom1_spl_m1,mom1_spl_p1,mom2_spl_m1,mom2_spl_p1
   contains
      procedure,public :: Init => radialinteg_init
      procedure,public :: SetFromInput => radialinteg_SetFromInput
      procedure,public :: Eval_len => radialinteg_Eval_len
      procedure,public :: Eval_mom => radialinteg_Eval_mom
   end type radialinteg_t

   type radialinteg_gen_t   
      integer  :: Lmax
      real(dp) :: kmin,kmax
      type(spline1d_t),allocatable,dimension(:) :: spl
   contains
      procedure,public :: Init => radialinteg_gen_init
      procedure,public :: Eval => radialinteg_gen_eval
   end type radialinteg_gen_t

   type :: integ_wrapper_t
      integer  :: l
      real(dp) :: k
      type(scattwf_t) :: swf
      type(radialwf_t) :: rwf
   end type integ_wrapper_t

   type(integ_wrapper_t) :: iwrap
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine InitIntegWrapper(rwf,swf)
      type(radialwf_t),intent(in) :: rwf
      type(scattwf_t),intent(in)  :: swf

      call iwrap%rwf%Copy(rwf)
      call iwrap%swf%Init(swf%wf_type,swf%Z)

   end subroutine InitIntegWrapper
!--------------------------------------------------------------------------------------
   subroutine CleanIntegWrapper()
      call iwrap%rwf%Clean()
   end subroutine CleanIntegWrapper
!--------------------------------------------------------------------------------------
   real(dp) function radfunc_len(r) result(f)
      real(dp),intent(in) :: r

      f = r**3 * iwrap%rwf%Eval(r) * iwrap%swf%Eval(iwrap%l, iwrap%k, r)

   end function radfunc_len
!--------------------------------------------------------------------------------------
   real(dp) function radfunc_mom1(r) result(f)
      real(dp),intent(in) :: r

      f = r**2 * iwrap%swf%Eval(iwrap%l, iwrap%k, r) * iwrap%rwf%Eval(r,idx=1)

   end function radfunc_mom1
!--------------------------------------------------------------------------------------
   real(dp) function radfunc_mom2(r) result(f)
      real(dp),intent(in) :: r

      f = r * iwrap%swf%Eval(iwrap%l, iwrap%k, r) * iwrap%rwf%Eval(r)

   end function radfunc_mom2
!--------------------------------------------------------------------------------------
   real(dp) function radfunc_gen(r) result(f)
      real(dp),intent(in) :: r

      f = r**2 * iwrap%rwf%Eval(r) * iwrap%swf%Eval(iwrap%l, iwrap%k, r)

   end function radfunc_gen
!--------------------------------------------------------------------------------------
   subroutine radialinteg_Init(me,l0,kmin,kmax,swf,rwf,nk,gauge)
      use scitools_utils,only: linspace
      integer,parameter   :: kx=4
      class(radialinteg_t)     :: me
      integer,intent(in)          :: l0
      real(dp),intent(in)         :: kmin,kmax
      type(scattwf_t),intent(in)  :: swf
      type(radialwf_t),intent(in) :: rwf
      integer,intent(in),optional :: nk
      integer,intent(in),optional :: gauge      
      integer :: nk_,gauge_
      integer :: ik,l_indx
      integer :: iflag
      real(dp),allocatable :: ks(:),rint(:),rint1(:),rint2(:)

      integer :: ir
      real(dp) :: rr

      nk_ = 50
      if(present(nk)) nk_ = nk

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      me%l0 = l0

      me%kmin = kmin; me%kmax = kmax
      ks = linspace(kmin, kmax, nk_)

      call InitIntegWrapper(rwf,swf)

      select case(gauge_)
      case(gauge_len)
         allocate(rint(nk_))

         l_indx = l0 - 1
         iwrap%l = l_indx
         if(l_indx >= 0) then
            do ik=1,nk_
               iwrap%k = ks(ik)
               call integral_1d(radfunc_len,0.0_dp,rwf%Rmax,quad_tol,rint(ik))  
            end do
         else
            rint = 0.0_dp
         end if

         iflag = 0
         call me%len_spl_m1%Init(ks,rint,kx,iflag)

         l_indx = l0 + 1
         iwrap%l = l_indx
         do ik=1,nk_
            iwrap%k = ks(ik)
            call integral_1d(radfunc_len,0.0_dp,rwf%Rmax,quad_tol,rint(ik))  
         end do

         iflag = 0
         call me%len_spl_p1%Init(ks,rint,kx,iflag)

         deallocate(rint)

      case(gauge_mom)
         allocate(rint1(nk_),rint2(nk_))
         l_indx = l0 - 1
         iwrap%l = l_indx
         if(l_indx >= 0) then
            do ik=1,nk_
               iwrap%k = ks(ik)
               call integral_1d(radfunc_mom1,0.0_dp,rwf%Rmax,quad_tol,rint1(ik),meth=12)  
               call integral_1d(radfunc_mom2,0.0_dp,rwf%Rmax,quad_tol,rint2(ik),meth=12)
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
         iwrap%l = l_indx
         do ik=1,nk_
            iwrap%k = ks(ik)
            call integral_1d(radfunc_mom1,0.0_dp,rwf%Rmax,quad_tol,rint1(ik),meth=12)  
            call integral_1d(radfunc_mom2,0.0_dp,rwf%Rmax,quad_tol,rint2(ik),meth=12)
         end do         

         iflag = 0
         call me%mom1_spl_p1%Init(ks,rint1,kx,iflag)        
         iflag = 0 
         call me%mom2_spl_p1%Init(ks,rint2,kx,iflag)

         deallocate(rint1,rint2)         
      end select

      call CleanIntegWrapper()

      deallocate(ks)

   end subroutine radialinteg_Init
!--------------------------------------------------------------------------------------
   subroutine radialinteg_SetFromInput(me,l0,scatt_input,iorb,gauge)
      integer,parameter   :: kx=4
      class(radialinteg_t)           :: me
      integer,intent(in)             :: l0
      type(scatt_input_t),intent(in) :: scatt_input
      integer,intent(in)             :: iorb
      integer,intent(in),optional    :: gauge  
      integer :: gauge_
      integer :: l,iflag
      real(dp),allocatable :: ks(:),rzero(:)

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      me%l0 = l0
      allocate(ks(scatt_input%nE))
      allocate(rzero(scatt_input%nE)); rzero = 0.0_dp
      ks(:) = sqrt(2.0_dp * scatt_input%Ex)
      me%kmin = minval(ks)
      me%kmax = maxval(ks)

      select case(gauge_)
      case(gauge_len)
         l = l0 - 1
         if(l >= 0) then
            iflag = 0
            call me%len_spl_m1%Init(ks,scatt_input%radint(:,l,iorb),kx,iflag)
         else
            iflag = 0
            call me%len_spl_m1%Init(ks,rzero,kx,iflag)
         end if

         l = l0 + 1
         iflag = 0
         call me%len_spl_p1%Init(ks,scatt_input%radint(:,l,iorb),kx,iflag)
      case default
         call stop_error("radial integral input not compatible with chosen gauge")
      end select

   end subroutine radialinteg_SetFromInput
!--------------------------------------------------------------------------------------
   subroutine radialinteg_Eval_len(me,k,rint)
      class(radialinteg_t)     :: me
      real(dp),intent(in)         :: k
      real(dp),intent(out)        :: rint(2)
      integer :: iflag,inbvx

      if(k < me%kmin .or. k > me%kmax) then
         write(output_unit,fmt700) "Eval_len: out of bounds"
         rint = 0.0_dp
         return
      end if

      print*, "[Eval_len]", me%len_spl_m1%nx, me%len_spl_p1%nx

      inbvx = 1
      rint(1) = me%len_spl_m1%Eval(k,0,iflag,inbvx)
      inbvx = 1
      rint(2) = me%len_spl_p1%Eval(k,0,iflag,inbvx)     

   end subroutine radialinteg_Eval_len
!--------------------------------------------------------------------------------------
   subroutine radialinteg_Eval_mom(me,k,rint)
      class(radialinteg_t)     :: me
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

   end subroutine radialinteg_Eval_mom
!--------------------------------------------------------------------------------------
 
!--------------------------------------------------------------------------------------
   subroutine radialinteg_gen_Init(me,Lmax,kmin,kmax,swf,rwf,nk)
      use scitools_utils,only: linspace
      integer,parameter   :: kx=4
      class(radialinteg_gen_t)    :: me
      integer,intent(in)          :: Lmax
      real(dp),intent(in)         :: kmin,kmax
      type(scattwf_t),intent(in)  :: swf
      type(radialwf_t),intent(in) :: rwf
      integer,intent(in),optional :: nk
      integer :: nk_
      integer :: l,n,ik
      integer :: iflag
      real(dp),allocatable :: ks(:),rint(:)

      nk_ = 40
      if(present(nk)) nk_ = nk

      me%Lmax = Lmax

      me%kmin = kmin; me%kmax = kmax
      ks = linspace(kmin, kmax, nk_)

      allocate(me%spl(0:Lmax))
      allocate(rint(nk_))

      call InitIntegWrapper(rwf,swf)

      do l=0,me%Lmax
         iwrap%l = l
         do ik=1,nk_
            iwrap%k = ks(ik)
            call integral_1d(radfunc_gen,0.0_dp,rwf%Rmax,quad_tol,rint(ik))  
         end do         
         iflag = 0
         call me%spl(l)%Init(ks,rint,kx,iflag)         
      end do

      call CleanIntegWrapper()

      deallocate(rint)
      deallocate(ks)

   end subroutine radialinteg_gen_Init
!--------------------------------------------------------------------------------------
   function radialinteg_gen_eval(me,l,k) result(integ)
      class(radialinteg_gen_t)    :: me
      integer,intent(in)          :: l
      real(dp),intent(in)         :: k
      real(dp) :: integ
      integer :: iflag,inbvx

      if(k < me%kmin .or. k > me%kmax) then
         write(output_unit,fmt700) "Eval: out of bounds"
         integ = 0.0_dp
         return
      end if

      if(l < 0 .or. l > me%Lmax) then
         write(output_unit,fmt700) "Eval: l out of bounds"
         integ = 0.0_dp
         return         
      end if

      inbvx = 1
      integ = me%spl(l)%Eval(k,0,iflag,inbvx)

   end function radialinteg_gen_eval
!--------------------------------------------------------------------------------------

!======================================================================================
end module pes_radialintegral
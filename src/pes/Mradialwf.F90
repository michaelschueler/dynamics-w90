module Mradialwf
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp
   use Mbsplines,only: spline1d_t
   use Matomic,only: SlaterWF_radial,SlaterWF_radial_deriv
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: radialwf_t

   integer,parameter :: wf_slater=0, wf_grid=1
!--------------------------------------------------------------------------------------   
   type radialwf_t   
      integer  :: wf_type = wf_slater
      integer  :: neff=1,L=0
      real(dp) :: Zeff=1.0_dp,Rmax=100.0_dp
      type(spline1d_t) :: spl
   contains
      procedure,public :: InitGrid
      procedure,public :: InitSlater
      procedure,public :: Eval
      procedure,public :: Clean
   end type radialwf_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine InitGrid(me,rs,wf)
      integer,parameter   :: kx=4
      class(radialwf_t)   :: me
      real(dp),intent(in) :: rs(:),wf(:)
      integer :: iflag

      me%wf_type = wf_grid
      me%Rmax = maxval(rs)

      iflag = 0
      call me%spl%Init(rs,wf,kx,iflag)

   end subroutine InitGrid
!--------------------------------------------------------------------------------------
   subroutine InitSlater(me,Zeff,neff,L)
      class(radialwf_t)   :: me   
      real(dp),intent(in) :: Zeff
      integer,intent(in)  :: neff,L
      integer :: nr

      me%wf_type = wf_slater
      me%Zeff = Zeff
      me%neff = neff
      me%L = L

      me%Rmax = (me%neff - me%L) * (me%L + 1) * 30.0_dp/me%Zeff

   end subroutine InitSlater
!--------------------------------------------------------------------------------------
   subroutine Clean(me)
      class(radialwf_t)   :: me   

      if(me%wf_type == wf_grid) then
         call me%spl%Clean()
      end if

   end subroutine Clean
!--------------------------------------------------------------------------------------
   function Eval(me,r,idx) result(Rr)
      class(radialwf_t)           :: me  
      real(dp),intent(in)         :: r
      integer,intent(in),optional :: idx
      real(dp)                    :: Rr
      integer :: idx_
      integer :: iflag,inbvx

      idx_ = 0
      if(present(idx)) idx_ = idx

      select case(me%wf_type)
      case(wf_slater)
         select case(idx_)
         case(0)
            Rr = SlaterWF_radial(r,me%Zeff,me%neff,me%L)
         case(1)
            Rr = SlaterWF_radial_deriv(r,me%Zeff,me%neff,me%L)
         case default
            Rr = 0.0_dp
         end select
      case(wf_grid)
          if(r < me%Rmax) then
            inbvx = 1
            Rr = me%spl%Eval(r,idx_,iflag,inbvx)
         else
            Rr = 0.0_dp
         end if
      case default
         Rr = 0.0_dp
      end select

   end function Eval
!--------------------------------------------------------------------------------------



!======================================================================================
end module Mradialwf
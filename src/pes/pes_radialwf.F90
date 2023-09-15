module pes_radialwf
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use scitools_def,only: dp
   use scitools_bsplines,only: spline1d_t
   use pes_atomic,only: SlaterWF_radial,SlaterWF_radial_deriv
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: radialwf_t, radialwf_multi_t

   integer,parameter :: wf_slater=0, wf_grid=1
!--------------------------------------------------------------------------------------   
   type radialwf_t   
      integer  :: wf_type = wf_slater
      integer  :: neff=1,L=0
      real(dp) :: Zeff=1.0_dp,Rmax=100.0_dp
      type(spline1d_t) :: spl
   contains
      procedure,public :: InitGrid => radialwf_InitGrid
      procedure,public :: InitSlater => radialwf_InitSlater
      procedure,public :: Copy => radialwf_Copy
      procedure,public :: Eval => radialwf_Eval
      procedure,public :: Clean => radialwf_Clean
   end type radialwf_t
!--------------------------------------------------------------------------------------   
   type radialwf_multi_t   
      integer  :: l0
      real(dp) :: Rmax
      type(spline1d_t),allocatable,dimension(:) :: xspl,yspl
   contains
      procedure,public :: Init => radialwf_multi_Init
      procedure,public :: Eval => radialwf_multi_Eval
      procedure,public :: Clean => radialwf_multi_Clean
   end type radialwf_multi_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine radialwf_InitGrid(me,rs,wf)
      integer,parameter   :: kx=4
      class(radialwf_t)   :: me
      real(dp),intent(in) :: rs(:),wf(:)
      integer :: iflag

      me%wf_type = wf_grid
      me%Rmax = maxval(rs)

      iflag = 0
      call me%spl%Init(rs,wf,kx,iflag)

   end subroutine radialwf_InitGrid
!--------------------------------------------------------------------------------------
   subroutine radialwf_InitSlater(me,Zeff,neff,L)
      class(radialwf_t)   :: me   
      real(dp),intent(in) :: Zeff
      integer,intent(in)  :: neff,L
      integer :: nr

      me%wf_type = wf_slater
      me%Zeff = Zeff
      me%neff = neff
      me%L = L

      me%Rmax = (me%neff - me%L) * (me%L + 1) * 30.0_dp/me%Zeff

   end subroutine radialwf_InitSlater
!--------------------------------------------------------------------------------------
   subroutine radialwf_Copy(me,rwf)
      class(radialwf_t)   :: me 
      type(radialwf_t),intent(in) :: rwf

      me%wf_type = rwf%wf_type
      me%Rmax = rwf%Rmax

      if(rwf%wf_type == wf_slater) then
         me%Zeff = rwf%Zeff
         me%neff = rwf%neff
         me%L = rwf%L     
      else
         if(.not.allocated(me%spl%tx)) &
            allocate(me%spl%tx(size(rwf%spl%tx))) 

         if(.not.allocated(me%spl%Cff)) &
            allocate(me%spl%Cff( size(rwf%spl%Cff) ) )

         me%spl%nx = rwf%spl%nx
         me%spl%kx = rwf%spl%kx
         me%spl%tx = rwf%spl%tx
         me%spl%Cff = rwf%spl%Cff
      end if

   end subroutine radialwf_Copy
!--------------------------------------------------------------------------------------
   subroutine radialwf_Clean(me)
      class(radialwf_t)   :: me   

      if(me%wf_type == wf_grid) then
         call me%spl%Clean()
      end if

   end subroutine radialwf_Clean
!--------------------------------------------------------------------------------------
   function radialwf_Eval(me,r,idx) result(Rr)
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

   end function radialwf_Eval
!--------------------------------------------------------------------------------------
   subroutine radialwf_multi_Init(me,rs,l0,wf)
      integer,parameter   :: kx=4
      class(radialwf_multi_t)   :: me
      real(dp),intent(in)       :: rs(:)
      integer,intent(in)        :: l0
      complex(dp),intent(in)    :: wf(:,:)
      integer :: m,i
      integer :: iflag

      me%Rmax = maxval(rs)
      me%l0 = l0

      allocate(me%xspl(2*me%l0 + 1),me%yspl(2*me%l0 + 1))
      do m=-l0,l0
         i = m + l0 + 1
         iflag = 0
         call me%xspl(i)%Init(rs,dble(wf(:,i)),kx,iflag)
         call me%yspl(i)%Init(rs,aimag(wf(:,i)),kx,iflag)
      end do

   end subroutine radialwf_multi_Init
!--------------------------------------------------------------------------------------
   subroutine radialwf_multi_Eval(me,r,wf_lm)
      class(radialwf_multi_t)   :: me
      real(dp),intent(in) :: r
      complex(dp),intent(inout) :: wf_lm(:)
      integer :: m,i
      integer :: inbvx,iflag
      real(dp) :: wx,wy

      call assert(size(wf_lm) >= 2*me%l0+1)

      do m=-me%l0, me%l0
         i = m + me%l0 + 1
         inbvx = 1
         wx = me%xspl(i)%Eval(r,0,iflag,inbvx)
         inbvx = 1
         wy = me%yspl(i)%Eval(r,0,iflag,inbvx)
         wf_lm(i) = cmplx(wx,wy,kind=dp)
      end do

   end subroutine radialwf_multi_Eval
!--------------------------------------------------------------------------------------
   subroutine radialwf_multi_Clean(me)
      class(radialwf_multi_t)   :: me
      integer :: n,i

      n = size(me%xspl)
      do i=1,n
         call me%xspl(i)%Clean()
         call me%yspl(i)%Clean()
      end do
      deallocate(me%xspl,me%yspl)

   end subroutine radialwf_multi_Clean
!--------------------------------------------------------------------------------------

!======================================================================================
end module pes_radialwf
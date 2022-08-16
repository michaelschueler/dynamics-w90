module Mscattwf
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp
   use Mspecial,only: spherical_bessel_jn
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: scattwf_t

   integer,parameter :: wf_pw=0, wf_coul=1
!--------------------------------------------------------------------------------------   
   type scattwf_t   
      integer  :: wf_type = wf_pw
      real(dp) :: Z=0.0_dp
   contains
      procedure,public :: Init
      procedure,public :: Eval
      procedure,public :: Phase      
   end type scattwf_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,wf_type,Z)
      class(scattwf_t)   :: me
      integer,intent(in) :: wf_type
      real(dp),intent(in),optional :: Z

      me%wf_type = wf_type

      if(present(Z)) me%Z = Z

   end subroutine InitGrid
!--------------------------------------------------------------------------------------
   function Eval(me,l,k,r) result(Rr)
      class(scattwf_t),intent(in) :: me  
      integer,intent(in)          :: l
      real(dp),intent(in)         :: k     
      real(dp),intent(in)         :: r
      real(dp)                    :: Rr
      integer :: IFAIL
      real(dp) :: eta
      real(dp),dimension(0:l) :: FC,GC,FCP,GCP

      idx_ = 0
      if(present(idx)) idx_ = idx

      select case(me%wf_type)
      case(wf_pw)
         Rr = spherical_bessel_jn(l,k * r)
      case(wf_coul)
         eta = -me%Z / k
         call COUL90(k*r, eta, 0.0_dp, l, FC, GC, FCP, GCP, 0, IFAIL)
         Rr = FC(l)
      case default
         Rr = 0.0_dp
      end select

   end function Eval
!--------------------------------------------------------------------------------------
   complex(dp) function Phase(me,l,k)
      class(scattwf_t),intent(in) :: me  
      integer,intent(in)          :: l
      real(dp),intent(in)         :: k   
      real(dp) :: eta
      complex(dp) :: zeta

      select case(me%wf_type)
      case(wf_coul)
         eta = -me%Z / k
         zeta = lacz_gamma(l + 1.0_dp + iu*eta)
         phase = zeta/abs(zeta)
      case default
         phase = one
      end select

      phase = iu**l * phase

   end function Phase
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
end module Mscattwf
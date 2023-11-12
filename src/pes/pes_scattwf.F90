module pes_scattwf
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use scitools_def,only: dp, one, iu
   use scitools_special,only: spherical_bessel_jn
   use scitools_vector_bsplines,only: real_vector_spline_t
   use pes_scatt_input,only: scatt_input_t
   implicit none
   include "../units_inc.f90"
!--------------------------------------------------------------------------------------
   private
   public :: scattwf_t

   integer,parameter :: wf_pw=0, wf_coul=1, wf_input=2
!--------------------------------------------------------------------------------------   
   type scattwf_t   
      logical  :: phase_from_input=.false.
      integer  :: wf_type = wf_pw
      integer  :: lmax
      real(dp) :: Z=0.0_dp
      type(real_vector_spline_t) :: phase_spl
   contains
      procedure,public :: Init
      procedure,public :: Eval
      procedure,public :: Phase      
   end type scattwf_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,wf_type,Z,scatt_input,scatt_iorb)
      class(scattwf_t)   :: me
      integer,intent(in) :: wf_type
      real(dp),intent(in),optional :: Z
      type(scatt_input_t),intent(in),optional :: scatt_input
      integer,intent(in),optional :: scatt_iorb

      me%wf_type = wf_type

      if(present(Z)) me%Z = Z

      if(me%wf_type == wf_input .and. present(scatt_input) .and. present(scatt_iorb)) then

         call me%phase_spl%Init(scatt_input%Ex, scatt_input%phase(:,0:,scatt_iorb), &
            scatt_input%lmax+1)

         me%lmax = scatt_input%lmax

         me%phase_from_input = .true.
      end if

   end subroutine Init
!--------------------------------------------------------------------------------------
   function Eval(me,l,k,r) result(Rr)
      real(dp),parameter :: small=1.0e-6_dp
      class(scattwf_t),intent(in) :: me  
      integer,intent(in)          :: l
      real(dp),intent(in)         :: k     
      real(dp),intent(in)         :: r
      real(dp)                    :: Rr
      integer :: IFAIL
      real(dp) :: eta,x
      real(dp),dimension(0:l) :: FC,GC,FCP,GCP

      select case(me%wf_type)
      case(wf_pw)
         Rr = spherical_bessel_jn(l,k * r)
      case(wf_coul)
         eta = -me%Z / k
         x = k * r
         if(x < small) then
            Rr = Coulomb_smallx(x, eta, l)
         else
            call COUL90(k*r, eta, 0.0_dp, l, FC, GC, FCP, GCP, 0, IFAIL)
            Rr = FC(l)
         end if

         Rr = Rr / (r + small)

      case default
         Rr = 0.0_dp
      end select

   end function Eval
!--------------------------------------------------------------------------------------
   complex(dp) function Phase(me,l,k)
      class(scattwf_t),intent(in) :: me  
      integer,intent(in)          :: l
      real(dp),intent(in)         :: k   
      real(dp) :: eta,Ek
      real(dp) :: phase_l(me%lmax+1)
      complex(dp) :: zeta

      select case(me%wf_type)
      case(wf_coul)
         eta = -me%Z / k
         zeta = lacz_gamma(l + 1.0_dp + iu*eta)
         phase = zeta/abs(zeta)
      case(wf_input) 
         phase = one
         if(me%phase_from_input) then
            Ek = 0.5_dp * k**2
            print*, "[Phase]", Ek, me%phase_spl%xlim(1), me%phase_spl%xlim(2)
            phase_l(:) = me%phase_spl%Eval(Ek)
            if(l <= me%lmax) phase = phase_l(l-1)
         end if
      case default
         phase = one
      end select

      phase = iu**l * phase

   end function Phase
!--------------------------------------------------------------------------------------
   recursive function lacz_gamma(a) result(g)
      complex(dp),intent(in) :: a
      complex(dp) :: g
      real(dp),parameter :: sq2p=sqrt(0.5D0*Qpi),pi=0.25D0*QPi
      integer,parameter :: cg = 7

      ! these precomputed values are taken by the sample code in Wikipedia,
      ! and the sample itself takes them from the GNU Scientific Library
      double precision, dimension(0:8), parameter :: p = &
      (/ 0.99999999999980993d0, 676.5203681218851d0, -1259.1392167224028d0, &
         771.32342877765313d0, -176.61502916214059d0, 12.50734327868690d05, &
         -0.13857109526572012d0, 9.9843695780195716d-6, 1.5056327351493116d-7 /)

      complex(dp) :: t, w, x
      integer :: i

      x = a

      if ( dble(x) < 0.5D0 ) then
         g = pi/( sin(pi*x) * lacz_gamma(1D0-x) )
      else
         x = x - 1D0
         t = p(0)
         do i=1, cg+1
            t = t + p(i)/(x+real(i,kind=dp))
         end do
         w = x + real(cg,8) + 0.5D0
         g = sq2p* w**(x+0.5D0) * exp(-w) * t
      end if
   end function lacz_gamma
!--------------------------------------------------------------------------------------
   function Coulomb_smallx(x,eta,l) result(Fl)
      real(dp),intent(in) :: x
      real(dp),intent(in) :: eta
      integer,intent(in)  :: l  
      real(dp) :: Fl
      real(dp) :: Ceta

      Ceta = Coulomb_prefac(eta,l)

      Fl = Ceta * x**(l+1)

   end function Coulomb_smallx
!------------------------------------------------------------------------
   pure real(dp) function fac10(n)
      integer, intent(in) :: n
      integer :: i
      real(dp) :: q

      if (n == 0) then
         fac10 = 1.0D0
      else
         fac10 = 1.0D0
         q = 1.0D0
         do i = 1, n
            fac10 = fac10 * q / 10.0D0
            q = q + 1.0D0
         end do
      end if

   end function fac10
!--------------------------------------------------------------------------------------
   function Coulomb_prefac(eta,l) result(ceta)
      real(dp),intent(in) :: eta
      integer,intent(in)  :: l
      real(dp) :: ceta
      real(dp) :: agm,lnfac
      complex(dp) :: zl,gm

      zl = cmplx(dble(l) + 1.0_dp, eta ,kind=dp)
      gm = lacz_gamma(zl)
      agm = abs(gm)

      lnfac = log(agm) + l * log(2.0_dp) - 0.5_dp*Pi*eta &
         -(fac10(2*l+1) + (2*l+1)*log(10.0_dp))

      ceta = exp(lnfac)

   end function Coulomb_prefac
!--------------------------------------------------------------------------------------

!======================================================================================
end module pes_scattwf
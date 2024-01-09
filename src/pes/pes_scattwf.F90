module pes_scattwf
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use scitools_def,only: dp, one, iu
   use scitools_special,only: spherical_bessel_jn
   use scitools_vector_bsplines,only: real_vector_spline_t
   use gamma_function,only: lacz_gamma, cdgamma
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
      real(dp),allocatable :: ks(:)

      me%wf_type = wf_type

      if(present(Z)) me%Z = Z

      if(me%wf_type == wf_input .and. present(scatt_input) .and. present(scatt_iorb)) then

         allocate(ks(scatt_input%nE))
         ks = sqrt(2.0_dp * scatt_input%Ex)
         call me%phase_spl%Init(ks, scatt_input%phase(:,0:,scatt_iorb), &
            scatt_input%lmax+1)

         me%lmax = scatt_input%lmax

         me%phase_from_input = .true.
         deallocate(ks)
      end if

   end subroutine Init
!--------------------------------------------------------------------------------------
   function Eval(me,l,k,r) result(Rr)
      real(dp),parameter :: small=5.0e-3_dp,eps=1.0e-8_dp
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
            ! print*, l, eta, x, Rr
         else
            call COUL90(x, eta, 0.0_dp, l, FC, GC, FCP, GCP, 0, IFAIL)
            Rr = FC(l)  / (x + eps)
         end if

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
      real(dp) :: phi_l(me%lmax+1)
      complex(dp) :: zeta

      select case(me%wf_type)
      case(wf_coul)
         eta = -me%Z / k
         zeta = lacz_gamma(l + 1.0_dp + iu*eta)
         phase = zeta/abs(zeta)
      case(wf_input) 
         phase = one
         if(me%phase_from_input) then
            phi_l(:) = me%phase_spl%Eval(k)
            if(l <= me%lmax) phase = exp(iu*phi_l(l-1))
         end if
      case default
         phase = one
      end select

      phase = iu**l * phase

   end function Phase

!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
   function Coulomb_smallx(x,eta,l) result(Fl)
      real(dp),intent(in) :: x
      real(dp),intent(in) :: eta
      integer,intent(in)  :: l  
      real(dp) :: Fl
      real(dp) :: Ceta

      Ceta = Coulomb_prefac(eta,l)

      Fl = Ceta * x**l

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
      real(dp) :: agm,lnfac,log10_fact
      complex(dp) :: zl,gm

      zl = cmplx(l + 1.0_dp, eta ,kind=dp)
      ! gm = lacz_gamma(zl)
      gm = cdgamma(zl)
      agm = abs(gm)
      log10_fact = fac10(2*l+1) / log10(exp(1.0_dp))

      lnfac = log(agm) + l * log(2.0_dp) - 0.5_dp*Pi*eta - log10_fact

      ceta = exp(lnfac)

      print*, eta, l, agm, log10_fact, lnfac, ceta

   end function Coulomb_prefac
!--------------------------------------------------------------------------------------

!======================================================================================
end module pes_scattwf
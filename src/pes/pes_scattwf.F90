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
            Rr = Coulomb_smallx(x, k, eta, l)
            ! print*, l, eta, x, Rr
         else
            call COUL90(x, eta, 0.0_dp, l, FC, GC, FCP, GCP, 0, IFAIL)
            Rr = FC(l)  / (r + eps)
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
! complex Gamma function in double precision
!
   function cdgamma(x)
      complex(dp) cdgamma, x
      real(dp),parameter :: 
      parameter (
     &    pi = 3.14159265358979324d+00, 
     &    pv = 7.31790632447016203d+00, 
     &    pu = 3.48064577727581257d+00, 
     &    pr = 3.27673720261526849d-02, 
     &    p1 = 1.05400280458730808d+01, 
     &    p2 = 4.73821439163096063d+01, 
     &    p3 = 9.11395751189899762d+01, 
     &    p4 = 6.62756400966213521d+01, 
     &    p5 = 1.32280130755055088d+01, 
     &    p6 = 2.93729529320536228d-01)
      parameter (
     &    q1 = 9.99999999999975753d-01, 
     &    q2 = 2.00000000000603851d+00, 
     &    q3 = 2.99999999944915534d+00, 
     &    q4 = 4.00000003016801681d+00, 
     &    q5 = 4.99999857982434025d+00, 
     &    q6 = 6.00009857740312429d+00)
      xr = real(x, kind=dp)
      xi = aimag(x)
      if (xr .lt. 0) then
          wr = 1.0d0 - xr
          wi = -xi
      else
          wr = xr
          wi = xi
      end if
      ur = wr + q6
      vr = ur * (wr + q5) - wi * wi
      vi = wi * (wr + q5) + ur * wi
      yr = p6 + (p5 * ur + p4 * vr)
      yi = p5 * wi + p4 * vi
      ur = vr * (wr + q4) - vi * wi
      ui = vi * (wr + q4) + vr * wi
      vr = ur * (wr + q3) - ui * wi
      vi = ui * (wr + q3) + ur * wi
      yr = yr + (p3 * ur + p2 * vr)
      yi = yi + (p3 * ui + p2 * vi)
      ur = vr * (wr + q2) - vi * wi
      ui = vi * (wr + q2) + vr * wi
      vr = ur * (wr + q1) - ui * wi
      vi = ui * (wr + q1) + ur * wi
      yr = yr + (p1 * ur + vr)
      yi = yi + (p1 * ui + vi)
      ur = vr * wr - vi * wi
      ui = vi * wr + vr * wi
      t = ur * ur + ui * ui
      vr = (yr * ur + yi * ui) + pr * t
      vi = yi * ur - yr * ui
      yr = wr + pv
      ur = 0.5d0 * log(yr * yr + wi * wi) - 1.0d0
      ui = atan2(wi, yr)
      yr = exp(ur * (wr - 0.5d0) - ui * wi - pu) / t
      yi = ui * (wr - 0.5d0) + ur * wi
      ur = yr * cos(yi)
      ui = yr * sin(yi)
      yr = ur * vr - ui * vi
      yi = ui * vr + ur * vi
      if (xr .lt. 0) then
          wr = pi * xr
          wi = exp(pi * xi)
          vi = 1.0d0 / wi
          ur = (vi + wi) * sin(wr)
          ui = (vi - wi) * cos(wr)
          vr = ur * yr + ui * yi
          vi = ui * yr - ur * yi
          ur = 2.0d0 * pi / (vr * vr + vi * vi)
          yr = ur * vr
          yi = ur * vi
      end if
      cdgamma = cmplx(yr, yi)
   end function cdgamma
!--------------------------------------------------------------------------------------
   function Coulomb_smallx(x,k,eta,l) result(Fl)
      real(dp),intent(in) :: x
      real(dp),intent(in) :: k
      real(dp),intent(in) :: eta
      integer,intent(in)  :: l  
      real(dp) :: Fl
      real(dp) :: Ceta

      Ceta = Coulomb_prefac(eta,l)

      Fl = Ceta * k * x**l

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

      zl = cmplx(dble(l) + 1.0_dp, eta ,kind=dp)
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
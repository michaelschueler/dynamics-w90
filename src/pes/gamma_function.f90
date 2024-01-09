module gamma_function
!======================================================================================
   use Mdebug
   use scitools_def,only: dp
   implicit none
   include "../units_inc.f90"
!--------------------------------------------------------------------------------------
   private
   public :: lacz_gamma, cdgamma
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   complex(dp) function cdgamma(x)
      real(dp),parameter :: pi = 3.14159265358979324d+00
      real(dp),parameter :: pv = 7.31790632447016203d+00
      real(dp),parameter :: pu = 3.48064577727581257d+00
      real(dp),parameter :: pr = 3.27673720261526849d-02
      real(dp),parameter :: p1 = 1.05400280458730808d+01
      real(dp),parameter :: p2 = 4.73821439163096063d+01
      real(dp),parameter :: p3 = 9.11395751189899762d+01
      real(dp),parameter :: p4 = 6.62756400966213521d+01
      real(dp),parameter :: p5 = 1.32280130755055088d+01
      real(dp),parameter :: p6 = 2.93729529320536228d-01
      real(dp),parameter :: q1 = 9.99999999999975753d-01
      real(dp),parameter :: q2 = 2.00000000000603851d+00
      real(dp),parameter :: q3 = 2.99999999944915534d+00
      real(dp),parameter :: q4 = 4.00000003016801681d+00
      real(dp),parameter :: q5 = 4.99999857982434025d+00
      real(dp),parameter :: q6 = 6.00009857740312429d+00
      complex(dp),intent(in) :: x
      real(dp) :: ur,ui,vr,vi,wr,wi,xr,xi,yr,yi,t

     !  complex(dp) cdgamma, x
     !  parameter (
     ! &    pi = 3.14159265358979324d+00, 
     ! &    pv = 7.31790632447016203d+00, 
     ! &    pu = 3.48064577727581257d+00, 
     ! &    pr = 3.27673720261526849d-02, 
     ! &    p1 = 1.05400280458730808d+01, 
     ! &    p2 = 4.73821439163096063d+01, 
     ! &    p3 = 9.11395751189899762d+01, 
     ! &    p4 = 6.62756400966213521d+01, 
     ! &    p5 = 1.32280130755055088d+01, 
     ! &    p6 = 2.93729529320536228d-01)
     !  parameter (
     ! &    q1 = 9.99999999999975753d-01, 
     ! &    q2 = 2.00000000000603851d+00, 
     ! &    q3 = 2.99999999944915534d+00, 
     ! &    q4 = 4.00000003016801681d+00, 
     ! &    q5 = 4.99999857982434025d+00, 
     ! &    q6 = 6.00009857740312429d+00)

      xr = real(x, kind=dp)
      xi = aimag(x)
      if (xr .lt. 0.d0) then
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
      if (xr .lt. 0.d0) then
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
      cdgamma = cmplx(yr, yi, kind=dp)

   end function cdgamma
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
         w = x + real(cg,kind=dp) + 0.5D0
         g = sq2p* w**(x+0.5D0) * exp(-w) * t
      end if
   end function lacz_gamma
!--------------------------------------------------------------------------------------

!======================================================================================
end module gamma_function
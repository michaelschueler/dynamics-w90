module Mwignerd
!======================================================================================
   !! this module contains Wigner-D matrices and spherical harmonics
   use Mdef,only: dp,iu
   implicit none
   include '../units_inc.f90'
!-------------------------------------------------------------------------------------- 
   private 
   public :: djmn,ylm,ylm_cart,Xlm_cart
!-------------------------------------------------------------------------------------- 
contains
!-------------------------------------------------------------------------------------- 
   pure real(dp) function fac10(n)
   !! computes n! / 10^n
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
!------------------------------------------------------------------------
   pure real(dp) function djmn(j,m,n,beta)
   !! Computes the Wigner's (small) d-matrix \(d^j_{m n}(\beta)}\) directly
   !! from the definition.
   !! See [Wikipedia](https://en.wikipedia.org/wiki/Wigner_D-matrix).
   ! .. formal arguments ..
      real(dp), intent(in) :: j,m,n !! angular momentum indices
      real(dp), intent(in) :: beta !! Euler angle beta (radian)
      ! .. local variables ..
      integer :: itmin1,itmin2,itmin,itmax1,itmax2,itmax,ij1,ij2,it,iphase,ia,ib,ic
      real(dp) :: cosb2,sinb2,sqrt_fac,sumt,denom,term

      cosb2=cos(0.5D0*beta)
      sinb2=sin(0.5D0*beta)

      !--------------------------------------------------------------------
      ! determine lower and upper limits for summation index it; these
      ! are derived from the requirement that all factorials n! in the
      ! denominator are restricted to values with n >=0.
      !--------------------------------------------------------------------

      itmin1 = 0
      itmin2 = nint(m-n)
      itmin = max(itmin1,itmin2)
      itmax1 = nint(j+m)
      itmax2 = nint(j-n)
      itmax = min(itmax1,itmax2)
      ij1 = nint(j-m)
      ij2 = nint(j+n)
      sqrt_fac = sqrt( fac10(itmax1) * fac10(ij1) * fac10(ij2) * fac10(itmax2) )

      sumt = 0.0D0
      do it = itmin, itmax
         iphase = (-1)**it
         ia = itmax1 - it
         ib = itmax2 - it
         ic = it + nint(n-m)
         denom = fac10(ia) * fac10(ib) * fac10(it) * fac10(ic)
         term = iphase * cosb2**(ia+ib) * sinb2**(it+ic) / denom
         sumt = sumt + term
      end do
      djmn = sqrt_fac * sumt

   end function djmn
!------------------------------------------------------------------------
   pure complex(dp) function ylm(l, m, thrad, phirad)
   !! Computes the spherical harmonic Y_lm (theta,phi) using the
   !! reduced rotation matrix \(d^l_{m 0} (\theta))\ and using the
   !! external function [[fac10]].
   !! Reference: D.M. Brink and G.R. Satchler, Angular Momentum,
   !!            second edition, Oxford University Press, p.22 and p. 145
      integer, intent(in)  :: l, m !! angular momentum quantum numbers l, m
      real(dp), intent(in) :: thrad, phirad !! angles theta and phi (radian)
      !--------------------------------------------------------------------
      !   local variables
      !--------------------------------------------------------------------
      integer :: itmin1,itmin2,itmin,itmax1,itmax2,itmax,it,iphase,ia,ib,ic
      real(dp) :: sqrt_fac, sumt, denom, term, dlm0, const, cosb2, sinb2
      complex(dp) :: exphi
      !--------------------------------------------------------------------
      !  program starts here
      !  first calculate d^l_{m 0} (theta)
      !--------------------------------------------------------------------
      cosb2 = cos(thrad/2.0D0)
      sinb2 = sin(thrad/2.0D0)
      !--------------------------------------------------------------------
      ! determine lower and upper limits for summation index it; these
      ! are derived from the requirement that all factorials n! in the
      ! denominator are restricted to values with n >=0.
      !--------------------------------------------------------------------
      itmin1 = 0
      itmin2 = m
      itmin = max(itmin1,itmin2)
      itmax1 = l+m
      itmax2 = l
      itmax = min(itmax1,itmax2)
      !  write (6,'(10X,A,2I6)') ' itmin, itmax = ', itmin, itmax
      sqrt_fac = sqrt( fac10(l+m) * fac10(l-m) * fac10(l) * fac10(l) )
      !
      sumt = 0.0D0
      do it = itmin, itmax
         iphase = (-1)**it
         ia = l + m - it
         ib = l - it
         ic = it - m
         denom = fac10(ia) * fac10(ib) * fac10(it) * fac10(ic)
         term = iphase * cosb2**(ia+ib) * sinb2**(it+ic) / denom
         sumt = sumt + term
      end do
      dlm0 = sqrt_fac * sumt
      !--------------------------------------------------------------------
      !  now compute Y_{l m} (theta,phi) from d^l_{m 0} (theta)
      !--------------------------------------------------------------------
      const = sqrt( (2.0D0 *l + 1.0D0) / (4.0D0 * pi) )
      exphi = exp( iu * m * phirad )
      ylm = const * exphi * dlm0

   end function ylm
!------------------------------------------------------------------------
   pure complex(dp) function Ylm_cart(l,m,rv)
   !! Computes the complex spherical  harmonic Y_lm(r) for a vector r
   !! in cartesian coordinates. The vector will be normalized internally.
   !! For l <= 2, the result is tabulated, otherwise [[ylm]] will be called.
      integer,intent(in)  :: l,m !! angular momentum quantum numbers l, m
      real(dp),intent(in) :: rv(3) !! the vector r
      real(dp) :: r,xr,yr,zr,theta,phi

      r = norm2(rv)
      xr = rv(1)/r; yr = rv(2)/r; zr = rv(3)/r

      select case(l)
      case(0)
         Ylm_cart = 1.0_dp/sqrt(QPI)
      case(1)
         select case(m)
         case(-1)
            Ylm_cart = 0.5_dp*sqrt(3.0_dp/DPI)*(xr-iu*yr)
         case(0)
            Ylm_cart = 0.5_dp*sqrt(6.0_dp/DPI) * zr
         case(1)
            Ylm_cart = -0.5_dp*sqrt(3.0_dp/DPI)*(xr+iu*yr)
         end select
      case(2)
         select case(m)
         case(-2)
            Ylm_cart = 0.25_dp*sqrt(15.0_dp/DPI) * (xr - iu*yr)**2
         case(-1)
            Ylm_cart = 0.5_dp*sqrt(15.0_dp/DPI) * (xr - iu*yr)*zr
         case(0)
            Ylm_cart = 0.25_dp*sqrt(10.0_dp/DPI) * (2*zr**2 - xr**2 - yr**2)
         case(1)
            Ylm_cart = -0.5_dp*sqrt(15.0_dp/DPI) * (xr + iu*yr)*zr
         case(2)
            Ylm_cart = 0.25_dp*sqrt(15.0_dp/DPI) * (xr + iu*yr)**2
         end select
      case default
         if(abs(zr) > 0.99999_dp) then
            theta = 0.0_dp
            phi = 0.0_dp
         else
            theta = acos(zr)
            phi = atan2(yr,xr)
         end if
         Ylm_cart = Ylm(l,m,theta,phi)
      end select 

  end function Ylm_cart
!------------------------------------------------------------------------
   pure real(dp) function Xlm_cart(l,m,rv)
   !! Computes the real spherical  harmonic X_lm(r) for a vector r
   !! in cartesian coordinates. The vector will be normalized internally
      real(dp),parameter :: sq2=sqrt(2.0_dp)
      integer,intent(in)  :: l,m !! angular momentum quantum numbers l, m
      real(dp),intent(in) :: rv(3) !! the vector r
      real(dp) :: r,xr,yr,zr,theta,phi
      complex(dp) :: y1, y2


      r = norm2(rv)
      xr = rv(1)/r; yr = rv(2)/r; zr = rv(3)/r

      select case(l)
      case(0)
         Xlm_cart = 1.0_dp/sqrt(QPI)
      case(1)
         select case(m)
         case(-1)
            Xlm_cart = 0.5_dp*sqrt(3.0_dp/DPI) * yr
         case(0)
            Xlm_cart = 0.5_dp*sqrt(6.0_dp/DPI) * zr
         case(1)
            Xlm_cart = 0.5_dp*sqrt(3.0_dp/DPI) * xr
         end select
      case(2)
         select case(m)
         case(-2)
            Xlm_cart = 0.5_dp*sqrt(15.0_dp/DPI) * xr * yr
         case(-1)
            Xlm_cart = 0.5_dp*sqrt(15.0_dp/DPI) * yr * zr
         case(0)
            Xlm_cart = 0.25_dp*sqrt(5.0_dp/DPI) * (2*zr**2 - xr**2 - yr**2)
         case(1)
            Xlm_cart = 0.5_dp*sqrt(15.0_dp/DPI) * xr * xr
         case(2)
            Xlm_cart = 0.25_dp*sqrt(15.0_dp/DPI) * (xr**2 - yr**2)
         end select
      case default
         if(abs(zr) > 0.99999_dp) then
            theta = 0.0_dp
            phi = 0.0_dp
         else
            theta = acos(zr)
            phi = atan2(yr,xr)
         end if
         y1 = Ylm(l,m,theta,phi)
         y2 = Ylm(l,-m,theta,phi)
         if(m==0) then
            Xlm_cart = dble(y1)
         elseif(m < 0) then
            Xlm_cart = dble(iu/sq2*(y1 - (-1)**m * y2))
         else
            Xlm_cart = dble(y2 + (-1)**m * y1) / sq2
         end if
      end select 

  end function Xlm_cart
!------------------------------------------------------------------------

!======================================================================================
end module Mwignerd

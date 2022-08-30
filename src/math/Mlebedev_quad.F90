module Mlebedev_quad
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero
   use Mlebedev_weights,only: GetLebedevQuad
   include "../units_inc.f90"
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: Lebedev_integral, Lebedev_integral_fixed

   interface Lebedev_integral
      module procedure DLebedev_integral, ZLebedev_integral
   end interface Lebedev_integral

   interface Lebedev_integral_fixed
      module procedure DLebedev_integral_fixed, ZLebedev_integral_fixed
   end interface Lebedev_integral_fixed
!--------------------------------------------------------------------------------------
   integer,parameter :: max_order=65, max_np=1454
   integer,parameter,dimension(21) :: lorders=[3,5,7,9,11,13,15,17,19,21,23,25,27,29,&
      31,35,41,47,53,59,65]
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   function DLebedev_integral_fixed(f,order,unit_deg) result(res)
      interface 
         function f(theta,phi) result(ftp)
            use Mdef,only: dp
            real(dp),intent(in) :: theta,phi
            real(dp) :: ftp
         end function f
      end interface
      integer,intent(in)  :: order
      logical,intent(in)  :: unit_deg
      real(dp) :: res
      integer :: np,ip
      real(dp),allocatable :: ths(:),phs(:),ws(:)

      allocate(ths(max_np),phs(max_np),ws(max_np))
      call GetLebedevQuad(order,np,ths,phs,ws,unit_deg=unit_deg)
      res = 0.0_dp
      do ip=1,np
         res = res + ws(ip) * f(ths(ip),phs(ip))
      end do

      deallocate(ths,phs,ws)

   end function DLebedev_integral_fixed
!--------------------------------------------------------------------------------------
   function ZLebedev_integral_fixed(f,order,unit_deg) result(res)
      interface 
         function f(theta,phi) result(ftp)
            use Mdef,only: dp
            real(dp),intent(in) :: theta,phi
            complex(dp) :: ftp
         end function f
      end interface
      integer,intent(in)  :: order
      logical,intent(in)  :: unit_deg
      complex(dp) :: res
      integer :: np,ip
      real(dp),allocatable :: ths(:),phs(:),ws(:)

      allocate(ths(max_np),phs(max_np),ws(max_np))
      call GetLebedevQuad(order,np,ths,phs,ws,unit_deg=unit_deg)
      res = zero
      do ip=1,np
         res = res + ws(ip) * f(ths(ip),phs(ip))
      end do

      deallocate(ths,phs,ws)

   end function ZLebedev_integral_fixed
!--------------------------------------------------------------------------------------
   function DLebedev_integral(f,epsabs,unit_deg) result(res)
      interface 
         function f(theta,phi) result(ftp)
            use Mdef,only: dp
            real(dp),intent(in) :: theta,phi
            real(dp) :: ftp
         end function f
      end interface
      real(dp),intent(in) :: epsabs
      logical,intent(in)  :: unit_deg
      real(dp) :: res
      integer :: i,p,np,ip
      real(dp) :: err,res_old
      real(dp),allocatable :: ths(:),phs(:),ws(:)

      allocate(ths(max_np),phs(max_np),ws(max_np))

      i = 2
      err = 1.0e8_dp

      p = lorders(i)
      call GetLebedevQuad(p,np,ths,phs,ws,unit_deg=unit_deg)
      res_old = 0.0_dp
      do ip=1,np
         res_old = res_old + ws(ip) * f(ths(ip),phs(ip))
      end do

      do while(i < 21 .and. err > epsabs)
         i = i + 1
         p = lorders(i)
         call GetLebedevQuad(p,np,ths,phs,ws,unit_deg=unit_deg)
         res = 0.0_dp
         do ip=1,np
            res = res + ws(ip) * f(ths(ip),phs(ip))
         end do
         err = abs(res - res_old)
         res_old = res
      end do

      deallocate(ths,phs,ws)

   end function DLebedev_integral
!--------------------------------------------------------------------------------------
   function ZLebedev_integral(f,epsabs,unit_deg) result(res)
      interface 
         function f(theta,phi) result(ftp)
            use Mdef,only: dp
            real(dp),intent(in) :: theta,phi
            complex(dp) :: ftp
         end function f
      end interface
      real(dp),intent(in) :: epsabs
      logical,intent(in)  :: unit_deg
      complex(dp) :: res
      integer :: i,p,np,ip
      real(dp) :: err
      complex(dp) :: res_old
      real(dp),allocatable :: ths(:),phs(:),ws(:)

      allocate(ths(max_np),phs(max_np),ws(max_np))

      i = 2
      err = 1.0e8_dp

      p = lorders(i)
      call GetLebedevQuad(p,np,ths,phs,ws,unit_deg=unit_deg)
      res_old = zero
      do ip=1,np
         res_old = res_old + ws(ip) * f(ths(ip),phs(ip))
      end do

      do while(i < 21 .and. err > epsabs)
         i = i + 1
         p = lorders(i)
         call GetLebedevQuad(p,np,ths,phs,ws,unit_deg=unit_deg)
         res = zero
         do ip=1,np
            res = res + ws(ip) * f(ths(ip),phs(ip))
         end do
         err = abs(res - res_old)
         res_old = res
      end do

      deallocate(ths,phs,ws)

   end function ZLebedev_integral
!--------------------------------------------------------------------------------------



!======================================================================================
end module Mlebedev_quad
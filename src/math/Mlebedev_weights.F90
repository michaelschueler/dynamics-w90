module Mlebedev_weights
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp
   include "../units_inc.f90"
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private 
   public :: GetLebedevQuad
!--------------------------------------------------------------------------------------
   real(dp), dimension(6) :: lebedev3_theta
   real(dp), dimension(6) :: lebedev3_phi
   real(dp), dimension(6) :: lebedev3_w
   real(dp), dimension(14) :: lebedev5_theta
   real(dp), dimension(14) :: lebedev5_phi
   real(dp), dimension(14) :: lebedev5_w
   real(dp), dimension(26) :: lebedev7_theta
   real(dp), dimension(26) :: lebedev7_phi
   real(dp), dimension(26) :: lebedev7_w
   real(dp), dimension(38) :: lebedev9_theta
   real(dp), dimension(38) :: lebedev9_phi
   real(dp), dimension(38) :: lebedev9_w
   real(dp), dimension(50) :: lebedev11_theta
   real(dp), dimension(50) :: lebedev11_phi
   real(dp), dimension(50) :: lebedev11_w
   real(dp), dimension(74) :: lebedev13_theta
   real(dp), dimension(74) :: lebedev13_phi
   real(dp), dimension(74) :: lebedev13_w
   real(dp), dimension(86) :: lebedev15_theta
   real(dp), dimension(86) :: lebedev15_phi
   real(dp), dimension(86) :: lebedev15_w
   real(dp), dimension(110) :: lebedev17_theta
   real(dp), dimension(110) :: lebedev17_phi
   real(dp), dimension(110) :: lebedev17_w
   real(dp), dimension(146) :: lebedev19_theta
   real(dp), dimension(146) :: lebedev19_phi
   real(dp), dimension(146) :: lebedev19_w
   real(dp), dimension(170) :: lebedev21_theta
   real(dp), dimension(170) :: lebedev21_phi
   real(dp), dimension(170) :: lebedev21_w
   real(dp), dimension(194) :: lebedev23_theta
   real(dp), dimension(194) :: lebedev23_phi
   real(dp), dimension(194) :: lebedev23_w
   real(dp), dimension(230) :: lebedev25_theta
   real(dp), dimension(230) :: lebedev25_phi
   real(dp), dimension(230) :: lebedev25_w
   real(dp), dimension(266) :: lebedev27_theta
   real(dp), dimension(266) :: lebedev27_phi
   real(dp), dimension(266) :: lebedev27_w
   real(dp), dimension(302) :: lebedev29_theta
   real(dp), dimension(302) :: lebedev29_phi
   real(dp), dimension(302) :: lebedev29_w
   real(dp), dimension(350) :: lebedev31_theta
   real(dp), dimension(350) :: lebedev31_phi
   real(dp), dimension(350) :: lebedev31_w
   real(dp), dimension(434) :: lebedev35_theta
   real(dp), dimension(434) :: lebedev35_phi
   real(dp), dimension(434) :: lebedev35_w
   real(dp), dimension(590) :: lebedev41_theta
   real(dp), dimension(590) :: lebedev41_phi
   real(dp), dimension(590) :: lebedev41_w
   real(dp), dimension(770) :: lebedev47_theta
   real(dp), dimension(770) :: lebedev47_phi
   real(dp), dimension(770) :: lebedev47_w
   real(dp), dimension(974) :: lebedev53_theta
   real(dp), dimension(974) :: lebedev53_phi
   real(dp), dimension(974) :: lebedev53_w
   real(dp), dimension(1202) :: lebedev59_theta
   real(dp), dimension(1202) :: lebedev59_phi
   real(dp), dimension(1202) :: lebedev59_w
   real(dp), dimension(1454) :: lebedev65_theta
   real(dp), dimension(1454) :: lebedev65_phi
   real(dp), dimension(1454) :: lebedev65_w
!--------------------------------------------------------------------------------------
   include "./lebedev/lebedev_3_inc.f90"
   include "./lebedev/lebedev_5_inc.f90"
   include "./lebedev/lebedev_7_inc.f90"
   include "./lebedev/lebedev_9_inc.f90"
   include "./lebedev/lebedev_11_inc.f90"
   include "./lebedev/lebedev_13_inc.f90"
   include "./lebedev/lebedev_15_inc.f90"
   include "./lebedev/lebedev_17_inc.f90"
   include "./lebedev/lebedev_19_inc.f90"
   include "./lebedev/lebedev_21_inc.f90"
   include "./lebedev/lebedev_23_inc.f90"
   include "./lebedev/lebedev_25_inc.f90"
   include "./lebedev/lebedev_27_inc.f90"
   include "./lebedev/lebedev_29_inc.f90"
   include "./lebedev/lebedev_31_inc.f90"
   include "./lebedev/lebedev_35_inc.f90"
   include "./lebedev/lebedev_41_inc.f90"
   include "./lebedev/lebedev_47_inc.f90"
   include "./lebedev/lebedev_53_inc.f90"
   include "./lebedev/lebedev_59_inc.f90"
   include "./lebedev/lebedev_65_inc.f90"
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine GetLebedevQuad(p,np,theta,phi,w,unit_deg)
      integer,intent(in)  :: p
      integer,intent(out) :: np
      real(dp),intent(inout) :: theta(:),phi(:),w(:)
      logical,intent(in),optional :: unit_deg
      logical :: unit_deg_

      unit_deg_ = .true.
      if(present(unit_deg)) unit_deg_ = unit_deg
 
      select case(p)
      case(3)
         np = size(lebedev3_theta)
         theta(1:np) = lebedev3_theta(1:np)
         phi(1:np) = lebedev3_phi(1:np)
         w(1:np) = lebedev3_w(1:np)
      case(5)
         np = size(lebedev5_theta)
         theta(1:np) = lebedev5_theta(1:np)
         phi(1:np) = lebedev5_phi(1:np)
         w(1:np) = lebedev5_w(1:np)
      case(7)
         np = size(lebedev7_theta)
         theta(1:np) = lebedev7_theta(1:np)
         phi(1:np) = lebedev7_phi(1:np)
         w(1:np) = lebedev7_w(1:np)
      case(9)
         np = size(lebedev9_theta)
         theta(1:np) = lebedev9_theta(1:np)
         phi(1:np) = lebedev9_phi(1:np)
         w(1:np) = lebedev9_w(1:np)
      case(11)
         np = size(lebedev11_theta)
         theta(1:np) = lebedev11_theta(1:np)
         phi(1:np) = lebedev11_phi(1:np)
         w(1:np) = lebedev11_w(1:np)
      case(13)
         np = size(lebedev13_theta)
         theta(1:np) = lebedev13_theta(1:np)
         phi(1:np) = lebedev13_phi(1:np)
         w(1:np) = lebedev13_w(1:np)
      case(15)
         np = size(lebedev15_theta)
         theta(1:np) = lebedev15_theta(1:np)
         phi(1:np) = lebedev15_phi(1:np)
         w(1:np) = lebedev15_w(1:np)
      case(17)
         np = size(lebedev17_theta)
         theta(1:np) = lebedev17_theta(1:np)
         phi(1:np) = lebedev17_phi(1:np)
         w(1:np) = lebedev17_w(1:np)
      case(19)
         np = size(lebedev19_theta)
         theta(1:np) = lebedev19_theta(1:np)
         phi(1:np) = lebedev19_phi(1:np)
         w(1:np) = lebedev19_w(1:np)
      case(21)
         np = size(lebedev21_theta)
         theta(1:np) = lebedev21_theta(1:np)
         phi(1:np) = lebedev21_phi(1:np)
         w(1:np) = lebedev21_w(1:np)
      case(23)
         np = size(lebedev23_theta)
         theta(1:np) = lebedev23_theta(1:np)
         phi(1:np) = lebedev23_phi(1:np)
         w(1:np) = lebedev23_w(1:np)
      case(25)
         np = size(lebedev25_theta)
         theta(1:np) = lebedev25_theta(1:np)
         phi(1:np) = lebedev25_phi(1:np)
         w(1:np) = lebedev25_w(1:np)
      case(27)
         np = size(lebedev27_theta)
         theta(1:np) = lebedev27_theta(1:np)
         phi(1:np) = lebedev27_phi(1:np)
         w(1:np) = lebedev27_w(1:np)
      case(29)
         np = size(lebedev29_theta)
         theta(1:np) = lebedev29_theta(1:np)
         phi(1:np) = lebedev29_phi(1:np)
         w(1:np) = lebedev29_w(1:np)
      case(31)
         np = size(lebedev31_theta)
         theta(1:np) = lebedev31_theta(1:np)
         phi(1:np) = lebedev31_phi(1:np)
         w(1:np) = lebedev31_w(1:np)
      case(35)
         np = size(lebedev35_theta)
         theta(1:np) = lebedev35_theta(1:np)
         phi(1:np) = lebedev35_phi(1:np)
         w(1:np) = lebedev35_w(1:np)
      case(41)
         np = size(lebedev41_theta)
         theta(1:np) = lebedev41_theta(1:np)
         phi(1:np) = lebedev41_phi(1:np)
         w(1:np) = lebedev41_w(1:np)
      case(47)
         np = size(lebedev47_theta)
         theta(1:np) = lebedev47_theta(1:np)
         phi(1:np) = lebedev47_phi(1:np)
         w(1:np) = lebedev47_w(1:np)
      case(53)
         np = size(lebedev53_theta)
         theta(1:np) = lebedev53_theta(1:np)
         phi(1:np) = lebedev53_phi(1:np)
         w(1:np) = lebedev53_w(1:np)
      case(59)
         np = size(lebedev59_theta)
         theta(1:np) = lebedev59_theta(1:np)
         phi(1:np) = lebedev59_phi(1:np)
         w(1:np) = lebedev59_w(1:np)
      case(65)
         np = size(lebedev65_theta)
         theta(1:np) = lebedev65_theta(1:np)
         phi(1:np) = lebedev65_phi(1:np)
         w(1:np) = lebedev65_w(1:np)
      end select

      if(.not.unit_deg_) then
         theta(1:np) = pi * theta(1:np) / 180.0d0
         phi(1:np) = pi * phi(1:np) / 180.0d0
      end if

   end subroutine GetLebedevQuad
!--------------------------------------------------------------------------------------

!======================================================================================
end module Mlebedev_weights
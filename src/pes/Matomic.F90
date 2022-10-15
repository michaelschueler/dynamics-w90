module Matomic
!======================================================================================
   use Mdebug
   use scitools_def,only: dp,zero,iu,one
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: SlaterWF_radial, SlaterWF_radial_deriv
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   elemental function SlaterWF_radial(r,Z,n,l,lam_z) result(Rrad)
      character(len=1),parameter ::  angular(0:3)=['s','p','d','f']
      real(dp),intent(in)          :: r
      real(dp),intent(in)          :: Z
      integer,intent(in)           :: n,l
      real(dp),intent(in),optional :: lam_z
      real(dp)                     :: Rrad
      real(dp) :: lam_z_
      real(dp) :: rho,Expr
      character(len=2) :: orb

      lam_z_ = 0.0d0
      if(present(lam_z)) lam_z_ = lam_z

      write(orb,'(I0,A1)') n, angular(l)

      rho = 2d0 * r * Z / dble(n)
      Expr = Exp(-rho/2d0 + lam_z_)

      select case ( orb )
      case('1s')
         Rrad = 2d0 * (Z**1.5d0) * Expr
      case('2s')
         Rrad = 1d0/(2d0*sqrt(2d0)) * ( 2d0 - rho ) * (Z**1.5d0) * Expr
      case('2p')
         Rrad = 1d0/(2d0*sqrt(6d0)) * rho * (Z**1.5d0) * Expr
      case('3s')
         Rrad = 1d0/(9d0*sqrt(3d0)) * ( 6d0 - 6d0*rho + rho**2) * (Z**1.5d0) * Expr
      case('3p')
         Rrad = 1d0/(9d0*sqrt(6d0)) * rho * ( 4d0 - rho ) * (Z**1.5d0) * Expr
      case('3d')
         Rrad = 1d0/(9d0*sqrt(30d0)) * (rho**2) * (Z**1.5d0) * Expr
      case('4s')
         Rrad = 1d0/(96d0) * (24d0 - 36d0*rho + 12d0*(rho**2) - rho**3) * (Z**1.5d0) * Expr
      case('4p')
         Rrad = 1d0/(32d0*sqrt(15d0)) * rho* (20d0 - 10d0*rho + rho**2 ) * (Z**1.5d0) * Expr
      case('4d')
         Rrad = 1d0/(96d0*sqrt(5d0)) * (rho**2) * (6d0 - rho) * (Z**1.5d0) * Expr
      case('4f')
         Rrad = 1d0/(96d0*sqrt(35d0)) * rho**3 * Expr
      case('5s')
         Rrad = 1d0/(300d0*sqrt(5d0)) * (120d0 - 240d0*rho + 120d0*(rho**2) &
            - 20d0*(rho**3)+rho**4) * (Z**1.5d0) * Expr
      case('5p')
         Rrad = 1d0/(150d0*sqrt(30d0)) * rho * (120d0 - 90d0*rho + 18d0*(rho**2) - rho**3) * (Z**1.5d0) * Expr
      case('5d')
         Rrad = 1d0/(150d0*sqrt(70d0)) * (rho**2) * (42d0 - 14d0*rho + rho**2) * (Z**1.5d0) * Expr
      case('6s')
         Rrad = 1d0/(2160d0*sqrt(6d0)) * (720d0 - 1800d0*rho + 1200d0*(rho**2) - 300d0*(rho**3) &
            + 30d0*(rho**4) - rho**5) * (Z**1.5d0) * Expr
      case('6p')
         Rrad = 1d0/(432d0*sqrt(210d0)) * rho * (840d0 - 840d0*rho + 252d0*(rho**2) &
            - 28d0*(rho**3) + (rho**4)) * (Z**1.5d0) * Expr
      case('6d')
         Rrad = 1d0/(864d0*sqrt(105d0)) * (rho**2) * (336d0 - 168d0*rho + 24d0*(rho**2) &
            - (rho**3) ) * (Z**1.5d0) * Expr
      case default
         Rrad = 0d0
      end select

   end function SlaterWF_radial
!--------------------------------------------------------------------------------------
   elemental function SlaterWF_radial_deriv(r,Z,n,l,lam_z) result(Rrad)
      character(len=1),parameter ::  angular(0:3)=['s','p','d','f']
      real(dp),intent(in)          :: r
      real(dp),intent(in)          :: Z
      integer,intent(in)           :: n,l
      real(dp),intent(in),optional :: lam_z
      real(dp)                     :: Rrad
      real(dp) :: lam_z_
      real(dp) :: rho,Prho,dPdrho,Expr
      character(len=2) :: orb

      lam_z_ = 0.0d0
      if(present(lam_z)) lam_z_ = lam_z

      write(orb,'(I0,A1)') n, angular(l)
      rho = 2d0 * r * Z / dble(n)
      Expr = Exp(-rho/2d0 + lam_z_)

      select case ( orb )
      case('1s')
         dPdrho = 0.0d0
         Prho = 2.0d0
      case('2s')
         dPdrho = -1d0/(2d0*sqrt(2d0))
         Prho = 1d0/(2d0*sqrt(2d0)) * ( 2d0 - rho )
      case('2p')
         dPdrho = 1d0/(2d0*sqrt(6d0))
         Prho = 1d0/(2d0*sqrt(6d0)) * rho 
      case('3s')
         dPdrho = 1d0/(9d0*sqrt(3d0)) * (-6d0 + 2d0*rho)
         Prho = 1d0/(9d0*sqrt(3d0)) * ( 6d0 - 6d0*rho + rho**2)
      case('3p')
         dPdrho = 1d0/(9d0*sqrt(6d0)) * ( 4d0 - 2d0*rho)
         Prho = 1d0/(9d0*sqrt(6d0)) * rho * ( 4d0 - rho )
      case('3d')
         dPdrho = 1d0/(9d0*sqrt(30d0)) * 2d0*rho
         Prho = 1d0/(9d0*sqrt(30d0)) * (rho**2)
      case('4s')
         dPdrho = 1d0/(96d0) * (- 36d0 + 24d0*rho - 3d0*rho**2)
         Prho = 1d0/(96d0) * (24d0 - 36d0*rho + 12d0*(rho**2) - rho**3)
      case('4p')
         dPdrho = 1d0/(32d0*sqrt(15d0)) * ((20d0 - 10d0*rho + rho**2 ) &
            + rho * (-10d0 + 2d0*rho) )
         Prho = 1d0/(32d0*sqrt(15d0)) * rho* (20d0 - 10d0*rho + rho**2 )
      case('4d')
         dPdrho = 1d0/(96d0*sqrt(5d0)) * (2d0*rho * (6d0 - rho) - (rho**2) )
         Prho = 1d0/(96d0*sqrt(5d0)) * (rho**2) * (6d0 - rho) 
      case('4f')
         dPdrho = 1d0/(96d0*sqrt(35d0)) * 3d0 * rho**2
         Prho = 1d0/(96d0*sqrt(35d0)) * rho**3
      case('5s')
         dPdrho = 1d0/(300d0*sqrt(5d0)) * (- 240d0 + 240d0*rho &
            - 60d0*(rho**2) + 4d0*rho**3)
         Prho = 1d0/(300d0*sqrt(5d0)) * (120d0 - 240d0*rho + 120d0*(rho**2) &
            - 20d0*(rho**3)+rho**4)
      case('5p')
         dPdrho = 1d0/(150d0*sqrt(30d0)) * (120d0 - 180d0*rho + 54d0*rho**2 - 4d0*rho**3)
         Prho = 1d0/(150d0*sqrt(30d0)) * rho * (120d0 - 90d0*rho + 18d0*(rho**2) - rho**3)
      case('5d')
         dPdrho = 1d0/(150d0*sqrt(70d0)) * (84d0*rho - 42d0*rho**2 + 4d0*rho**3)
         Prho = 1d0/(150d0*sqrt(70d0)) * (rho**2) * (42d0 - 14d0*rho + rho**2)
      case('6s')
         dPdrho = 1d0/(2160d0*sqrt(6d0)) * (- 1800d0 + 2400d0*rho - 900d0*(rho**2) &
            + 120d0*(rho**3) - 5d0*rho**4)
         Prho =  1d0/(2160d0*sqrt(6d0)) * (720d0 - 1800d0*rho + 1200d0*(rho**2) - 300d0*(rho**3) &
            + 30d0*(rho**4) - rho**5)
      case('6p')
         dPdrho = 1d0/(432d0*sqrt(210d0)) * (840d0 - 2*840d0*rho + 3*252d0*(rho**2) &
            - 4*28d0*(rho**3) + 5d0*(rho**4))
         Prho = 1d0/(432d0*sqrt(210d0)) * rho * (840d0 - 840d0*rho + 252d0*(rho**2) &
            - 28d0*(rho**3) + (rho**4))
      case('6d')
         dPdrho = 1d0/(864d0*sqrt(105d0)) * ( 2*336d0*rho - 3*168d0*rho**2 + 4*24d0*(rho**3) &
            - 5d0*(rho**4))
         Prho = 1d0/(864d0*sqrt(105d0)) * (rho**2) * (336d0 - 168d0*rho + 24d0*(rho**2) &
            - (rho**3) )
      case default
         dPdrho = 0d0
         Prho = 0d0
      end select

      Rrad = 2d0 * Z / dble(n) * (Z**1.5d0) * Expr * (dPdrho - 0.5d0*Prho) 

   end function SlaterWF_radial_deriv
!--------------------------------------------------------------------------------------

!======================================================================================
end module Matomic

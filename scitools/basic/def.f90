module scitools_def
!! Provides the basic definition of double precision, some constants and 
!! some useful functions.
!======================================================================================
   use,intrinsic::iso_fortran_env,dp=>real64
   implicit none
   include "../units_inc.f90"
!--------------------------------------------------------------------------------------
   ! in case Fortran 2008 standards are not available use this  
   !  integer,parameter::dp=kind(1.0D0)
   real(dp),parameter    :: expmax=1.0e2_dp,expmin=1.0e-10_dp
   complex(dp),parameter :: iu=(0D0,1D0),one=(1D0,0D0),zero=(0D0,0D0)

   integer,parameter :: large_int=8
   !! use ZGEMM instead of MATMUL for matrix rank > large_int
!--------------------------------------------------------------------------------------
   private
   public :: dp, iu, one, zero, large_int
   public :: nfermi, nnfermi, nbose
   public :: Lorentz, Gauss, GaussOne
   public :: save_exp
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   pure elemental real(dp) function save_exp(x)
   !! evaluates exp(x) avoiding overflow or underflow
      real(dp),intent(in) :: x 

      if(x > 30.0_dp) then
         save_exp = Huge(1.0_dp)
      elseif(x < -30.0_dp) then
         save_exp = 0.0_dp
      else
         save_exp = exp(x)
      end if

   end function save_exp
!--------------------------------------------------------------------------------------
   pure elemental real(dp) function nfermi(b,x)
      !! Fermi distribution function
      real(dp),intent(in) :: b !! inverse temperature
      real(dp),intent(in) :: x !! energy
      real(dp)::y

      y=b*x

      if(abs(y)>expmax) then
         if(y>0.0_dp) then
            nfermi=0.0_dp
         else
            nfermi=1.0_dp
         end if
         return
      end if

      nfermi=1.0_dp/(exp(y)+1.0_dp)

   end function nfermi
!--------------------------------------------------------------------------------------
   pure elemental real(dp) function nnfermi(b,x)
      !! Returns 1 - f, where f is the Fermi distribution
      real(dp),intent(in) :: b !! inverse temperature
      real(dp),intent(in) :: x !! energy

      if(x > 0.0_dp) then
         nnfermi = nfermi(b,-x)
      else
         nnfermi = exp(b*x)*nfermi(b,x)
      end if

   end function nnfermi
!--------------------------------------------------------------------------------------  
   pure recursive function nbose(b,x) result(nb)
      !! stable evaluation of the Bose distribution function
      real(dp),intent(in) :: b !! inverse temperature
      real(dp),intent(in) :: x !! energy
      real(dp) :: nb
      real(dp) :: y

      y=b*x
      if(y<0.0_dp) then
         nb=-1.0_dp-nbose(b,-x)
         return
      end if

      if(abs(y)>expmax) then
         nb=0.0_dp
         return
      end if
      if(y<expmin) then
         nb=1.0_dp/y
         return
      end if

      nb=1.0_dp/(exp(b*x)-1.0_dp)

   end function nbose
!-------------------------------------------------------------------------------------- 
   pure elemental real(dp) function Lorentz(a,x)
      !! Lorentz peak function
      real(dp),intent(in) :: a !! width
      real(dp),intent(in) :: x !! argument

      Lorentz=a/pi/(a**2+x**2)

   end function Lorentz
!-------------------------------------------------------------------------------------- 
   pure elemental real(dp) function Gauss(a,x)
      !! Gauss peak function
      real(dp),intent(in) :: a !! width
      real(dp),intent(in) :: x !! argument

      if(abs(x/a) > 10) then
         Gauss = 0.0D0
         return
      end if

      Gauss = 1D0/sqrt(dpi)*exp(-0.5_dp*(x/a)**2)/a

   end function Gauss
!-------------------------------------------------------------------------------------- 
   pure real(dp) function GaussOne(a,x)
      !! Unnormalized Gauss peak function
      real(dp),intent(in) :: a !! width
      real(dp),intent(in) :: x !! argument

      if(abs(x/a) > 10) then
         GaussOne = 0.0D0
         return
      end if

      GaussOne=exp(-0.5_dp*(x/a)**2)

   end function GaussOne
!--------------------------------------------------------------------------------------  

!======================================================================================
end module scitools_def


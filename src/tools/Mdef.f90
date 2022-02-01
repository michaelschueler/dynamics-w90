module MDef
!======================================================================================
!****m* src/Mdef.f90
!
! NAME
!  Mdef.f90
!
! AUTHOR
!  Michael Schueler
!
! DESCRIPTION
!  Contains some basic functions.
!
! CONTAINS
!
!   Mdef.f90/nfermi
!   Mdef.f90/nbose
!   Mdef.f90/Lorentz
!   Mdef.f90/Gauss
!  
! DEFINITIONS
!  dp  ...  kind(double precision)
!****
!======================================================================================
  use,intrinsic::iso_fortran_env,dp=>real64
  implicit none
!--------------------------------------------------------------------------------------
! in case Fortran 2008 standards are not available use this  
!  integer,parameter::dp=kind(1.0D0)
  real(dp),parameter::expmax=1.0e2_dp,expmin=1.0e-10_dp
  complex(dp),parameter::iu=(0D0,1D0),one=(1D0,0D0),zero=(0D0,0D0)

  private
  public &
       dp,&
       iu,&
       one,&
       zero,&
       nfermi,&
       nnfermi,&
       nbose,&
       Lorentz,&
       Gauss,&
       GaussOne
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
!****f* Mdef.f90/nfermi
!
! NAME 
!  nfermi
!
! AUTHOR
!  Michael Schueler
! 
! DESCRIPTION
!  Fermi distribution function
!
! INPUTS
!  b  -  inverse temperature
!  x  -  energy
!  
! OUTPUT
!  Occupation according to Fermi distribtion.
!
! SOURCE
! 
!--------------------------------------------------------------------------------------
  pure elemental real(dp) function nfermi(b,x)
    real(dp),intent(in) :: b,x
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
    real(dp),intent(in) :: b,x

    if(x > 0.0_dp) then
      nnfermi = nfermi(b,-x)
    else
      nnfermi = exp(b*x)*nfermi(b,x)
    end if

  end function nnfermi
!--------------------------------------------------------------------------------------
!***
!-------------------------------------------------------------------------------------- 
!****f* Mdef.f90/nbose
!
! NAME 
!  nbose
!
! AUTHOR
!  Michael Schueler
! 
! DESCRIPTION
!  Bose distribution function
!
! INPUTS
!  b  -  inverse temperature
!  x  -  energy
!  
! OUTPUT
!  Occupation according to Bose distribtion.
!
! SOURCE
! 
!-------------------------------------------------------------------------------------- 
  pure recursive function nbose(b,x) result(nb)
    real(dp),intent(in) :: b,x
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
!***
!-------------------------------------------------------------------------------------- 
!****f* Mdef.f90/Lorentz
!
! NAME 
!  Lorentz
!
! AUTHOR
!  Michael Schueler
! 
! DESCRIPTION
!  Lorentzian function
!
! INPUTS
!  a  -  width
!  x  -  energy
!  
! OUTPUT
!  Lorentzian, normalized
!
! SOURCE
! 
!-------------------------------------------------------------------------------------- 
  pure elemental real(dp) function Lorentz(a,x)
    real(dp),parameter :: pi=4.0_dp*atan(1.0_dp)
    real(dp),intent(in) :: a,x

    Lorentz=a/pi/(a**2+x**2)

  end function Lorentz
!-------------------------------------------------------------------------------------- 
!***
!-------------------------------------------------------------------------------------- 
!****f* Mdef.f90/Gauss
!
! NAME 
!  Gauss
!
! AUTHOR
!  Michael Schueler
! 
! DESCRIPTION
!  Gaussian function
!
! INPUTS
!  a  -  width
!  x  -  energy
!  
! OUTPUT
!  Gaussian, normalized
!
! SOURCE
! 
!-------------------------------------------------------------------------------------- 
  pure elemental real(dp) function Gauss(a,x)
    real(dp),parameter :: dpi=8.0_dp*atan(1.0_dp)
    real(dp),intent(in) :: a,x

    if(abs(x/a) > 10) then
       Gauss = 0.0D0
       return
    end if
    
    Gauss=1D0/sqrt(dpi*a**2)*exp(-0.5_dp*(x/a)**2)

  end function Gauss
!-------------------------------------------------------------------------------------- 
  pure real(dp) function GaussOne(a,x)
    real(dp),parameter :: dpi=8.0_dp*atan(1.0_dp)
    real(dp),intent(in) :: a,x

    if(abs(x/a) > 10) then
       GaussOne = 0.0D0
       return
    end if
    
    GaussOne=exp(-0.5_dp*(x/a)**2)

  end function GaussOne
!--------------------------------------------------------------------------------------  
!*** 
!-------------------------------------------------------------------------------------- 

!======================================================================================
end module MDef


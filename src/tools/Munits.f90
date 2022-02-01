!****m* src/Munits.f90
!
! NAME
!  Munits
!
! AUTHOR
!  Yaroslav Pavlyukh
!
! DESCRIPTION
!   Contains fundamental constants and unit conversion factors. The data is taken
!   from http://physics.nist.gov/cuu/index.html.
!   
!   Physical constants
!   - Electric constant [ F/m ] =[s^4 A^2 / m^3 kg] 
!   - Speed of light [ m / s ]
!   - Planck constant over 2 pi [J s] = [ m^2 kg /s^2]
!   - Electron charge [C] = [s A]
!   - Fine-structure constant
!
!  Mathematical constants
!   - 2*Pi,4*Pi
!
!  Unit conversion factors
!   - HrSI: energy [a.u.] to [J]
!   - BohrSI: length [a.u] to [m]
!   - E0SI: electric field [a.u.] to [V/m]
!   - t0SI: time [a.u.] to [s]
!   - I0SIa: intensity [a.u.] to [J / s m^2]
!   - HreV:  energy [a.u.] to [eV]
!
!  Unit prefixes
!   - Nano=10^-9
!   - Angstrom=10^-10
!   - Femto=10^-15
!   - Cent=10^-2
!
!****
Module Munits
  !-------------------------------------------------
  ! Data from http://physics.nist.gov/cuu/index.html
  !---------------Fundamental Constants-------------
  ! Electric constant [ F/m ] =[s 4 A 2 / m 3 kg] 
  double precision,parameter::Epsilon0SI=8.854187817D-12
  ! Speed of light [ m / s ]
  double precision,parameter::CSI=2.99792458D8
  ! Planck constant over 2 pi [J s] = [ m 2 kg /s 2]
  double precision,parameter::HbarSI=1.054571628D-34
  ! Electron charge [C] = [s A]
  double precision,parameter::eSI=1.602176487D-19
  ! Electron mass [kg]
  double precision,parameter::mSI=9.10938215D-31
  ! Pi
  double precision,parameter::Pi=3.141592653589793238462643383279D0
  ! 2 Pi
  double precision,parameter::DPi=6.28318530717958647688D0
  ! 4 Pi
  double precision,parameter::QPi=12.56637061435917295376D0
  ! Fine-structure constant
  double precision,parameter::Alpha=eSI**2/(QPi*Epsilon0SI*HbarSi*cSI)
  !---------------Units conversion------------------
  double precision,parameter::HrSI=mSI*(eSI**2/(QPi*Epsilon0SI*HbarSI))**2
  double precision,parameter::BohrSI=(QPi*Epsilon0SI/mSI)*(HbarSI/eSI)**2
  double precision,parameter::BohrAngstrom=0.52917720859D0
  double precision,parameter::E0SI=eSI/(QPi*Epsilon0SI*BohrSI**2)
  double precision,parameter::t0SI=HbarSI/HrSI
  double precision,parameter::I0SIa=0.5D0*Epsilon0SI*cSI*E0SI**2
  double precision,parameter::I0SIb=(HrSI/BohrSI)**2/HbarSI
  double precision,parameter::HreV=HrSI/eSI
  double precision,parameter::Nano=1.0D-9,Angstrom=1.0D-10,Femto=1.0D-15,Cent=1.0D-2
end Module Munits

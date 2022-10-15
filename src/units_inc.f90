  !-------------------------------------------------
  ! Data from http://physics.nist.gov/cuu/index.html
  !---------------Fundamental Constants-------------
  double precision,parameter::Epsilon0SI=8.854187817D-12 !! Electric constant [ F/m ] =[s 4 A 2 / m 3 kg] 
  double precision,parameter::CSI=2.99792458D8 !! Speed of light [ m / s ]
  double precision,parameter::HbarSI=1.054571628D-34 !! Planck constant over 2 pi [J s] = [ m 2 kg /s 2]
  double precision,parameter::eSI=1.602176487D-19 !! Electron charge [C] = [s A]
  double precision,parameter::mSI=9.10938215D-31 !! Electron mass [kg]
  double precision,parameter::Pi=3.141592653589793238462643383279D0 !! Pi
  double precision,parameter::DPi=6.28318530717958647688D0 !! 2 Pi
  double precision,parameter::QPi=12.56637061435917295376D0 !! 4 Pi
  double precision,parameter::Alpha=eSI**2/(QPi*Epsilon0SI*HbarSi*cSI) !! Fine-structure constant
  !---------------Units conversion------------------
  double precision,parameter::HrSI=mSI*(eSI**2/(QPi*Epsilon0SI*HbarSI))**2 !! Units conversion: Hartree to SI
  double precision,parameter::BohrSI=(QPi*Epsilon0SI/mSI)*(HbarSI/eSI)**2 !! Units conversion: Bohr to SI
  double precision,parameter::BohrAngstrom=0.52917720859D0 !! Units conversion: Bohr Angstroms
  double precision,parameter::E0SI=eSI/(QPi*Epsilon0SI*BohrSI**2) !! Unit
  double precision,parameter::t0SI=HbarSI/HrSI
  double precision,parameter::I0SIa=0.5D0*Epsilon0SI*cSI*E0SI**2 !! I
  double precision,parameter::I0SIb=(HrSI/BohrSI)**2/HbarSI
  double precision,parameter::HreV=HrSI/eSI !! Units conversion: Hartree to eV
  double precision,parameter::Nano=1.0D-9,Angstrom=1.0D-10,Femto=1.0D-15,Cent=1.0D-2
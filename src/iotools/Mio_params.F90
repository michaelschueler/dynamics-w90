module Mio_params
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero
   implicit none
   include "../formats.h"  
!--------------------------------------------------------------------------------------
   private
   public :: WannierCalcParams_t, HamiltonianParams_t, PESParams_t
!--------------------------------------------------------------------------------------
   integer,parameter :: field_mode_positions=0,field_mode_dipole=1,field_mode_berry=2
   integer,parameter :: gauge_len=0, gauge_mom=1
   integer,parameter :: wf_pw=0, wf_coul=1
!--------------------------------------------------------------------------------------
   type :: WannierCalcParams_t
   !! Options and parameters for the Wannier calculation performed by [[wann_calc]]
      logical :: calc_orbweight=.false. !! triggers output of orbital weight
      logical :: calc_spin=.false. !! triggers output of spin texture
      logical :: calc_berry=.false. !! triggers output of Berry curvature
      logical :: calc_spin_berry=.false. !! triggers output of spin Berry curvature
      logical :: calc_oam=.false. !! triggers outout of orbital angular momentum (OAM)
      logical :: calc_metric=.false. !! triggers output of the quantum metric
      logical :: calc_evecs=.false. !! triggers output of eigenvectors 
      logical :: berry_valence=.false. !! if .true. we compute Berry-phase properties only for 
                                       !! the occupiend states, while sums over intermediate
                                       !! states are performed over unoccupied states only
      logical :: write_kpts=.false. !! triggers output of the k-points used in the calculation
      integer :: gauge=0 !! gauge for calculating Berry-phase properties: velocity gauge = 0,
                         !! dipole gauge = 1
   contains
      procedure, public :: ReadFromFile => wannier_calc_ReadFromFile
   end type WannierCalcParams_t

   type :: HamiltonianParams_t
   !! Options and parameters for the Wannier Hamiltonian
      character(len=256) :: file_ham="" !! file name for the Wannier Hamiltonian
      character(len=256) :: file_xyz="" !! file name for the Wannier centers
      character(len=256) :: file_lam="" !! file name for SOC constants
      character(len=256) :: file_soc="" !! file name for SOC Hamiltonian
      logical            :: w90_with_soc=.false. !! if .true., we assume the SOC is already 
                                                 !! included in the Hamiltonian from `Wannier90`
      logical            :: apply_field=.false. !! Option to include a static electric field in
                                                !! a non-periodic direction (e.g. out of plane)
      logical            :: slab_mode=.false. !! Option to construct a slab in z direction from
                                              !! a bulk Wannier Hamiltonian
      integer            :: field_mode=field_mode_positions !! How to include the effects of the   
                                                            !! static electric field
      real(dp)           :: MuChem=0.0_dp !! Chemical potential for separating occupied vs. unoccupied bands
      real(dp)           :: energy_thresh=0.0_dp !! Hopping amplitudes smaller than `energy_thresh` will 
                                                 !! be disregarded when compressing the Hamiltonian
      real(dp)           :: Efield(3)=[0.0_dp,0.0_dp,0.0_dp] !! static electric field
      logical            :: use_degen_pert=.false. !! [Expert] triggers the use of degenerate 
                                                   !! perturbation theory when calculating Berry phase
                                                   !! properties
      logical            :: force_herm=.true. !! [Expert] if .true. only the hermitian part of the
                                              !! velocity and dipole matrix is computed
      logical            :: force_antiherm=.true. !! [Export] if .true., only the anti-hermitian
                                                  !! part of the Berry connection is computed
      real(dp)           :: degen_thresh=1.0e-5_dp  !! threshold for considering to bands degenerate
      ! .. slab parameters ..
      integer            :: slab_nlayer=0 !! number of layers in a slab calculation, triggered by
                                          !! `slab_mode = .true.`
      integer            :: slab_max_zhop=10 !! hopping cutoff of hopping in the out-of-plane
                                             !! direction for slab calculation
   contains
      procedure, public :: ReadFromFile => Ham_ReadFromFile  
   end type HamiltonianParams_t

   type :: PESParams_t
      character(len=256) :: file_orbs="" !! file name for Wannier orbitals
      logical            :: kpts_reduced=.true. !! if .true. we assume the input k-points are in 
                                                !! reduced coordinates
      logical            :: lambda_orbital_term=.false. !! triggers the calculation of atomic matrix
                                                        !! elements with complex wave-vector
      integer            :: gauge=gauge_len !! Gauge for dipole operator \(\hat{\Delta}\).
                                            !! 0: dipole gauge \(\hat{\Delta} = \mathbf{r}\), 
                                            !! 1: velocity gauge \(\hat{\Delta} = \mathbf{p}\).
      integer            :: scatt_type=wf_pw !! Type of final states. 0: plane waves, 1: Coulomb waves
      integer            :: Nepe=1 !! number of energy points for outputting spectra
      integer            :: radint_numpoints_k=40 !! Number of points for interpolating radial integrals. 
      integer            :: radint_numpoints_r=256 !! Number of radial points for applying the dipole operator.
                                                   !! Only relevant for `lambda_orbital_term=.true.`
      integer            :: expansion_lmax=8 !! Angular momentum cutoff for expanding the expontential of
                                             !! the complex wave-vector in the atomic matrix elements.
                                             !! Only relevant for `lambda_orbital_term=.true.`
      real(dp)           :: wphot=1.0_dp !! The photon energy.
      real(dp)           :: Eshift=0.0_dp !! All band energies are shifted by 
                                          !! \(\varepsilon_i(\mathbf{k}) \rightarrow \varepsilon_i(\mathbf{k}) + \Delta E\),
                                          !! where \(\Delta E\) is given by `Eshift`. 
      real(dp)           :: Epe_min,Epe_max  !! Interval of final state energies.
      real(dp)           :: lambda_esc=0.0_dp  !! Escape depth. For `lambda_orbital_term=.true.` this also determines
                                               !! the imaginary part of the wave-vector
      real(dp)           :: eta_smear=1.0e-3_dp !! Gaussian smearing of the energy conservation in Fermis Golden rule.
      complex(dp)        :: polvec(3) !! Polarization vector of the photons.
   contains
      procedure, public :: ReadFromFile => PES_ReadFromFile  
   end type PESParams_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine wannier_calc_ReadFromFile(me,fname)
      class(WannierCalcParams_t)   :: me
      character(len=*),intent(in)  :: fname
      integer :: unit_inp
      logical :: calc_orbweight=.false.
      logical :: calc_spin=.false.
      logical :: calc_berry=.false.
      logical :: calc_spin_berry=.false.
      logical :: calc_oam=.false.
      logical :: calc_metric=.false. 
      logical :: calc_evecs=.false.
      logical :: berry_valence=.false.
      logical :: write_kpts=.false.
      integer :: gauge=0
      namelist/CALCOPT/calc_orbweight,calc_spin,calc_berry,calc_spin_berry,&
         calc_oam,calc_metric,calc_evecs,berry_valence,write_kpts,gauge

      open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
      read(unit_inp,nml=CALCOPT)
      close(unit_inp)

      me%calc_orbweight = calc_orbweight
      me%calc_spin = calc_spin
      me%calc_berry = calc_berry
      me%calc_spin_berry = calc_spin_berry
      me%calc_oam = calc_oam
      me%calc_metric = calc_metric
      me%calc_evecs = calc_evecs
      me%berry_valence = berry_valence
      me%write_kpts = write_kpts  


   end subroutine wannier_calc_ReadFromFile
!--------------------------------------------------------------------------------------
   subroutine Ham_ReadFromFile(me,fname)
      class(HamiltonianParams_t)   :: me
      character(len=*),intent(in)  :: fname
      character(len=256) :: file_ham=""
      character(len=256) :: file_xyz=""
      character(len=256) :: file_lam=""
      character(len=256) :: file_soc=""
      logical            :: w90_with_soc=.false.
      logical            :: apply_field=.false.
      logical            :: slab_mode=.false.
      integer            :: field_mode=field_mode_positions
      real(dp)           :: MuChem=0.0_dp
      real(dp)           :: Efield(3)=[0.0_dp,0.0_dp,0.0_dp]
      real(dp)           :: energy_thresh=0.0_dp
      logical            :: use_degen_pert=.false.
      logical            :: force_herm=.true.
      logical            :: force_antiherm=.true.
      real(dp)           :: degen_thresh=1.0e-5_dp  
      namelist/HAMILTONIAN/file_ham,file_xyz,file_lam,file_soc,slab_mode,w90_with_soc,&
         energy_thresh,use_degen_pert,force_herm,force_antiherm,degen_thresh,apply_field,&
         field_mode,Efield,MuChem
      integer :: slab_nlayer=0
      integer :: slab_max_zhop=10
      namelist/SLAB/slab_nlayer,slab_max_zhop

      integer :: unit_inp

      open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
      read(unit_inp,nml=HAMILTONIAN)
      close(unit_inp)

      me%file_ham = file_ham
      me%file_xyz = file_xyz
      me%file_lam = file_lam
      me%file_soc = file_soc
      me%w90_with_soc = w90_with_soc
      me%slab_mode = slab_mode
      me%apply_field = apply_field
      me%field_mode = field_mode
      me%Efield = Efield
      me%MuChem = MuChem
      me%energy_thresh = energy_thresh
      me%use_degen_pert = use_degen_pert
      me%force_herm = force_herm
      me%force_antiherm = force_antiherm
      me%degen_thresh = degen_thresh

      if(me%slab_mode) then
         open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
         read(unit_inp,nml=SLAB)
         close(unit_inp)

         me%slab_nlayer = slab_nlayer
         me%slab_max_zhop = slab_max_zhop
      end if

   end subroutine Ham_ReadFromFile
!--------------------------------------------------------------------------------------
   subroutine PES_ReadFromFile(me,fname)
      class(PESParams_t)  :: me
      character(len=*),intent(in)  :: fname
      character(len=256) :: file_orbs=""
      logical            :: kpts_reduced=.true.
      logical            :: lambda_orbital_term=.false.
      integer            :: gauge=gauge_len
      integer            :: scatt_type=wf_pw
      integer            :: Nepe=1
      integer            :: radint_numpoints_k=40
      integer            :: radint_numpoints_r=256
      integer            :: expansion_lmax=8
      real(dp)           :: wphot=1.0_dp
      real(dp)           :: Eshift=0.0_dp
      real(dp)           :: Epe_min,Epe_max  
      real(dp)           :: lambda_esc=0.0_dp      
      real(dp)           :: eta_smear=1.0e-3_dp
      real(dp)           :: polvec_real(3),polvec_imag(3)
      namelist/PESPARAMS/file_orbs,gauge,Nepe,wphot,Eshift,Epe_min,Epe_max,lambda_esc,&
         eta_smear,polvec_real,polvec_imag,kpts_reduced,scatt_type,radint_numpoints_k,&
         radint_numpoints_r,lambda_orbital_term,expansion_lmax
      integer :: unit_inp

      open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
      read(unit_inp,nml=PESPARAMS)
      close(unit_inp)

      me%file_orbs = file_orbs
      me%kpts_reduced = kpts_reduced
      me%lambda_orbital_term = lambda_orbital_term
      me%gauge = gauge
      me%scatt_type = scatt_type
      me%Nepe = Nepe
      me%radint_numpoints_k = radint_numpoints_k
      me%radint_numpoints_r = radint_numpoints_r
      me%expansion_lmax = expansion_lmax
      me%wphot = wphot
      me%Eshift = Eshift
      me%Epe_min = Epe_min
      me%Epe_max = Epe_max
      me%lambda_esc = lambda_esc
      me%eta_smear = eta_smear
      me%polvec = polvec_real + iu * polvec_imag

   end subroutine PES_ReadFromFile
!--------------------------------------------------------------------------------------


!====================================================================================== 
end module Mio_params
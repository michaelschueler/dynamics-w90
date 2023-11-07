module io_params
!! This module provides high-level tools for reading input parameters for the various programs.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero
   use scitools_utils,only: stop_error
   implicit none
   include "../formats.h"  
!--------------------------------------------------------------------------------------
   private
   public :: WannierCalcParams_t, HamiltonianParams_t, PESParams_t, TimeParams_t
!--------------------------------------------------------------------------------------
   integer,parameter :: field_mode_positions=0,field_mode_dipole=1,field_mode_berry=2
   integer,parameter :: gauge_len=0, gauge_mom=1
   integer,parameter :: wf_pw=0, wf_coul=1
   integer,parameter :: prop_unitary=0, prop_rk4=1, prop_rk5=2, prop_hybrid=3
!--------------------------------------------------------------------------------------
   type :: WannierCalcParams_t
   !! Options and parameters for the Wannier calculation performed by [[wann_calc]]
      logical  :: calc_orbweight=.false. !! triggers output of orbital weight
      logical  :: calc_spin=.false. !! triggers output of spin texture
      logical  :: calc_berry=.false. !! triggers output of Berry curvature
      logical  :: calc_spin_berry=.false. !! triggers output of spin Berry curvature
      logical  :: calc_oam=.false. !! triggers outout of orbital angular momentum (OAM)
      logical  :: calc_metric=.false. !! triggers output of the quantum metric
      logical  :: calc_evecs=.false. !! triggers output of eigenvectors 
      logical  :: write_velocity=.false. !! if `.true.` the velocity matrix elements are written to file
      logical  :: berry_valence=.false. !! if `.true.` we compute Berry-phase properties only for 
                                       !! the occupiend states, while sums over intermediate
                                       !! states are performed over unoccupied states only
      logical  :: calc_dos=.false. !! triggers calculation of density of states
      logical  :: calc_pdos=.false. !! triggers calculation of partial density of states
      logical  :: write_kpts=.false. !! triggers output of the k-points used in the calculation
      integer  :: gauge=0 !! gauge for calculating Berry-phase properties: velocity gauge = 0,
                         !! dipole gauge = 1
      integer  :: Nomega=0 !! number of frequency points for calculation of DOS
      real(dp) :: omega_min=0.0_dp,omega_max=1.0_dp !! interval of frequency points for calculation of DOS
      real(dp) :: DeltaE=0.01_dp !! Gaussian broadening for calculation of DOS
   contains
      procedure, public :: ReadFromFile => wannier_calc_ReadFromFile
   end type WannierCalcParams_t

   type :: TimeParams_t
   !! Options and parameters for the time evolution by [[wann_evol]]
      character(len=256) :: file_field="" !! file name for external electric field
      character(len=256) :: file_dens="" !! file for density matrix
      logical  :: restart_evolution=.false. !! if `.true.`, the density matrix is read from `file_dens` and the
                                            !! time evolution is restarted
      integer  :: propagator=prop_unitary !! method for time evolution: 0 ... unitary, 1 ... RK4, 2 ... RK5
                                          !! for propagator=1,2, the equation with
                                          !! phenomenological damping `T1_relax` and decoherence `T2_relax` will be solved.
      integer  :: Nt=100 !! Number of time steps
      integer  :: output_step=1 !! output will we written to file evert `output_step` time steps
      real(dp) :: Tmax=1.0_dp !! propagation time, defining the time step `dt=Tmax / Nt`
      real(dp) :: T1_relax=1.0e10_dp  !! Diagonal relaxation towards the instantaneous equilibrium density matrix.
      real(dp) :: T2_relax=1.0e10_dp  !! Decay time of off-diagonal elements of the density matrix.
   contains
      procedure, public :: ReadFromFile => Time_ReadFromFile  
   end type TimeParams_t

   type :: HamiltonianParams_t
   !! Options and parameters for the Wannier Hamiltonian
      character(len=256) :: file_ham="" !! file name for the Wannier Hamiltonian
      character(len=256) :: file_xyz="" !! file name for the Wannier centers
      character(len=256) :: file_ovlp="" !! file name for the ovlap matrices
      character(len=256) :: file_lam="" !! file name for SOC constants
      character(len=256) :: file_soc="" !! file name for SOC Hamiltonian
      character(len=256) :: file_elpot="" !! file name for scalar electrostatic potential
      logical            :: w90_with_soc=.false. !! if .true., we assume the SOC is already 
                                                 !! included in the Hamiltonian from `Wannier90`
      logical            :: apply_field=.false. !! Option to include a static electric field in
                                                !! a non-periodic direction (e.g. out of plane)
      logical            :: slab_mode=.false. !! Option to construct a slab in z direction from
                                              !! a bulk Wannier Hamiltonian
      logical            :: exclude_orbitals=.false. !! if `.true.`, selected orbitals can be excluded
                                                     !! from the calculation
      logical            :: exclude_layers=.false. !! if `.true.`, selected layers can be excluded
                                                     !! from the calculation (only for slab mode)
      integer,allocatable,dimension(:) :: orbs_excl !! indices or excluded orbitals
      integer,allocatable,dimension(:) :: lays_excl !! indices or excluded layers
      integer            :: norb_exc=0 !! number of excluded orbitals
      integer            :: nlay_exc=0 !! number of excluded layers
      integer            :: field_mode=field_mode_positions !! How to include the effects of the   
                                                            !! static electric field
      real(dp)           :: Beta=1000.0_dp !! inverse temperature
      real(dp)           :: filling=1.0_dp !! number of electrons per unit cell (per spin in the case without SOC)
      real(dp)           :: MuChem=0.0_dp !! Chemical potential for separating occupied vs. unoccupied bands
      logical            :: FixMuChem=.true. !! if `.true.` the occupations will be computed with respect to
                                             !! the input chemical potential `MuChem`
      real(dp)           :: energy_thresh=0.0_dp !! Hopping amplitudes smaller than `energy_thresh` will 
                                                 !! be disregarded when compressing the Hamiltonian
      real(dp)           :: ovlp_thresh=0.0_dp !! overlap amplitudes smaller than `ovlp_thresh` will 
                                                 !! be disregarded when compressing the overlaps
      real(dp)           :: Efield(3)=[0.0_dp,0.0_dp,0.0_dp] !! static electric field
      logical            :: use_degen_pert=.false. !! [Expert] triggers the use of degenerate 
                                                   !! perturbation theory when calculating Berry phase
                                                   !! properties
      logical            :: force_herm=.true. !! [Expert] if .true. only the hermitian part of the
                                              !! velocity and dipole matrix is computed
      logical            :: force_antiherm=.true. !! [Expert] if .true., only the anti-hermitian
                                                  !! part of the Berry connection is computed
      real(dp)           :: degen_thresh=1.0e-5_dp  !! threshold for considering to bands degenerate
      ! .. light-matter coupling ..
      integer            :: lm_gauge=0 !! gauge of light-matter coupling 
      ! .. slab parameters ..
      integer            :: slab_nlayer=0 !! number of layers in a slab calculation, triggered by
                                          !! `slab_mode = .true.`
   contains
      procedure, public :: ReadFromFile => Ham_ReadFromFile  
   end type HamiltonianParams_t

   type :: PESParams_t
      character(len=256) :: file_orbs="" !! file name for Wannier orbitals
      character(len=256) :: file_scatt="" !! file name for scattering input
      logical            :: kpts_reduced=.true. !! if .true. we assume the input k-points are in 
                                                !! reduced coordinates
      logical            :: dipole_approximation=.true. !! if `.false.`, the finite momentum of the 
                                                        !! photons is taken into account.
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
      real(dp)           :: angle_rot_z=0.0_dp !! angle of rotation of the system around the z axis 
      real(dp)           :: lambda_esc=0.0_dp  !! Escape depth. For `lambda_orbital_term=.true.` this also determines
                                               !! the imaginary part of the wave-vector
      real(dp)           :: eta_smear=1.0e-3_dp !! Gaussian smearing of the energy conservation in Fermis Golden rule.
      real(dp)           :: qmom_phot(3)=[0.0_dp,0.0_dp,0.0_dp] !! photon momentum
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
      logical :: calc_dos=.false.
      logical :: calc_pdos=.false.
      logical :: write_velocity=.false.
      logical :: write_kpts=.false.
      integer :: gauge=0
      integer  :: Nomega=0 
      real(dp) :: omega_min=0.0_dp,omega_max=1.0_dp 
      real(dp) :: DeltaE=0.01_dp
      namelist/CALCOPT/calc_orbweight,calc_spin,calc_berry,calc_spin_berry,&
         calc_oam,calc_metric,calc_evecs,berry_valence,write_kpts,write_velocity,gauge,&
         calc_dos,calc_pdos,Nomega,omega_min,omega_max,DeltaE

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
      me%calc_dos = calc_dos
      me%calc_pdos = calc_pdos
      me%berry_valence = berry_valence
      me%write_kpts = write_kpts  
      me%write_velocity = write_velocity

      me%Nomega = Nomega
      me%omega_min = omega_min
      me%omega_max = omega_max
      me%DeltaE = DeltaE

   end subroutine wannier_calc_ReadFromFile
!--------------------------------------------------------------------------------------
   subroutine Time_ReadFromFile(me,fname)
      class(TimeParams_t)          :: me
      character(len=*),intent(in)  :: fname
      integer :: unit_inp
      character(len=256) :: file_field="" 
      character(len=256) :: file_dens=""
      logical  :: restart_evolution=.false.
      integer  :: propagator=0
      integer  :: Nt=100 
      integer  :: output_step=1 
      real(dp) :: Tmax=1.0_dp 
      real(dp) :: T1_relax=1.0e10_dp  
      real(dp) :: T2_relax=1.0e10_dp  
      namelist/TIMEPARAMS/Nt,Tmax,output_step,propagator,T1_relax,T2_relax,&
         file_field,file_dens,restart_evolution

      open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
      read(unit_inp,nml=TIMEPARAMS)
      close(unit_inp)

      me%Nt = Nt
      me%Tmax = Tmax
      me%output_step = output_step
      me%propagator = propagator
      me%T1_relax = T1_relax
      me%T2_relax = T2_relax
      me%file_field = file_field
      me%file_dens = file_dens
      me%restart_evolution = restart_evolution

   end subroutine Time_ReadFromFile
!--------------------------------------------------------------------------------------
   subroutine Ham_ReadFromFile(me,fname)
      class(HamiltonianParams_t)   :: me
      character(len=*),intent(in)  :: fname
      character(len=256) :: file_ham=""
      character(len=256) :: file_xyz=""
      character(len=256) :: file_ovlp=""
      character(len=256) :: file_lam=""
      character(len=256) :: file_soc=""
      character(len=256) :: file_elpot=""
      logical            :: w90_with_soc=.false.
      logical            :: apply_field=.false.
      logical            :: slab_mode=.false.
      integer            :: norb_exc=0 
      integer            :: field_mode=field_mode_positions
      real(dp)           :: Beta=1000.0_dp 
      real(dp)           :: filling=1.0_dp 
      real(dp)           :: MuChem=0.0_dp 
      logical            :: FixMuChem=.true. 
      real(dp)           :: Efield(3)=[0.0_dp,0.0_dp,0.0_dp]
      real(dp)           :: energy_thresh=0.0_dp
      logical            :: use_degen_pert=.false.
      logical            :: force_herm=.true.
      logical            :: force_antiherm=.true.
      real(dp)           :: degen_thresh=1.0e-5_dp  
      integer            :: lm_gauge=0
      character(len=1000) :: exclude_orbitals="",exclude_layers=""
      namelist/HAMILTONIAN/file_ham,file_xyz,file_ovlp,file_lam,file_soc,file_elpot,slab_mode,w90_with_soc,&
         energy_thresh,use_degen_pert,force_herm,force_antiherm,degen_thresh,apply_field,&
         field_mode,Efield,Beta,Filling,MuChem,FixMuChem,lm_gauge,exclude_orbitals,exclude_layers
      integer :: slab_nlayer=0
      namelist/SLAB/slab_nlayer

      integer :: unit_inp

      open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
      read(unit_inp,nml=HAMILTONIAN)
      close(unit_inp)

      me%file_ham = file_ham
      me%file_xyz = file_xyz
      me%file_ovlp = file_ovlp
      me%file_lam = file_lam
      me%file_soc = file_soc
      me%file_elpot = file_elpot
      me%w90_with_soc = w90_with_soc
      me%slab_mode = slab_mode
      me%apply_field = apply_field
      me%field_mode = field_mode
      me%Efield = Efield
      me%Beta = Beta
      me%MuChem = MuChem
      me%Filling = Filling
      me%FixMuChem = FixMuChem
      me%energy_thresh = energy_thresh
      me%use_degen_pert = use_degen_pert
      me%force_herm = force_herm
      me%force_antiherm = force_antiherm
      me%degen_thresh = degen_thresh
      me%lm_gauge = lm_gauge

      if(len_trim(exclude_orbitals) > 0) then
         call ReadList(exclude_orbitals, me%orbs_excl, "exclude_orbitals")
         me%exclude_orbitals = all(me%orbs_excl .ne. 0)
         me%norb_exc = size(me%orbs_excl)
      end if

      if(me%slab_mode) then
         open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
         read(unit_inp,nml=SLAB)
         close(unit_inp)

         me%slab_nlayer = slab_nlayer

         if(len_trim(exclude_layers) > 0) then
            call ReadList(exclude_layers, me%lays_excl, "exclude_layers")
            me%exclude_layers = all(me%lays_excl .ne. 0)
            me%nlay_exc = size(me%lays_excl)    
         end if        
      end if

   end subroutine Ham_ReadFromFile
!--------------------------------------------------------------------------------------
   subroutine PES_ReadFromFile(me,fname)
      class(PESParams_t)  :: me
      character(len=*),intent(in)  :: fname
      character(len=256) :: file_orbs=""
      character(len=256) :: file_scatt=""
      logical            :: kpts_reduced=.true.
      logical            :: dipole_approximation=.true.
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
      real(dp)           :: angle_rot_z=0.0_dp
      real(dp)           :: lambda_esc=0.0_dp      
      real(dp)           :: eta_smear=1.0e-3_dp
      real(dp)           :: polvec_real(3),polvec_imag(3)
      real(dp)           :: qmom_phot(3)=[0.0_dp,0.0_dp,0.0_dp] 
      namelist/PESPARAMS/file_orbs,gauge,Nepe,wphot,Eshift,Epe_min,Epe_max,lambda_esc,&
         eta_smear,polvec_real,polvec_imag,kpts_reduced,scatt_type,radint_numpoints_k,&
         radint_numpoints_r,lambda_orbital_term,expansion_lmax,dipole_approximation,&
         qmom_phot,angle_rot_z,file_scatt
      integer :: unit_inp

      open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
      read(unit_inp,nml=PESPARAMS)
      close(unit_inp)

      me%file_orbs = file_orbs
      me%file_scatt = file_scatt
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
      me%angle_rot_z = angle_rot_z
      me%lambda_esc = lambda_esc
      me%eta_smear = eta_smear
      me%polvec = polvec_real + iu * polvec_imag

      me%dipole_approximation = dipole_approximation
      me%qmom_phot = qmom_phot

   end subroutine PES_ReadFromFile
!--------------------------------------------------------------------------------------
   subroutine ReadList(str,list,func_tag)
      character(len=*),intent(in) :: str
      integer,allocatable,intent(out) :: list(:)
      character(len=*),intent(in) ::    func_tag
      integer :: nl
      integer :: iost
      integer :: ncomma,ndash,i,istart,iend,ppos

      ncomma = count(transfer(str, 'a', len(str)) == ",")
      ndash = count(transfer(str, 'a', len(str)) == "-")

      if(ncomma == 0 .and. ndash == 0) then
         nl = 1
         allocate(list(nl)); list = 0
         read(str, *, iostat=iost) list
         if(iost .ne. 0) then
            write(output_unit,fmt700) trim(func_tag)//": invalid input"
            list = 0
         end if
      elseif(ncomma > 0 .and. ndash == 0) then
         nl = ncomma + 1
         allocate(list(nl)); list = 0
         read(str, *, iostat=iost) list
         if(iost .ne. 0) then
            write(output_unit,fmt700) trim(func_tag)//": invalid input"
            list = 0
         end if
      elseif(ncomma ==0 .and. ndash == 1) then
         ppos = scan(trim(str),"-", BACK= .false.)
         read(str(1:ppos-1), *, iostat=iost) istart
         if(iost .ne. 0) call stop_error(trim(func_tag)//": invalid input")

         read(str(ppos+1:), *, iostat=iost) iend
         if(iost .ne. 0) call stop_error(trim(func_tag)//": invalid input")

         nl = iend - istart + 1
         allocate(list(nl))
         do i=1,nl
            list(i) = i + istart - 1
         end do
      else
         call stop_error(trim(func_tag)//": invalid input")
      end if

   end subroutine ReadList
!--------------------------------------------------------------------------------------
!====================================================================================== 
end module io_params
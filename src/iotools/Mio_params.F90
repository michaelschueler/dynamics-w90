module Mio_params
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero
   implicit none
   include "../formats.h"  
!--------------------------------------------------------------------------------------
   private
   public :: WannierCalcParams_t,HamiltonianParams_t
!--------------------------------------------------------------------------------------
   integer,parameter :: field_mode_positions=0,field_mode_dipole=1,field_mode_berry=2
!--------------------------------------------------------------------------------------
   type :: WannierCalcParams_t
      ! .. calculation options ..
      logical :: slab_mode=.false.
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
      ! .. slab parameters ..
      integer :: slab_nlayer=0
      integer :: slab_max_zhop=10
   contains
      procedure, public :: ReadFromFile => wannier_calc_ReadFromFile
   end type WannierCalcParams_t

   type :: HamiltonianParams_t
      character(len=256) :: file_ham=""
      character(len=256) :: file_xyz=""
      logical            :: w90_with_soc=.false.
      logical            :: apply_field=.false.
      integer            :: field_mode=field_mode_positions
      real(dp)           :: energy_thresh=0.0_dp
      real(dp)           :: Efield(3)=[0.0_dp,0.0_dp,0.0_dp]
      logical            :: use_degen_pert=.false.
      logical            :: force_herm=.true.
      logical            :: force_antiherm=.true.
      real(dp)           :: degen_thresh=1.0e-5_dp  
   contains
      procedure, public :: ReadFromFile => Ham_ReadFromFile  
   end type HamiltonianParams_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine wannier_calc_ReadFromFile(me,fname)
      class(WannierCalcParams_t)   :: me
      character(len=*),intent(in)  :: fname
      integer :: unit_inp
      logical :: slab_mode=.false.
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
      namelist/CALCOPT/slab_mode,calc_orbweight,calc_spin,calc_berry,calc_spin_berry,&
         calc_oam,calc_metric,calc_evecs,berry_valence,write_kpts,gauge
      integer :: slab_nlayer=0
      integer :: slab_max_zhop=10
      namelist/SLAB/slab_nlayer,slab_max_zhop

      open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
      read(unit_inp,nml=CALCOPT)
      close(unit_inp)

      me%slab_mode = slab_mode
      me%calc_orbweight = calc_orbweight
      me%calc_spin = calc_spin
      me%calc_berry = calc_berry
      me%calc_spin_berry = calc_spin_berry
      me%calc_oam = calc_oam
      me%calc_metric = calc_metric
      me%calc_evecs = calc_evecs
      me%berry_valence = berry_valence
      me%write_kpts = write_kpts  

      if(me%slab_mode) then
         open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
         read(unit_inp,nml=SLAB)
         close(unit_inp)

         me%slab_nlayer = slab_nlayer
         me%slab_max_zhop = slab_max_zhop
      end if


   end subroutine wannier_calc_ReadFromFile
!--------------------------------------------------------------------------------------
   subroutine Ham_ReadFromFile(me,fname)
      class(HamiltonianParams_t)   :: me
      character(len=*),intent(in)  :: fname
      character(len=256) :: file_ham=""
      character(len=256) :: file_xyz=""
      logical            :: w90_with_soc=.false.
      logical            :: apply_field=.false.
      integer            :: field_mode=field_mode_positions
      real(dp)           :: Efield(3)=[0.0_dp,0.0_dp,0.0_dp]
      real(dp)           :: energy_thresh=0.0_dp
      logical            :: use_degen_pert=.false.
      logical            :: force_herm=.true.
      logical            :: force_antiherm=.true.
      real(dp)           :: degen_thresh=1.0e-5_dp  
      namelist/HAMILTONIAN/file_ham,file_xyz,w90_with_soc,energy_thresh,use_degen_pert,&
         force_herm,force_antiherm,degen_thresh,apply_field,field_mode,Efield
      integer :: unit_inp

      open(newunit=unit_inp,file=trim(fname),status='OLD',action='READ')
      read(unit_inp,nml=HAMILTONIAN)
      close(unit_inp)

      me%file_ham = file_ham
      me%file_xyz = file_xyz
      me%w90_with_soc = w90_with_soc
      me%apply_field = apply_field
      me%field_mode = field_mode
      me%Efield = Efield
      me%energy_thresh = energy_thresh
      me%use_degen_pert = use_degen_pert
      me%force_herm = force_herm
      me%force_antiherm = force_antiherm
      me%degen_thresh = degen_thresh

   end subroutine Ham_ReadFromFile
!--------------------------------------------------------------------------------------

!====================================================================================== 
end module Mio_params
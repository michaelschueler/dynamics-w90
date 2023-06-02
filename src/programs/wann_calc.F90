program wann_calc
!! Computes various Berry-phase related quantities from a Wannier Hamiltonian.
!! The Hamiltonian can also be extended to a slab geometry.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,zero
   use scitools_time,only: Timer_act, Timer_Tic, Timer_Toc
   use scitools_utils,only: print_title, print_header
   use wan_latt_kpts,only: Read_Kpoints,kpoints_t
   use Mwannier_calc,only: wannier_calc_t
   use io_params,only: HamiltonianParams_t, WannierCalcParams_t
   use io_obs,only: WannierCalcOutput_t
   use io_hamiltonian,only: ReadHamiltonian, WriteHamiltonian
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! -- constants --
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   ! -- for reading i/o --
   logical :: PrintToFile = .false.
   integer :: Narg,unit_inp
   character(len=255) :: FlIn,FlOutPref
   ! -- internal variables --
   real(dp),allocatable,dimension(:,:,:)      :: orb_weight,spin,berry,oam
   real(dp),allocatable,dimension(:,:,:,:)    :: metric,spin_berry
   complex(dp),allocatable,dimension(:,:,:,:) :: velok
   type(kpoints_t)                            :: kp
   type(wannier_calc_t)                       :: wann
   type(HamiltonianParams_t)                  :: par_ham
   type(WannierCalcParams_t)                  :: par_calc
   type(WannierCalcOutput_t)                  :: output
!--------------------------------------------------------------------------------------
   call Timer_Act()
   call Timer_Tic('total',1)
   call print_title(output_unit,"Wannier90 post-processing")
!--------------------------------------------------------------------------------------
!                               ++  Read input ++
!--------------------------------------------------------------------------------------
   call Timer_Tic('Initialize Wannier', 2)
   Narg=command_argument_count()
   if(Narg>=1) then
      call get_command_argument(1,FlIn)
      call par_ham%ReadFromFile(FlIn)
      call par_calc%ReadFromFile(FlIn)
      call Read_Kpoints(FlIn,kp,print_info=.true.)
   else
      write(error_unit,fmt900) 'Please provide a namelist input file.'
      stop
   end if

   if(Narg>=2) then
      call get_command_argument(2,FlOutPref)
      PrintToFile=.true.
   end if
   if(.not.PrintToFile) then
      write(output_unit,fmt_info) 'No output prefix given. No output will be produced.'
   end if

   call wann%Init(par_ham,par_calc,kp)
 
   call output%AddEpsk(wann%epsk)
   if(par_calc%write_kpts) call output%AddKpts(kp%kpts)
   if(par_calc%calc_evecs) call output%AddEvecs(wann%vectk)

   write(output_unit,*)
   call Timer_Toc(N=2)
!--------------------------------------------------------------------------------------
!                               ++  Calculation ++
!--------------------------------------------------------------------------------------
   if(par_calc%calc_orbweight) then
      call Timer_Tic('orbital weight', 2)
      call wann%GetOrbitalWeight(orb_weight)
      call output%AddOrbweight(orb_weight)
      call Timer_Toc(N=2)
   end if

   if(par_calc%calc_spin) then
      call Timer_Tic('spin', 2)
      call wann%GetSpin(spin)
      call output%AddSpin(spin)
      call Timer_Toc(N=2)
   end if

   if(par_calc%calc_berry) then
      call Timer_Tic('Berry curvature', 2)
      call wann%GetBerryCurvature(berry,gauge=par_calc%gauge)
      call output%AddBerry(berry)
      call Timer_Toc(N=2)
   end if

   if(par_calc%calc_spin_berry) then
      call Timer_Tic('Spin-Berry curvature', 2)
      call wann%GetSpinBerryCurvature(spin_berry,gauge=par_calc%gauge)
      call output%AddSpinBerry(spin_berry)
      call Timer_Toc(N=2)
   end if

   if(par_calc%calc_oam) then
      call Timer_Tic('OAM', 2)
      call wann%GetOAM(oam,gauge=par_calc%gauge)
      call output%AddOAM(oam)
      call Timer_Toc(N=2)
   end if

   if(par_calc%calc_metric) then
      call Timer_Tic('quantum metric', 2)
      call wann%GetMetric(metric)
      call output%AddMetric(metric)
      call Timer_Toc(N=2)
   end if   

   if(par_calc%write_velocity) then
      call Timer_Tic('velocity matrix elements', 2)
      call wann%GetVelocity(velok)
      call output%AddVelocity(velok)
      call Timer_Toc(N=2)
   end if     

   if(PrintToFile) call output%SaveToFile(FlOutPref)
!--------------------------------------------------------------------------------------
   write(output_unit,*)
   call Timer_Toc(N=1)
   write(output_unit,fmt72)
!--------------------------------------------------------------------------------------

!======================================================================================
end program wann_calc

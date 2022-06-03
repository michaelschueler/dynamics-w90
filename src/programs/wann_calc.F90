program wann_calc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,zero
   use Mtime,only: Timer_act, Timer_Tic, Timer_Toc
   use Mutils,only: print_title, print_header
   use Mwannier_calc,only: wannier_calc_t
   use Mio_obs,only: WannierCalcOutput_t
   use Mio_hamiltonian,only: ReadHamiltonian, WriteHamiltonian
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! -- constants --
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   ! -- for reading i/o --
   logical :: PrintToFile = .false.
   integer :: Narg,unit_inp
   character(len=255) :: FlIn,FlOutPref
   ! -- input variables --
   logical :: slab_mode=.false.
   logical :: calc_orbweight=.false.
   logical :: calc_spin=.false.
   logical :: calc_berry=.false.
   logical :: calc_oam=.false.
   logical :: calc_evecs=.false.
   logical :: write_kpts=.false.
   integer :: gauge=0
   namelist/CALCOPT/slab_mode,calc_orbweight,calc_spin,calc_berry,calc_oam,calc_evecs,&
      write_kpts,gauge
   ! -- internal variables --
   real(dp),allocatable,dimension(:,:,:) :: orb_weight,spin,berry,oam
   type(wannier_calc_t)                  :: wann
   type(WannierCalcOutput_t)             :: output
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
      open(newunit=unit_inp,file=trim(FlIn),status='OLD',action='READ')
      read(unit_inp,nml=CALCOPT); rewind(unit_inp)
      close(unit_inp)
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

   call wann%Init(FlIn,slab_mode)
 
   call output%AddEpsk(wann%epsk)
   if(write_kpts) call output%AddKpts(wann%kpts)
   if(calc_evecs) call output%AddEvecs(wann%vectk)

   write(output_unit,*)
   call Timer_Toc(N=2)
!--------------------------------------------------------------------------------------
!                               ++  Calculation ++
!--------------------------------------------------------------------------------------
   if(calc_orbweight) then
      call Timer_Tic('orbital weight', 2)
      call wann%GetOrbitalWeight(orb_weight)
      call output%AddOrbweight(orb_weight)
      call Timer_Toc(N=2)
   end if

   if(calc_spin) then
      call Timer_Tic('spin', 2)
      call wann%GetSpin(spin)
      call output%AddSpin(spin)
      call Timer_Toc(N=2)
   end if

   if(calc_berry) then
      call Timer_Tic('Berry curvature', 2)
      call wann%GetBerryCurvature(berry)
      call output%AddBerry(berry)
      call Timer_Toc(N=2)
   end if

   if(calc_oam) then
      call Timer_Tic('OAM', 2)
      call wann%GetOAM(oam)
      call output%AddOAM(oam)
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

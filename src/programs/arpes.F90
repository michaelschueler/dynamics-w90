program arpes
!! Computes angle-resolved photoemission spectrum from Wannier functions.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,zero
   use scitools_time,only: Timer_act, Timer_Tic, Timer_Toc
   use scitools_utils,only: print_title, print_header
   use Mlatt_kpts,only: Read_Kpoints
   use Marpes_calc,only: arpes_calc_t
   use Mio_params,only: HamiltonianParams_t, PESParams_t
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
   real(dp),allocatable,dimension(:,:) :: kpts
   type(arpes_calc_t)                  :: calc
   type(HamiltonianParams_t)           :: par_ham
   type(PESParams_t)                   :: par_pes
!--------------------------------------------------------------------------------------
   call Timer_Act()
   call Timer_Tic('total',1)
   call print_title(output_unit,"Model ARPES from Wannier functions")
!--------------------------------------------------------------------------------------
!                               ++  Read input ++
!--------------------------------------------------------------------------------------
   call Timer_Tic('Initialize calculation', 2)
   Narg=command_argument_count()
   if(Narg >= 1) then
      call get_command_argument(1,FlIn)
      call par_ham%ReadFromFile(FlIn)
      call par_pes%ReadFromFile(FlIn)
      call Read_Kpoints(FlIn,kpts,print_info=.true.)
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

   call calc%Init(par_ham,par_pes,kpts)

   write(output_unit,*)
   call Timer_Toc(N=2)
!--------------------------------------------------------------------------------------
!                               ++  Calculation ++
!--------------------------------------------------------------------------------------
   call Timer_Tic('integrals', 2)
   call calc%CalcIntegrals()
   call Timer_Toc(N=2)

   call Timer_Tic('spectrum', 2)
   call calc%CalcPES()
   call Timer_Toc(N=2)

   ! .. output ...
   if(PrintToFile) call calc%WriteSpectrum(FlOutPref)
!--------------------------------------------------------------------------------------
   write(output_unit,*)
   call Timer_Toc(N=1)
   write(output_unit,fmt72)
!--------------------------------------------------------------------------------------

!======================================================================================
end program arpes

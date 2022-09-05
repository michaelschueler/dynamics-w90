program arpes_mpi
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use mpi
   use Mdebug
   use Mdef,only: dp,zero
   use Mtime,only: Timer_act, Timer_Tic, Timer_Toc, PrintTime
   use Mutils,only: print_title, print_header
   use Mlatt_kpts,only: Read_Kpoints
   use Marpes_calc_mpi,only: arpes_calc_t
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
   ! -- parallelization --
   logical  :: on_root
   integer  :: taskid,ntasks,nthreads,threadid,ierr
   real(dp) :: tic,toc
!--------------------------------------------------------------------------------------
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
   on_root = (taskid == 0)  
!--------------------------------------------------------------------------------------
   if(on_root) then 
      call Timer_Act()
      call Timer_Tic('total',1)
      call print_title(output_unit,"Model ARPES from Wannier functions")
   end if
!--------------------------------------------------------------------------------------
!                               ++  Read input ++
!--------------------------------------------------------------------------------------
   ! call Timer_Tic('Initialize calculation', 2)
   tic = MPI_Wtime()
   Narg=command_argument_count()
   if(Narg >= 1) then
      call get_command_argument(1,FlIn)
      call par_ham%ReadFromFile(FlIn)
      call par_pes%ReadFromFile(FlIn)
      call Read_Kpoints(FlIn,kpts,print_info=.true.,root_tag=on_root)
   else
      if(on_root) write(error_unit,fmt900) 'Please provide a namelist input file.'
      call MPI_Finalize(ierr); stop
   end if

   if(Narg>=2) then
      call get_command_argument(2,FlOutPref)
      PrintToFile=.true.
   end if
   if(.not.PrintToFile) then
      if(on_root) write(output_unit,fmt_info) 'No output prefix given. No output will be produced.'
   end if

   call calc%Init(par_ham,par_pes,kpts)

   toc = MPI_Wtime()
   if(on_root) then
      write(output_unit,*)
      call PrintTime('Initialize calculation',toc-tic,FlUnt=output_unit)
   end if
!--------------------------------------------------------------------------------------
!                               ++  Calculation ++
!--------------------------------------------------------------------------------------
   tic = MPI_Wtime()
   call calc%CalcIntegrals()
   toc = MPI_Wtime()
   if(on_root) call PrintTime('integrals',toc-tic,FlUnt=output_unit)

   tic = MPI_Wtime()
   call calc%CalcPES()
   toc = MPI_Wtime()
   if(on_root) call PrintTime('spectrum',toc-tic,FlUnt=output_unit)

   ! .. output ...
   if(PrintToFile) call calc%WriteSpectrum(FlOutPref)
!--------------------------------------------------------------------------------------
   if(on_root) then
      write(output_unit,*)
      call Timer_Toc(N=1)
      write(output_unit,fmt72)
   end if
   call MPI_Finalize(ierr)
!--------------------------------------------------------------------------------------

!======================================================================================
end program arpes_mpi

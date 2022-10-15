program wann_soc
!! Extends a Wannier Hamiltonian in spin space and adds atomic spin-orbit coupling (SOC).
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,zero
   use scitools_time,only: Timer_act, Timer_Tic, Timer_Toc
   use scitools_utils,only: print_title, print_header, get_file_ext, check_file_ext, str, &
      stop_error
   use Mham_w90,only: wann90_tb_t
   use Mwann_soc,only: ham_soc_t, AddSOC_Wannier
   use Mio_params,only: HamiltonianParams_t
   use Mio_hamiltonian,only: ReadHamiltonian, WriteHamiltonian, Read_SOC_Hamiltonian,&
      Read_SOC_lambda
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! -- constants --
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   ! -- for reading i/o --
   logical :: PrintToFile = .false.
   integer :: Narg,unit_inp
   character(len=255) :: FlIn,file_out
   ! -- input variables --
   type(HamiltonianParams_t)  :: par_ham
   ! -- internal variables --
   real(dp),allocatable       :: lam(:) 
   type(ham_soc_t)            :: soc
   type(wann90_tb_t)          :: w90_nosoc
   type(wann90_tb_t)          :: w90_soc
!--------------------------------------------------------------------------------------
   call Timer_Act()
   call Timer_Tic('total',1)
   call print_title(output_unit,"Wannier90 with atomic SOC")
!--------------------------------------------------------------------------------------
!                               ++  Read input ++
!--------------------------------------------------------------------------------------
   Narg=command_argument_count()
   if(Narg>=1) then
      call get_command_argument(1,FlIn)
      call par_ham%ReadFromFile(FlIn)
   else
      call stop_error('Please provide a namelist input file.')
   end if

   if(Narg>=2) then
      call get_command_argument(2,file_out)
      PrintToFile=.true.
   end if
   if(.not.PrintToFile) then
      write(output_unit,fmt_info) 'No output prefix given. No output will be produced.'
      write(output_unit,*)
   end if

   call ReadHamiltonian(par_ham%file_ham,w90_nosoc)

   call Read_SOC_lambda(par_ham%file_lam,lam)

   call Read_SOC_Hamiltonian(par_ham%file_soc,soc)

   call Timer_Tic('adding SOC', 2)
   call AddSOC_Wannier(soc,lam,w90_nosoc,w90_soc)
   call Timer_Toc(N=2)

   if(PrintToFile) then
      call WriteHamiltonian(w90_soc,file_out)
   end if
!--------------------------------------------------------------------------------------
   write(output_unit,*)
   call Timer_Toc(N=1)
   write(output_unit,fmt72)
!--------------------------------------------------------------------------------------

!======================================================================================
end program wann_soc

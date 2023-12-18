program wann_scissor
!! Creates a new Wannier Hamiltonian where a scissor operator is applied to the bands
!! 
!! ## How to run ##
!! Run the `wann_scissor` program by
!! ```
!!    ./exe/wann_scissor.x input_file output_file
!! ```
!! If no output file is provided as the second argument, the program will run, but
!! no output will be produced.
!!
!! ## Input variables ##
!! Input variables are read from `input_file` in Fortran name list format. The following
!! name list tags and variables are read:
!! ##### HAMILTONIAN #####
!! * `file_ham` : File name for reading the Hamiltonian. Can be hdf5 format (`*.h5`) or
!!                Wannier90 format (`*.dat` or `*.tb`)
!! * `scissor_index` : bands with index larger or equal `scissor_index` are shifted
!! * `scissor_energy`: Value of energy shift. Unit: atomic units.
!!
!! ## Output ##
!! The new Hamiltonian will be written to `output_file` (if given). The file extension
!! is used to determine the file format: `h5` triggers hdf5 output (if compiled with hdf5 support),
!! otherwise the Wannier90 plain text format is used.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,zero
   use scitools_time,only: Timer_act, Timer_Tic, Timer_Toc
   use scitools_utils,only: print_title, print_header, get_file_ext, check_file_ext, str
   use wan_hamiltonian,only: wann90_tb_t
   use wan_scissor,only: ScissorHamiltonian
   use io_params,only: HamiltonianParams_t
   use io_hamiltonian,only: ReadHamiltonian, WriteHamiltonian
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! -- constants --
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   ! -- for reading i/o --
   logical :: PrintToFile = .false.
   integer :: Narg,unit_inp
   character(len=255) :: FlIn,file_out
   ! -- internal variables --
   type(HamiltonianParams_t) :: par !! Collection of input variables
   type(wann90_tb_t) :: ham !! orginal Hamiltonian 
   type(wann90_tb_t) :: ham_sci !! Hamiltonian with scissor included
!--------------------------------------------------------------------------------------
   call Timer_Act()
   call Timer_Tic('total',1)
   call print_title(output_unit,"Wannier90 scissor operation")
!--------------------------------------------------------------------------------------
!                               ++  Read input ++
!--------------------------------------------------------------------------------------
   Narg=command_argument_count()
   if(Narg>=1) then
      call get_command_argument(1,FlIn)
      call par%ReadFromFile(FlIn)
   else
      write(error_unit,fmt900) 'Please provide a namelist input file.'
      stop
   end if

   if(Narg>=2) then
      call get_command_argument(2,file_out)
      PrintToFile=.true.
   end if
   if(.not.PrintToFile) then
      write(output_unit,fmt_info) 'No output prefix given. No output will be produced.'
      write(output_unit,*)
   end if

   if(len_trim(par%file_xyz) > 0) then
      call ReadHamiltonian(par%file_ham,Ham,file_xyz=par%file_xyz)
   else
      call ReadHamiltonian(par%file_ham,Ham)
   end if

   call ScissorHamiltonian(par%scissor_index,par%scissor_energy,Ham,Ham_sci)

   if(PrintToFile) then
      call WriteHamiltonian(Ham_sci,file_out)
   end if
!--------------------------------------------------------------------------------------
   write(output_unit,*)
   call Timer_Toc(N=1)
   write(output_unit,fmt72)
!--------------------------------------------------------------------------------------

!======================================================================================
end program wann_scissor

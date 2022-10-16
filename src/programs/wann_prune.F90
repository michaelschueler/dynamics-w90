program wann_prune
!! Compresses a Wannier Hamiltonian by cutting small hopping amplitudes. 
!! Output can be written in original Wannier90 or in hdf5 format.
!! 
!! ## How to run ##
!! Run the `wann_prune` program by
!! ```
!!    ./exe/wann_prune.x input_file output_file
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
!! * `file_xyz` : File name for reading the Wannier centers. Only relevant if the 
!!                Hamiltonian is given in Wannier90 format. If not provided,
!!                the Wannier centers are not read and the output Hamiltonian will not
!!                contain the Wannier centers.
!! * `energy_thresh`: Threshold for cutting hopping amplitudes. Unit: atomic units.
!!
!! ## Output ##
!! The pruned Hamiltonian will be written to `output_file` (if given). The file extension
!! is used to determine the file format: `h5` triggers hdf5 output (if compiled with hdf5 support),
!! otherwise the Wannier90 plain text format is used.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,zero
   use scitools_time,only: Timer_act, Timer_Tic, Timer_Toc
   use scitools_utils,only: print_title, print_header, get_file_ext, check_file_ext, str
   use wan_hamiltonian,only: wann90_tb_t
   use wan_compress,only: PruneHoppings
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
   real(dp) :: comp_rate !! compression rate
   type(HamiltonianParams_t) :: par !! Collection of input variables
   type(wann90_tb_t) :: wann !! orginal Hamiltonian 
   type(wann90_tb_t) :: wann_pruned !! Hamiltonian with pruned hoppings
!--------------------------------------------------------------------------------------
   call Timer_Act()
   call Timer_Tic('total',1)
   call print_title(output_unit,"Wannier90 compression")
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
      call ReadHamiltonian(par%file_ham,wann,file_xyz=par%file_xyz)
   else
      call ReadHamiltonian(par%file_ham,wann)
   end if

   call PruneHoppings(par%energy_thresh,wann,wann_pruned,comp_rate)

   write(output_unit,fmt148) "hoppings (before compression)", wann%nrpts
   write(output_unit,fmt148) "hoppings (after compression)", wann_pruned%nrpts
   write(output_unit,'(a)') "compression rate: "//str(nint(100 * comp_rate)) // "%"

   if(PrintToFile) then
      call WriteHamiltonian(wann_pruned,file_out)
   end if
!--------------------------------------------------------------------------------------
   write(output_unit,*)
   call Timer_Toc(N=1)
   write(output_unit,fmt72)
!--------------------------------------------------------------------------------------

!======================================================================================
end program wann_prune

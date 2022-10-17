program wann_soc
!! Extends a Wannier Hamiltonian in spin space and adds atomic spin-orbit coupling (SOC).
!! This is accomplished by adding an on-site term to the Wannier Hamiltonian 
!! \(H_{m\mathbf{R},m^\prime\mathbf{R}^\prime}\).
!!
!! ## Description ##
!! $$H^\mathrm{SOC}_{m\sigma\mathbf{R},m^\prime\sigma^\prime\mathbf{R}^\prime} = H_{m\mathbf{R},m^\prime\mathbf{R}^\prime} 
!! \delta_{\sigma \sigma^\prime}
!!  + \delta_{\mathbf{R},\mathbf{R}^\prime} \sum^{n_g}_{i=1} \lambda_i 
!! \langle m \sigma | \mathbf{L}^{(i)}\cdot \mathbf{S} | m^\prime \sigma^\prime \rangle .$$
!!
!! Here, we assume that the orbitals \(m\) of the original Hamiltonian are grouped into \(n_g\) groups that are represented
!! by the same orbital angular momentum operator \(\mathbf{L}^{(i)}\). Typically this corresponds to a including
!! all 3 p, 5 d, ... orbitals per atom in the Wannier basis.
!!
!! The program reads the original Wannier Hamiltonian, a file containg the orbital angular momentum operators
!!  \(\mathbf{L}^{(i)}\), and a plain text file with the \(n_g\) values of \(\lambda_i\). The helper script
!! `SOCInput.py` (for hdf5 format) or `SOCInput_txt.py` (plain text format) in the `python_utils` directory
!! can be used to create the corresponding input files.
!! 
!! ## How to run ##
!! Run the `wann_soc` program by
!! ```
!!    ./exe/wann_soc.x input_file output_file
!! ```
!! If no output file is provided as the second argument, the program will run, but
!! no output will be produced.
!!
!! ## Input variables ##
!! Input variables are read from `input_file` in Fortran name list format. The following
!! name list tags and variables are read:
!! ##### HAMILTONIAN #####
!! * `file_ham` : File name for reading the Hamiltonian without SOC. Can be hdf5 format (`*.h5`) or
!!                Wannier90 format (`*.dat` or `*.tb`)
!! * `file_xyz` : File name for reading the Wannier centers. Only relevant if the 
!!                Hamiltonian is given in Wannier90 format. If not provided,
!!                the Wannier centers are not read and the output Hamiltonian will not
!!                contain the Wannier centers.
!! * `file_lam` : File name for reading the SOC parameters \(\lambda_i\) from plain text format.
!! * `file_soc` : File name for reading the orbital angular momentum operators \(\mathbf{L}^{(i)}\). 
!!                hdf5 or plain text format.
!!
!! ## Output ##
!! The Hamiltonian including SOC will be written to `output_file` (if given). The file extension
!! is used to determine the file format: `h5` triggers hdf5 output (if compiled with hdf5 support),
!! otherwise the Wannier90 plain text format is used.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,zero
   use scitools_time,only: Timer_act, Timer_Tic, Timer_Toc
   use scitools_utils,only: print_title, print_header, get_file_ext, check_file_ext, str, &
      stop_error
   use wan_hamiltonian,only: wann90_tb_t
   use wan_soc,only: ham_soc_t, AddSOC_Wannier
   use io_params,only: HamiltonianParams_t
   use io_hamiltonian,only: ReadHamiltonian, WriteHamiltonian, Read_SOC_Hamiltonian,&
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

   if(len_trim(par%file_xyz) > 0) then
      call ReadHamiltonian(par%file_ham,w90_nosoc,file_xyz=par%file_xyz)
   else
      call ReadHamiltonian(par%file_ham,w90_nosoc)
   end if

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

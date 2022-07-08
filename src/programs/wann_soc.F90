program wann_soc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,zero
   use Mtime,only: Timer_act, Timer_Tic, Timer_Toc
   use Mutils,only: print_title, print_header, get_file_ext, check_file_ext, str
   use Mham_w90,only: wann90_tb_t
   use Mwann_soc,only: ham_soc_t, AddSOC_Wannier
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
   character(len=255) :: file_ham
   logical            :: w90_with_soc,expert_params
   real(dp)           :: energy_thresh
   namelist/HAMILTONIAN/file_ham,w90_with_soc,energy_thresh,expert_params
   character(len=255) :: file_lam,file_soc
   namelist/SPINORBIT/file_lam,file_soc
   ! -- internal variables --
   real(dp),allocatable :: lam(:) 
   type(ham_soc_t)      :: soc
   type(wann90_tb_t)    :: w90_nosoc
   type(wann90_tb_t)    :: w90_soc
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
      open(newunit=unit_inp,file=trim(FlIn),status='OLD',action='READ')
      read(unit_inp,nml=HAMILTONIAN); rewind(unit_inp)
      read(unit_inp,nml=SPINORBIT)
      close(unit_inp)
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

   call ReadHamiltonian(file_ham,w90_nosoc)

   call Read_SOC_lambda(file_lam,lam)

   call Read_SOC_Hamiltonian(file_soc,soc)

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

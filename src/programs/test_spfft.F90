program test_spfft
!! Test program for SpFFT implementation
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,zero
   use scitools_time,only: Timer_act, Timer_Tic, Timer_Toc
   use scitools_utils,only: print_title, print_header, stop_error
   use wan_latt_kpts,only: Read_Kpoints,kpoints_t
   use wan_hamiltonian,only: wann90_tb_t
   use wan_fft_ham,only: wann_fft_t
   use wan_spfft_ham,only: wann_spfft_t
   use io_params,only: HamiltonianParams_t
   use io_hamiltonian,only: ReadHamiltonian, WriteHamiltonian
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! -- constants --
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   character(len=*),parameter :: fmt_input='(" Input: [",a,"]")'
   ! -- for reading i/o --
   logical :: PrintToFile = .false.
   integer :: Narg,unit_inp
   character(len=255) :: FlIn,FlOutPref
   ! -- input variables --
   integer :: nthreads_fft=1,nthreads_orb=1
   namelist/PARALLELIZATION/nthreads_fft,nthreads_orb
   ! -- internal variables --
   logical :: file_ok
   integer :: nbnd,ik
   real(dp) :: kpt(3),Ared(3)
   complex(dp),allocatable    :: Hk_dft(:,:,:), Hk_fft(:,:,:), Hk_spfft(:,:,:)
   type(kpoints_t)            :: kp
   type(HamiltonianParams_t)  :: par_ham
   type(wann90_tb_t)          :: ham
   type(wann_fft_t)           :: ham_ft
   type(wann_spfft_t)         :: ham_sp
!--------------------------------------------------------------------------------------
   call Timer_Act()
   call Timer_Tic('total',1)
   call print_title(output_unit,"Test: SpFFT")
!--------------------------------------------------------------------------------------
!                               ++  Read input ++
!--------------------------------------------------------------------------------------
   call Timer_Tic('Initialize Wannier', 2)
   Narg=command_argument_count()
   if(Narg>=1) then
      call get_command_argument(1,FlIn)
      call par_ham%ReadFromFile(FlIn)
      call Read_Kpoints(FlIn,kp,print_info=.true.)
      open(newunit=unit_inp,file=trim(FlIn),STATUS='OLD',ACTION='READ')
      read(unit_inp,nml=PARALLELIZATION)
      close(unit_inp)
   else
      write(error_unit,fmt900) 'Please provide a namelist input file.'
      stop
   end if

   inquire(file=trim(par_ham%file_ham),exist=file_ok)
   if(.not.file_ok) then
      call stop_error('Input file does not exist: '//trim(par_ham%file_ham))
   end if
   write(output_unit,fmt_input) 'Hamiltionian from file: '//trim(par_ham%file_ham)
   call ReadHamiltonian(par_ham%file_ham,Ham)   

   write(output_unit,*)
   call Timer_Toc(N=2)
!--------------------------------------------------------------------------------------
   nbnd = Ham%num_wann
   allocate(Hk_dft(nbnd,nbnd,kp%Nk))
   allocate(Hk_fft(nbnd,nbnd,kp%Nk))
   allocate(Hk_spfft(nbnd,nbnd,kp%Nk))

   call ham_ft%InitFromW90(Ham, [kp%nk1, kp%nk2, kp%nk3], &
      nthreads_fft=nthreads_fft,nthreads_orb=nthreads_orb)
   call ham_sp%InitFromW90(Ham, [kp%nk1, kp%nk2, kp%nk3], nthreads=nthreads_fft)

   Ared = [0.1_dp, 0.0_dp, -0.4_dp]

   call Timer_Tic('DFT', 2)

   !$OMP PARALLEL PRIVATE(kpt)
   !$OMP DO
   do ik=1,kp%Nk
      kpt = kp%kpts(ik,:) - Ared
      Hk_dft(:,:,ik) = Ham%get_ham(kpt)
   end do
   !$OMP END DO
   !$OMP END PARALLEL

   write(output_unit,*)
   call Timer_Toc(N=2)

   call Timer_Tic('FFT', 2) 

   call Ham_ft%GetHam_Dressed(Ared,Hk_fft)

   write(output_unit,*)
   call Timer_Toc(N=2)

   call Timer_Tic('SpFFT', 2) 

   call Ham_sp%GetHam(Hk_spfft, Ar=Ared)

   write(output_unit,*)
   call Timer_Toc(N=2)

   ! do ik=1,kp%Nk
   !    print*, Hk_dft(1,1,ik), Hk_fft(1,1,ik), Hk_spfft(1,1,ik)
   ! end do

!--------------------------------------------------------------------------------------
   write(output_unit,*)
   call Timer_Toc(N=1)
   write(output_unit,fmt72)
!--------------------------------------------------------------------------------------


!======================================================================================
end program test_spfft
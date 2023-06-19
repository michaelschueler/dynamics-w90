program test_spfft
!! Test program for SpFFT implementation
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use mpi
   use Mdebug
   use scitools_def,only: dp,zero
   use scitools_time,only: Timer_act, Timer_Tic, Timer_Toc
   use scitools_utils,only: print_title, print_header, stop_error
   use scitools_array1d_dist,only: dist_array1d_t
   use wan_latt_kpts,only: Read_Kpoints,kpoints_t
   use wan_hamiltonian,only: wann90_tb_t
   use wan_fft_ham_mpi,only: wann_fft_t
   use wan_spfft_ham_mpi,only: wann_spfft_t
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
   integer :: nbnd,Nk_loc,ik,ik_glob
   real(dp) :: kpt(3),Ared(3)
   complex(dp),allocatable    :: Hk_dft(:,:,:), Hk_fft(:,:,:), Hk_spfft(:,:,:)
   type(dist_array1d_t)       :: kdist
   type(kpoints_t)            :: kp
   type(HamiltonianParams_t)  :: par_ham
   type(wann90_tb_t)          :: ham
   type(wann_fft_t)           :: ham_ft
   type(wann_spfft_t)         :: ham_sp
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
      call print_title(output_unit,"Test: SpFFT")
   end if
!--------------------------------------------------------------------------------------
!                               ++  Read input ++
!--------------------------------------------------------------------------------------
   tic = MPI_Wtime()
   Narg=command_argument_count()
   if(Narg>=1) then
      call get_command_argument(1,FlIn)
      call par_ham%ReadFromFile(FlIn)
      call Read_Kpoints(FlIn,kp,print_info=.true.,root_tag=on_root)
      open(newunit=unit_inp,file=trim(FlIn),STATUS='OLD',ACTION='READ')
      read(unit_inp,nml=PARALLELIZATION)
      close(unit_inp)
   else
      call stop_error('Please provide a namelist input file.', root_flag=on_root)
   end if

   inquire(file=trim(par_ham%file_ham),exist=file_ok)
   if(.not.file_ok) then
      call stop_error('Input file does not exist: '//trim(par_ham%file_ham), root_flag=on_root)
   end if
   if(on_root) write(output_unit,fmt_input) 'Hamiltionian from file: '//trim(par_ham%file_ham)
   call ReadHamiltonian(par_ham%file_ham,Ham)   

   toc = MPI_Wtime()

   if(on_root) then
      write(output_unit,*)
      write(output_unit,fmt555) "Reading input",toc-tic
   end if
!--------------------------------------------------------------------------------------
   nbnd = Ham%num_wann

   call kdist%Init(ntasks,taskid,kp%Nk)
   Nk_loc = kdist%N_loc(taskid)

   allocate(Hk_dft(nbnd,nbnd,Nk_loc))
   allocate(Hk_fft(nbnd,nbnd,Nk_loc))
   allocate(Hk_spfft(nbnd,nbnd,Nk_loc))

   call ham_ft%InitFromW90(Ham, [kp%nk1, kp%nk2, kp%nk3], kdist)
   call ham_sp%InitFromW90(Ham, [kp%nk1, kp%nk2, kp%nk3], kdist, nthreads=nthreads_fft)

   Ared = [0.1_dp, 0.0_dp, -0.4_dp]

   tic = MPI_Wtime()

   !$OMP PARALLEL PRIVATE(kpt,ik_glob)
   !$OMP DO
   do ik=1,Nk_loc
      ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
      kpt = kp%kpts(ik_glob,:) - Ared
      Hk_dft(:,:,ik) = Ham%get_ham(kpt)
   end do
   !$OMP END DO
   !$OMP END PARALLEL

   toc = MPI_Wtime()
   if(on_root) then
      write(output_unit,*)
      write(output_unit,fmt555) "DFT",toc-tic
   end if

   tic = MPI_Wtime()
   call Ham_ft%GetHam(kdist, Hk_fft, Ar=Ared)

   toc = MPI_Wtime()
   if(on_root) then
      write(output_unit,*)
      write(output_unit,fmt555) "FFT",toc-tic
   end if

   tic = MPI_Wtime()
   call Ham_sp%GetHam(kdist, Hk_spfft, Ar=Ared)

   toc = MPI_Wtime()
   if(on_root) then
      write(output_unit,*)
      write(output_unit,fmt555) "SpFFT",toc-tic
   end if

   ! if(on_root) then
   !    do ik=1,Nk_loc
   !       print*, Hk_dft(1,1,ik), Hk_fft(1,1,ik), Hk_spfft(1,1,ik)
   !    end do
   ! end if

!--------------------------------------------------------------------------------------
   if(on_root) then
      write(output_unit,*)
      call Timer_Toc(N=1)
      write(output_unit,fmt72)
   end if
   call MPI_Finalize(ierr)
!--------------------------------------------------------------------------------------


!======================================================================================
end program test_spfft
program test_fftw_mpi
!! Test program for SpFFT implementation
!======================================================================================
   use,intrinsic :: iso_fortran_env,only: output_unit,error_unit
   use,intrinsic :: ISO_C_binding
   use mpi
   use Mdebug
   use scitools_def,only: dp,zero
   use scitools_time,only: Timer_act, Timer_Tic, Timer_Toc
   use scitools_utils,only: print_title, print_header, stop_error
   use scitools_array1d_dist,only: dist_array1d_t
   use wan_latt_kpts,only: Read_Kpoints,kpoints_t
   use wan_hamiltonian,only: wann90_tb_t
   use wan_fft_ham_mpi,only: wann_fft_t
   use io_params,only: HamiltonianParams_t
   use io_hamiltonian,only: ReadHamiltonian, WriteHamiltonian
   implicit none
   include '../formats.h'
   include 'fftw3-mpi.f03'
!--------------------------------------------------------------------------------------
   ! -- constants --
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   character(len=*),parameter :: fmt_input='(" Input: [",a,"]")'
   integer,parameter :: kp_list=0, kp_path=1, kp_grid=2, kp_fft_grid_2d=3, kp_fft_grid_3d=4
   ! -- for reading i/o --
   logical :: PrintToFile = .false.
   integer :: Narg,unit_inp
   character(len=255) :: FlIn,FlOutPref
   ! -- internal variables --
   logical :: file_ok
   integer :: nbnd,Nk_loc,ik,ik_glob
   real(dp) :: kpt(3),Ared(3)
   complex(dp),allocatable    :: Hk_dft(:,:,:), Hk_fft(:,:,:)
   type(kpoints_t)            :: kp
   type(HamiltonianParams_t)  :: par_ham
   type(wann90_tb_t)          :: ham
   type(wann_fft_t)           :: ham_ft
   type(dist_array1d_t)       :: kdist
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
      call print_title(output_unit,"Test: FFTW_MPI")
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
   else
      call stop_error('Please provide a namelist input file.', root_flag=on_root)
   end if

   inquire(file=trim(par_ham%file_ham),exist=file_ok)
   if(.not.file_ok) then
      call stop_error('Input file does not exist: '//trim(par_ham%file_ham), root_flag=on_root)
   end if
   if(on_root) write(output_unit,fmt_input) 'Hamiltionian from file: '//trim(par_ham%file_ham)
   call ReadHamiltonian(par_ham%file_ham,Ham)   

   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   toc = MPI_Wtime()

   if(on_root) then
      write(output_unit,*)
      write(output_unit,fmt555) "Reading input",toc-tic
   end if
!--------------------------------------------------------------------------------------
   if(kp%kpoints_type == kp_fft_grid_2d) then
      call kdist%Init(ntasks,taskid,kp%nk2,dist_scheme=1,blocksize=kp%nk1)
   elseif(kp%kpoints_type == kp_fft_grid_3d) then
      call kdist%Init(ntasks,taskid,kp%nk3,dist_scheme=1,blocksize=kp%nk1*kp%nk2)
   else
      call stop_error("FFTW+MPI only available for FFT grids in 2D and 3D", root_flag=on_root)
   end if

   Nk_loc = kdist%N_loc(taskid)

   nbnd = Ham%num_wann
   allocate(Hk_dft(nbnd,nbnd,Nk_loc))
   allocate(Hk_fft(nbnd,nbnd,Nk_loc))

   call ham_ft%InitFromW90(Ham, [kp%nk1, kp%nk2, kp%nk3], kdist)

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   tic = MPI_Wtime()

   Ared = [0.0_dp, 0.0_dp, 0.0_dp]

   do ik=1,Nk_loc
      ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
      kpt = kp%kpts(ik_glob,:) - Ared
      Hk_dft(:,:,ik) = Ham%get_ham(kpt)
   end do

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   toc = MPI_Wtime()
   if(on_root) then
      write(output_unit,*)
      write(output_unit,fmt555) "DFT",toc-tic
   end if

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   tic = MPI_Wtime()
   call Ham_ft%GetHam(kdist, Hk_fft, Ar=Ared)

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   toc = MPI_Wtime()
   if(on_root) then
      write(output_unit,*)
      write(output_unit,fmt555) "FFT",toc-tic
   end if

   ! if(on_root) then
   !    do ik=1,Nk_loc
   !       print*, ik, Hk_dft(1,1,ik), Hk_fft(1,1,ik) 
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
end program test_fftw_mpi
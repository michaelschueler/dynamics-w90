program wann_evol_mpi
!! Computes the time-dependent dynamics of an electron system described by a 
!! Wannier Hamiltonian upon laser excitation and computes various observables.
!! MPI parallelized version.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use omp_lib
   use mpi
   use Mdebug
   use Mdef,only: dp,zero
   use Mtime,only: Timer_Act, Timer_Tic, Timer_Toc
   use Mutils,only: print_title, print_header, get_file_ext, check_file_ext
   use Mlaserpulse,only: LaserPulse_spline_t
   use Mham_w90,only: wann90_tb_t
   use Mwann_evol_mpi,only: wann_evol_t
   use Mlatt_kpts,only: Read_Kpoints
   use Mio_hamiltonian,only: ReadHamiltonian
   use Mio_obs,only: SaveTDObs, SaveTDOccupation
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! -- constants --
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   character(len=*),parameter :: fmt_input='(" Input: [",a,"]")'
   integer,parameter :: velocity_gauge=0,dipole_gauge=1,velo_emp_gauge=2,dip_emp_gauge=3
   ! -- for reading i/o --
   logical :: PrintToFile = .false.
   integer :: Narg,unit_inp
   character(len=255) :: FlIn,FlOutPref
   ! -- input variables --
   logical            :: Output_Occ_KPTS=.false.
   logical            :: FixMuChem=.true.
   integer            :: gauge
   real(dp)           :: MuChem,Beta,Filling
   character(len=255) :: file_ham
   namelist/SYSPARAMS/MuChem,FixMuChem,Filling,Beta,file_ham,gauge,Output_Occ_KPTS
   !......................................
   logical  :: relaxation_dynamics=.false.
   integer  :: Nt,output_step=1
   real(dp) :: Tmax,tstart=0.0_dp,T1_relax=1.0e10_dp,T2_relax=1.0e10_dp
   character(len=255) :: file_field=""
   namelist/TIMEPARAMS/Nt,Tmax,output_step,tstart,relaxation_dynamics,T1_relax,T2_relax,&
      file_field
   !......................................
   ! -- internal variables --
   logical  :: file_ok
   logical  :: ApplyField=.false.
   integer  :: Nsteps,tstp,tstp_max,step
   real(dp) :: dt,pulse_tmin,pulse_tmax
   real(dp),allocatable,dimension(:)     :: Etot,Ekin
   real(dp),allocatable,dimension(:,:)   :: kpts,BandOcc,Jcurr,Dip
   real(dp),allocatable,dimension(:,:)   :: Jpara,Jdia,JHk,Jpol,Jintra
   real(dp),allocatable,dimension(:,:,:) :: Occk
   type(LaserPulse_spline_t) :: pulse_x,pulse_y,pulse_z
   type(wann90_tb_t)         :: Ham
   type(wann_evol_t)         :: lattsys
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
      call print_title(output_unit,"Wannier90 dynamics")  

      !$OMP PARALLEL PRIVATE(nthreads,threadid)
      threadid = omp_get_thread_num()
      if(threadid == 0) then
         nthreads = omp_get_num_threads()
         write(output_unit,fmt148) 'number of threads',nthreads
         write(output_unit,*)
      end if
      !$OMP END PARALLEL
   end if
!--------------------------------------------------------------------------------------
!                               ++  Read input ++
!--------------------------------------------------------------------------------------
   tic = MPI_Wtime()

   Narg=command_argument_count()
   if(Narg>=1) then
      call get_command_argument(1,FlIn)
      open(newunit=unit_inp,file=trim(FlIn),status='OLD',action='READ')
      read(unit_inp,nml=SYSPARAMS) ; rewind(unit_inp)
      read(unit_inp,nml=TIMEPARAMS)
      close(unit_inp)
   else
      if(on_root) write(output_unit,fmt900) 'Please provide a namelist input file.'
      call MPI_Finalize(ierr); stop
   end if

   if(Narg>=2) then
      call get_command_argument(2,FlOutPref)
      PrintToFile=.true.
   end if
   if((.not.PrintToFile) .and. on_root) then
      write(output_unit,fmt_info) 'No output prefix given. No output will be produced.'
   end if

   inquire(file=trim(file_ham),exist=file_ok)
   if(.not.file_ok) then
      if(on_root) write(error_unit,fmt900) 'Input file does not exist: '//trim(file_ham)
      call MPI_Finalize(ierr); stop
   end if
   if(on_root) write(output_unit,fmt_input) 'Hamiltionian from file: '//trim(file_ham)
   call ReadHamiltonian(file_ham,Ham)

   call Read_Kpoints(FlIn,kpts,print_info=.true.,root_tag=on_root)

   ApplyField = len_trim(file_field) > 0
   if(ApplyField) then
      inquire(file=trim(file_field),exist=file_ok)
      if(.not.file_ok) then
         if(on_root) write(error_unit,fmt900) 'Input file does not exist: '//trim(file_field)
         call MPI_Finalize(ierr); stop
      end if
      if(on_root) write(output_unit,fmt_input) 'External field from file: '//trim(file_field)
      call pulse_x%Load_ElectricField(file_field,usecol=2,CalcAfield=.true.)
      call pulse_y%Load_ElectricField(file_field,usecol=3,CalcAfield=.true.)
      call pulse_z%Load_ElectricField(file_field,usecol=4,CalcAfield=.true.)
   end if

   toc = MPI_Wtime()

   if(on_root) then
      write(output_unit,*)
      write(output_unit,fmt555) "Reading input",toc-tic
   end if
!--------------------------------------------------------------------------------------
!                               ++  Equilibrium ++
!--------------------------------------------------------------------------------------
   if(on_root) then
      call print_header(output_unit,"Equilibrium","*")
   end if

   tic = MPI_Wtime()

   call lattsys%Init(Beta,MuChem,ham,kpts,gauge)
   call lattsys%SetLaserPulse(external_field)

   if(FixMuChem) then
      call lattsys%SolveEquilibrium()
      if(on_root) write(output_unit,fmt149) "number of electrons", lattsys%nelec
   else
      call lattsys%SolveEquilibrium(Filling)
      if(on_root) write(output_unit,fmt149) "chemical potential", lattsys%MuChem
   end if
   toc = MPI_Wtime()

   if(on_root) then
      write(output_unit,*)
      write(output_unit,fmt555) "equilibrium",toc-tic
      write(output_unit,*)
      write(output_unit,fmt_info,advance='no') "gauge: "
      select case(gauge)
      case(dipole_gauge)
         write(output_unit,*) "dipole gauge"
      case(dip_emp_gauge)
         write(output_unit,*) "Peierls substitution"
      case(velocity_gauge)
         write(output_unit,*) "velocity gauge"
      case(velo_emp_gauge)
         write(output_unit,*) "empirical velocity gauge"
      case default
         write(error_unit,fmt900) "unrecognized gauge"
      end select
   end if
!--------------------------------------------------------------------------------------
!                         ++  Time propagation ++
!--------------------------------------------------------------------------------------
   if(on_root) then
      call print_header(output_unit,"Time propagation","*")
   end if

   tic = MPI_Wtime()

   Nsteps = ceiling((Nt+1)/dble(output_step)) - 1

   dt = Tmax/dble(Nt)

   allocate(Ekin(0:Nsteps),Etot(0:Nsteps),Jcurr(3,0:Nsteps),Dip(3,0:Nsteps))
   allocate(BandOcc(Ham%num_wann,0:Nsteps))
   if(gauge == dipole_gauge .or. gauge == dip_emp_gauge) then
      allocate(JHk(3,0:Nsteps),Jpol(3,0:Nsteps))
   else
      allocate(Jpara(3,0:Nsteps),Jdia(3,0:Nsteps),Jintra(3,0:Nsteps))
   end if

   if(Output_Occ_KPTS) then
      allocate(Occk(lattsys%nbnd,lattsys%Nk,0:Nsteps))
   end if

   if(gauge == dipole_gauge .or. gauge == dip_emp_gauge) then
      call lattsys%CalcObservables_dip(0,dt,Ekin(0),Etot(0),Jcurr(:,0),JHk(:,0),Jpol(:,0),&
               Dip(:,0),BandOcc(:,0))
   else
      call lattsys%CalcObservables_velo(0,dt,Ekin(0),Etot(0),Jcurr(:,0),Jpara(:,0),Jdia(:,0),&
         Jintra(:,0),Dip(:,0),BandOcc(:,0))
   end if

   if(Output_Occ_KPTS) then
      call lattsys%GetOccupationKPTS(Occk(:,:,0))
   end if

   pulse_tmin = min(pulse_x%Tmin,pulse_y%Tmin,pulse_z%Tmin) - tstart
   pulse_tmax = max(pulse_x%Tmax,pulse_y%Tmax,pulse_z%Tmax) - tstart
   if(.not. ApplyField) pulse_tmax = 0.0_dp

   step = 0
   do tstp=0,Nt-1
      if(relaxation_dynamics) then
         call lattsys%Timestep_RelaxTime(T1_relax,T2_relax,tstp,dt)
      else
         call lattsys%Timestep(tstp,dt,field_tmax=pulse_tmax)
      end if

      if(mod(tstp+1, output_step) == 0) then
         step = step + 1

         if(gauge == dipole_gauge .or. gauge == dip_emp_gauge) then
            call lattsys%CalcObservables_dip(tstp+1,dt,Ekin(step),Etot(step),Jcurr(:,step),JHk(:,step),&
                  Jpol(:,step),Dip(:,step),BandOcc(:,step))
         else
            call lattsys%CalcObservables_velo(tstp+1,dt,Ekin(step),Etot(step),Jcurr(:,step),Jpara(:,step),&
                  Jdia(:,step),Jintra(:,step),Dip(:,step),BandOcc(:,step))
         end if

         if(Output_Occ_KPTS) then
            call lattsys%GetOccupationKPTS(Occk(:,:,step))
         end if

         if(on_root) write(output_unit,fmt145) "tstp",tstp+1

      end if

   end do

   toc = MPI_Wtime()
   if(on_root) then
      write(output_unit,fmt50)
      write(output_unit,fmt555) 'propagation',toc-tic
   end if
!--------------------------------------------------------------------------------------
   if(on_root .and. PrintToFile) then
      select case(gauge)
      case(velocity_gauge, velo_emp_gauge)
         call SaveTDObs(FlOutPref,Nsteps,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
            Jpara,Jdia,Jintra)
      case(dipole_gauge, dip_emp_gauge)
         call SaveTDObs(FlOutPref,Nsteps,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
            JHk,Jpol)         
      end select

      if(Output_Occ_KPTS) then
         call SaveTDOccupation(FlOutPref,Nsteps,output_step,dt,Occk)
      end if
   end if
!--------------------------------------------------------------------------------------
   if(on_root) then
      write(output_unit,*)
      call Timer_Toc(N=1)
      write(output_unit,fmt72)
   end if
   call MPI_Finalize(ierr)
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine external_field(t,AF,EF)
      real(dp),intent(in)   :: t
      real(dp),intent(out) :: AF(3),EF(3)

      AF = 0.0_dp; EF = 0.0_dp

      if(ApplyField) then
         AF(1) = pulse_x%Afield(t+tstart)
         AF(2) = pulse_y%Afield(t+tstart)
         AF(3) = pulse_z%Afield(t+tstart)
         EF(1) = pulse_x%Efield(t+tstart)
         EF(2) = pulse_y%Efield(t+tstart)
         EF(3) = pulse_z%Efield(t+tstart)
      end if

   end subroutine external_field
!--------------------------------------------------------------------------------------

!======================================================================================
end program wann_evol_mpi

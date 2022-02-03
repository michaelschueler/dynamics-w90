program wann_evol_mpi
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use omp_lib
   use mpi
   use Mdebug
   use Mdef,only: dp,zero
   use Mtime,only: Timer_Act, Timer_SetName, Timer_Run, Timer_stop, Timer_DTShow
   use Mlaserpulse,only: LaserPulse_spline_t
   use Mham_w90,only: wann90_tb_t
   use Mwann_evol_mpi,only: wann_evol_t
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! -- constants --
   character(len=*),parameter :: fmt_time='(a,"/",a,"/",a," at ",a,":",a,":",a)'
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   character(len=*),parameter :: fmt_input='(" Input: [",a,"]")'
   integer,parameter :: velocity_gauge=0,dipole_gauge=1,velo_emp_gauge=2,dip_emp_gauge=3
   ! -- timing --
   character(len=8) :: date
   character(len=10) :: time
   ! -- for reading i/o --
   logical :: PrintToFile = .false.
   integer :: Narg,unit_inp
   character(len=255) :: FlIn,FlOutPref
   ! -- input variables --
   logical            :: Output_Dens=.false.
   logical            :: FixMuChem=.true.
   integer            :: gauge
   integer            :: Nk1,Nk2,Nk3
   real(dp)           :: MuChem,Beta,Filling
   character(len=255) :: file_ham,file_dens=''
   namelist/SYSPARAMS/MuChem,FixMuChem,Filling,Beta,Nk1,Nk2,Nk3,&
      file_ham,gauge,Output_Dens,file_dens
   !......................................
   logical  :: relaxation_dynamics=.false.
   integer  :: Nt,output_step=1
   real(dp) :: Tmax,tstart=0.0_dp,tau_relax=1.0e10_dp
   namelist/TIMEPARAMS/Nt,Tmax,output_step,tstart,relaxation_dynamics,tau_relax
   !......................................
   logical            :: ApplyField=.false.
   character(len=255) :: file_field
   namelist/FIELDPARAMS/ApplyField,file_field
   !......................................
   ! -- internal variables --
   logical  :: file_ok,ReadDens
   integer  :: Nsteps,tstp,ppos,tstp_max,step
   real(dp) :: dt,pulse_tmin,pulse_tmax
   real(dp),allocatable :: Etot(:),Ekin(:),BandOcc(:,:),Jcurr(:,:),Dip(:,:)
   real(dp),allocatable :: Jpara(:,:),Jdia(:,:),JHk(:,:),Jpol(:,:),Jintra(:,:)
   real(dp),allocatable :: Jcurr_kpt(:,:,:)
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
      write(output_unit,'(A)') '+------------------------------------------------------+'
      write(output_unit,'(A)') '|       Lattice dynamics in 3D: free evolution         |'
      write(output_unit,'(A)') '+------------------------------------------------------+'
      write(output_unit,*)

      call date_and_time(date=date,time=time)
      write(output_unit,'(a)',advance='no') '  Calculation started on '
      write(output_unit,fmt_time) date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
      write(output_unit,*)

      call Timer_Act
      call Timer_SetName('total',1)
      call Timer_Run(N=1)

      !$OMP PARALLEL PRIVATE(nthreads,threadid)
      threadid = omp_get_thread_num()
      if(threadid == 0) then
         nthreads = omp_get_num_threads()
         write(output_unit,fmt148) 'number of ranks',ntasks
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
      read(unit_inp,nml=TIMEPARAMS) ; rewind(unit_inp)
      read(unit_inp,nml=FIELDPARAMS)
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
   ppos = scan(trim(file_ham),".", BACK= .true.)
   if(trim(file_ham(ppos+1:)) == "h5") then
      call Ham%ReadFromHDF5(file_ham)
   else
      call Ham%ReadFromW90(file_ham)
   end if

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
      write(output_unit,fmt50)
      write(output_unit,'(A)') '         Equilibrium'
      write(output_unit,fmt50)
   end if

   tic = MPI_Wtime()

   call lattsys%Init(Beta,MuChem,ham,Nk1,Nk2,Nk3,gauge)

   call lattsys%SetLaserPulse(external_field)

   if(FixMuChem) then
      call lattsys%SolveEquilibrium()
   else
      call lattsys%SolveEquilibrium(Filling)
   end if
   toc = MPI_Wtime()

   if(on_root) then
      write(output_unit,*)
      write(output_unit,fmt555) "equilibrium",toc-tic
      write(output_unit,*)
      select case(gauge)
      case(dipole_gauge)
         write(output_unit,*) "gauge: length"
      case(dip_emp_gauge)
         write(output_unit,*) "gauge: Peierls substitution"
      case(velocity_gauge)
         write(output_unit,*) "gauge: velocity"
      case(velo_emp_gauge)
         write(output_unit,*) "gauge: empirical velocity"
      case default
         write(error_unit,fmt900) "unrecognized gauge"
      end select
   end if

   ReadDens = len_trim(file_dens) > 0

   if(ReadDens) then
      inquire(file=trim(file_dens),exist=file_ok)
      if(.not.file_ok) then
         if(on_root) write(error_unit,fmt900) 'Input file does not exist: '//trim(file_dens)
         call MPI_Finalize(ierr); stop
      end if
      if(on_root) write(output_unit,fmt_input) 'Density matrix from file: '//trim(file_dens)
      call lattsys%ReadDensM(file_dens)
   end if
!--------------------------------------------------------------------------------------
!                         ++  Time propagation ++
!--------------------------------------------------------------------------------------
   if(on_root) then
      write(output_unit,fmt50)
      write(output_unit,'(A)') '         Time propagation'
      write(output_unit,fmt50)
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

   if(gauge == dipole_gauge .or. gauge == dip_emp_gauge) then
      call lattsys%CalcObservables_dip(0,dt,Ekin(0),Etot(0),Jcurr(:,0),JHk(:,0),Jpol(:,0),&
               Dip(:,0),BandOcc(:,0))
   else
      call lattsys%CalcObservables_velo(0,dt,Ekin(0),Etot(0),Jcurr(:,0),Jpara(:,0),Jdia(:,0),&
         Jintra(:,0),Dip(:,0),BandOcc(:,0))
   end if

   pulse_tmin = min(pulse_x%Tmin,pulse_y%Tmin,pulse_z%Tmin) - tstart
   pulse_tmax = max(pulse_x%Tmax,pulse_y%Tmax,pulse_z%Tmax) - tstart
   if(.not. ApplyField) pulse_tmax = 0.0_dp

   step = 0
   do tstp=0,Nt-1
      if(relaxation_dynamics) then
         call lattsys%Timestep_RelaxTime(tau_relax,tstp,dt)
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

         if(on_root) write(output_unit,fmt145) "tstp",tstp+1

      end if

   end do

   toc = MPI_Wtime()
   if(on_root) then
      write(output_unit,fmt50)
      write(output_unit,*)
      write(output_unit,fmt555) 'propagation',toc-tic
   end if
!--------------------------------------------------------------------------------------
   if(on_root .and. PrintToFile) then
      call SaveObservables(FlOutPref)
   end if

   if(Output_Dens) call lattsys%SaveDensM(trim(FlOutPref)//'_densm.h5')
!--------------------------------------------------------------------------------------
   if(on_root) then
      write(output_unit,*)
      call Timer_Stop(N=1)
      call Timer_DTShow(N=1)
      write(output_unit,'(A)') '+------------------------------------------------------+'
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
   subroutine SaveObservables(prefix)
      use Mhdf5_utils
      character(len=*),intent(in) :: prefix
      integer(HID_t) :: file_id
      character(len=255) :: Flname
      integer  :: tstp
      real(dp) :: AF(3)
      real(dp),allocatable :: ts(:)
      real(dp),allocatable :: EF(:,:)

      Flname = trim(prefix)//'_observables.h5'
      call hdf_open_file(file_id, trim(Flname), STATUS='NEW')
      call hdf_write_attribute(file_id,'','method', 0)
      if(ApplyField) then
         call hdf_write_attribute(file_id,'','applyfield', 1)
      else
         call hdf_write_attribute(file_id,'','applyfield', 0)
      end if
      call hdf_write_attribute(file_id,'','nbnd', lattsys%nbnd)
      call hdf_write_attribute(file_id,'','Nt', Nt)

      allocate(ts(0:Nsteps))
      forall(tstp=0:Nsteps) ts(tstp) = dt * tstp * output_step
      call hdf_write_dataset(file_id,'time',ts)

      call hdf_write_dataset(file_id,'etot',Etot)
      call hdf_write_dataset(file_id,'ekin',Ekin)
      call hdf_write_dataset(file_id,'occ',BandOcc)
      call hdf_write_dataset(file_id,'current',Jcurr)
      call hdf_write_dataset(file_id,'dipole',Dip)

      if(gauge == dipole_gauge .or. gauge == dip_emp_gauge) then
         call hdf_write_dataset(file_id,'current_hk',JHk)
         call hdf_write_dataset(file_id,'current_pol',Jpol)
      else
         call hdf_write_dataset(file_id,'current_para',Jpara)
         call hdf_write_dataset(file_id,'current_dia',Jdia)
         call hdf_write_dataset(file_id,'current_intra',Jintra)
      end if

      if(ApplyField) then
         allocate(EF(3,0:Nsteps))
         do tstp=0,Nsteps
            call external_field(ts(tstp),AF,EF(:,tstp))
         end do
         call hdf_write_dataset(file_id,'efield',EF)
         deallocate(EF)
      end if

      call hdf_close_file(file_id)

   end subroutine SaveObservables
!--------------------------------------------------------------------------------------

!======================================================================================
end program wann_evol_mpi

program wann_evol
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use omp_lib
   use Mdebug
   use Mdef,only: dp,zero
   use Mtime,only: Timer_Act, Timer_Tic, Timer_Toc
   use Mutils,only: print_title, print_header, get_file_ext, check_file_ext
   use Mlaserpulse,only: LaserPulse_spline_t
   use Mham_w90,only: wann90_tb_t
   use Mwann_evol,only: wann_evol_t
   use Mio_hamiltonian,only: ReadHamiltonian
   use Mio_obs,only: SaveTDObs
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! -- constants --
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
   logical            :: FixMuChem=.true.
   integer            :: gauge
   integer            :: Nk1,Nk2,Nk3
   real(dp)           :: MuChem,Beta,Filling
   character(len=255) :: file_ham
   namelist/SYSPARAMS/MuChem,FixMuChem,Filling,Beta,Nk1,Nk2,Nk3,&
      file_ham,gauge
   !......................................
   logical  :: relaxation_dynamics=.false.
   integer  :: Nt,output_step=1
   real(dp) :: Tmax,tstart=0.0_dp,tau_relax=1.0e10_dp
   namelist/TIMEPARAMS/Nt,Tmax,output_step,tstart,relaxation_dynamics,tau_relax
   !......................................
   character(len=255) :: file_field=""
   namelist/FIELDPARAMS/file_field
   !......................................
   ! -- internal variables --
   logical  :: file_ok
   logical  :: ApplyField=.false.
   integer  :: Nsteps,tstp,tstp_max,step
   real(dp) :: dt,pulse_tmin,pulse_tmax
   real(dp),allocatable :: Etot(:),Ekin(:),BandOcc(:,:),Jcurr(:,:),Dip(:,:)
   real(dp),allocatable :: Jpara(:,:),Jdia(:,:),JHk(:,:),Jpol(:,:),Jintra(:,:)
   real(dp),allocatable :: Jcurr_kpt(:,:,:)
   type(LaserPulse_spline_t) :: pulse_x,pulse_y,pulse_z
   type(wann90_tb_t)         :: Ham
   type(wann_evol_t)         :: lattsys
   ! -- parallelization --
   integer  :: nthreads,threadid
!--------------------------------------------------------------------------------------
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
!--------------------------------------------------------------------------------------
!                               ++  Read input ++
!--------------------------------------------------------------------------------------
   call Timer_Tic('Reading input', 2)

   Narg=command_argument_count()
   if(Narg>=1) then
      call get_command_argument(1,FlIn)
      open(newunit=unit_inp,file=trim(FlIn),status='OLD',action='READ')
      read(unit_inp,nml=SYSPARAMS) ; rewind(unit_inp)
      read(unit_inp,nml=TIMEPARAMS) ; rewind(unit_inp)
      read(unit_inp,nml=FIELDPARAMS)
      close(unit_inp)
   else
      write(output_unit,fmt900) 'Please provide a namelist input file.'
      stop
   end if

   if(Narg>=2) then
      call get_command_argument(2,FlOutPref)
      PrintToFile=.true.
   end if
   if(.not.PrintToFile) then
      write(output_unit,fmt_info) 'No output prefix given. No output will be produced.'
   end if

   inquire(file=trim(file_ham),exist=file_ok)
   if(.not.file_ok) then
      write(error_unit,fmt900) 'Input file does not exist: '//trim(file_ham)
      stop
   end if
   write(output_unit,fmt_input) 'Hamiltionian from file: '//trim(file_ham)
   call ReadHamiltonian(file_ham,Ham)

   ApplyField = len_trim(file_field) > 0
   if(ApplyField) then
      inquire(file=trim(file_field),exist=file_ok)
      if(.not.file_ok) then
         write(error_unit,fmt900) 'Input file does not exist: '//trim(file_field)
         stop
      end if
      write(output_unit,fmt_input) 'External field from file: '//trim(file_field)
      call pulse_x%Load_ElectricField(file_field,usecol=2,CalcAfield=.true.)
      call pulse_y%Load_ElectricField(file_field,usecol=3,CalcAfield=.true.)
      call pulse_z%Load_ElectricField(file_field,usecol=4,CalcAfield=.true.)
   end if

   write(output_unit,*)
   call Timer_Toc(N=2)
!--------------------------------------------------------------------------------------
!                               ++  Equilibrium ++
!--------------------------------------------------------------------------------------
   call print_header(output_unit,"Equilibrium","*")
   call Timer_Tic('equilibrium', 2)

   call lattsys%Init(Beta,MuChem,ham,Nk1,Nk2,Nk3,gauge)
   call lattsys%SetLaserPulse(external_field)

   if(FixMuChem) then
      call lattsys%SolveEquilibrium()
   else
      call lattsys%SolveEquilibrium(Filling)
   end if

   write(output_unit,*)
   call Timer_Toc(N=2)
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
!--------------------------------------------------------------------------------------
!                         ++  Time propagation ++
!--------------------------------------------------------------------------------------
   call print_header(output_unit,"Time propagation","*")
   call Timer_Tic('propagation', 2)

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

         write(output_unit,fmt145) "tstp",tstp+1

      end if

   end do

   write(output_unit,fmt50)
   call Timer_Toc(N=2)
!--------------------------------------------------------------------------------------
   if(PrintToFile) then
      select case(gauge)
      case(velocity_gauge, velo_emp_gauge)
         call SaveTDObs(FlOutPref,Nsteps,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
            Jpara,Jdia,Jintra)
      case(dipole_gauge, dip_emp_gauge)
         call SaveTDObs(FlOutPref,Nsteps,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
            JHk,Jpol)         
      end select
   end if
!--------------------------------------------------------------------------------------
   write(output_unit,*)
   call Timer_Toc(N=1)
   write(output_unit,fmt72)
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
end program wann_evol

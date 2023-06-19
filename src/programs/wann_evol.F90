program wann_evol
!! Computes the time-dependent dynamics of an electron system described by a 
!! Wannier Hamiltonian upon laser excitation and computes various observables.
!!
!! ## Description ##
!!
!! Computes the time evolution of the density matrix \(\rho(\mathbf{k},t)\) in the
!! presence of an external electric field \(\mathbf{E}(t)\) (dipole approximation).
!! The evolution can treated as fully coherent or as dissipative. 
!! The initial state \(\rho_\mathrm{eq}(\mathbf{k},t)\) is computed in thermal equilibrium.
!! 
!! In the coherent case the evolution is computed by the unitary time-stepping 
!! $$\rho(\mathbf{k},t_{n+1}) = U(\mathbf{k},t_{n+1},t_n) \rho(\mathbf{k},t_{n}) 
!! U^\dagger(\mathbf{k},t_{n+1},t_n), $$ 
!! where \(U(\mathbf{k},t_{n+1},t_n)\) is the time-evolution operator for time step 
!! \(t_{n+1} = t_n + \Delta t\). We use the 4th order commutator-free method from
!! [J. Comp. Phys. 230, 5930 (2011)](http://www.sciencedirect.com/science/article/pii/S0021999111002300).
!!
!! In the dissipative case we solve 
!! $$ \frac{d}{dt}\rho(\mathbf{k},t) = -i [H(\mathbf{k},t),\rho(\mathbf{k},t) ] + D[\rho(\mathbf{k},t)] $$
!! as an ODE using the 5th order Runge-Kutta-Fehlberg method. The dissipation term is given by
!! $$D[\rho(\mathbf{k},t)] = \frac{\rho(\mathbf{k},t) - \rho_\mathrm{eq}(\mathbf{k})}{T_1} 
!! + \left(\frac{1}{T_1} - \frac{1}{T_2}\right) \rho_\mathrm{off}(\mathbf{k},t) , $$ 
!! where \(\rho_\mathrm{off}(\mathbf{k},t)\) is a purely off-diagonal matrix in the basis of the
!! instantaneous Hamiltonian \(H(\mathbf{k},t)\) (Houston basis).
!!
!! The time-dependent Hamiltonian \(H(\mathbf{k},t)\) is constructed in velocity gauge (VG) or dipole gauge (DG).
!! In case of VG, both the Hamiltonian and the density matrix are presented in the band basis, while
!! the basis of Wannier orbitals is used in the DG. Check 
!! [Phys. Rev. B 103, 1155409 (2021)](https://link.aps.org/doi/10.1103/PhysRevB.103.155409) 
!! for more details.
!!
!! ## How to run ##
!! Run the `wann_evol` program by
!! ```
!!    ./exe/wann_evol.x input_file output_prefix
!! ```
!! If no output file is provided as the second argument, the program will run, but
!! no output will be produced.
!!
!! @Note OpenMP parallization is used for some internal routines. Set `OMP_NUM_THREADS` 
!!       to control the number of threads used.
!!
!! ## Input variables ##
!! Input variables are read from `input_file` in Fortran name list format. The following
!! name list tags and variables are read:
!! #### HAMILTONIAN ####
!! * `MuChem`: the chemical potential (a.u.)
!! * `FixMuChem`: if `.true.`, the input chemical potential will be used. Otherwise the 
!!    chemical potential will be recalculated to match the given filling.
!! * `Filling`: Filling of the bands, corresponding to the total number of electrons per unit cell 
!!   (per spin without SOC).
!! * `Beta`: inverse temperature (a.u.)
!! * `file_ham`: File with the Wannier Hamiltonian. This is the `_tb.dat` file obtained from Wannier90,
!!               or an hdf5 file procuded by [[wann_prune]] or [[wann_soc]].
!! * `lm_gauge`: Velocity gauge (`lm_gauge=0`), dipole gauge (`lm_gauge=1`), empirical velocity gauge (`lm_gauge=2`), 
!!    Peierls substitution (`lm_gauge=3`)

!!
!! #### TIMEPARAMS ####
!! * `Nt`: The number of time steps.
!! * `Tmax`: The propagation time (a.u.)
!! * `output_step`: Calculate observables every `output_step` time steps.
!! * `relaxation_dynamics`: If set to `.true`, the equation of motion for the density matrix with
!!                          phenomenological damping `T1_relax` and decoherence `T2_relax` will be solved.
!! * `T1_relax`: diagonal relaxation towards the instantaneous equilibrium density matrix.
!! * `T2_relax`: decay time of off-diagonal elements of the density matrix.
!! * `file_field`: Data file with three-dimensional electric field. The following format is expected: 
!!    time (1st column), Ex (2nd column), Ey (3rd column), Ez (4th column). 
!!    The time and field strength is expected in atomic units.
!!
!! !! #### OUTPUT ####
!! * `spin_current`: Triggers the calculation of spin currents for systems with SOC. The Hamiltonian
!!                   supports \(n_b = 2 n_\mathrm{orb}\) bands, where \(n_\mathrm{orb}\) is the
!!                   number of Wannier functions without spin. 
!! * `Output_Occ_KPTS`: If `.true.`, the momentum-dependent band occupation will be written to file.
!!                      This can produce large output file size; hdf5 support is recommended.
!!
!! #### KPOINTS ####
!! Specification of the k-points, which can be along a path, a Monkhorst-Pack grid, or a user-supplied. 
!! See [[Read_Kpoints]] for details.
!!
!! ## Output ##
!! If compiled with HDF5 support, after running the code the hdf5 file `output_prefix_observables.h5` will be produced. 
!! The observables can be plotted by using the script `python_utils/plot_obs.h5`.
!! Without HDF5 support there will be several output files `output_prefix_etot.txt`, `out/output_prefix_curr.txt` 
!! etc. that can directly be plotted or read by numpy's `loadtxt`.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use omp_lib
   use Mdebug
   use scitools_def,only: dp,zero
   use scitools_time,only: Timer_Act, Timer_Tic, Timer_Toc
   use scitools_utils,only: print_title, print_header, get_file_ext, check_file_ext,&
      stop_error
   use scitools_laserpulse,only: Laserpulse_3D_t
   use wan_hamiltonian,only: wann90_tb_t
   use wan_latt_kpts,only: Read_Kpoints,kpoints_t
   use io_params,only: HamiltonianParams_t, TimeParams_t
   use io_hamiltonian,only: ReadHamiltonian
   use io_obs,only: SaveTDObs, SaveTDOccupation, SaveSpinCurrent
   use io_density,only: ReadDensityMatrix
   use Mwann_evol,only: wann_evol_t
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! -- constants --
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   character(len=*),parameter :: fmt_input='(" Input: [",a,"]")'
   integer,parameter :: velocity_gauge=0,dipole_gauge=1,velo_emp_gauge=2,dip_emp_gauge=3
   integer,parameter :: prop_unitary=0, prop_rk4=1, prop_rk5=2, prop_hybrid=3
   ! -- for reading i/o --
   logical :: PrintToFile = .false.
   integer :: Narg,unit_inp
   character(len=255) :: FlIn,FlOutPref
   ! -- input variables --
   type(HamiltonianParams_t) :: par_ham
   type(TimeParams_t) :: par_time
   logical :: spin_current=.false.,Output_Occ_KPTS=.false.,save_dens=.false.
   namelist/OUTPUT/spin_current,Output_Occ_KPTS,save_dens
   ! -- internal variables --
   logical  :: file_ok
   logical  :: ApplyField=.false.
   integer  :: Nsteps,tstp,tstp_max,step
   real(dp) :: dt,tstart,tstop,pulse_tmin,pulse_tmax
   real(dp),allocatable,dimension(:)     :: Etot,Ekin
   real(dp),allocatable,dimension(:,:)   :: kpts,BandOcc,Jcurr,Dip
   real(dp),allocatable,dimension(:,:)   :: Jpara,Jdia,JHk,Jpol,Jintra
   real(dp),allocatable,dimension(:,:,:) :: Occk,Jspin
   complex(dp),allocatable,dimension(:,:,:) :: Rhok
   type(Laserpulse_3D_t)     :: pulse
   type(kpoints_t)           :: kp
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
      call par_ham%ReadFromFile(FlIn)
      call par_time%ReadFromFile(FlIn)
      open(newunit=unit_inp,file=trim(FlIn),STATUS='OLD',ACTION='READ')
      read(unit_inp,nml=OUTPUT); rewind(unit_inp)
      call stop_error('Please provide a namelist input file.')
   end if

   if(Narg>=2) then
      call get_command_argument(2,FlOutPref)
      PrintToFile=.true.
   end if
   if(.not.PrintToFile) then
      write(output_unit,fmt_info) 'No output prefix given. No output will be produced.'
   end if


   inquire(file=trim(par_ham%file_ham),exist=file_ok)
   if(.not.file_ok) then
      call stop_error('Input file does not exist: '//trim(par_ham%file_ham))
   end if
   write(output_unit,fmt_input) 'Hamiltionian from file: '//trim(par_ham%file_ham)
   call ReadHamiltonian(par_ham%file_ham,Ham)

   tstart = 0.0_dp
   if(par_time%restart_evolution) then
      inquire(file=trim(par_time%file_dens),exist=file_ok)
      if(.not.file_ok) then
         call stop_error('Input file does not exist: '//trim(par_time%file_dens))
      end if
      write(output_unit,fmt_input) 'Density matrix from file: '//trim(par_time%file_dens)
      call ReadDensityMatrix(par_time%file_dens,Rhok,tstart)  
   end if 

   call Read_Kpoints(FlIn,kp,print_info=.true.)

   ApplyField = len_trim(par_time%file_field) > 0
   if(ApplyField) then
      inquire(file=trim(par_time%file_field),exist=file_ok)
      if(.not.file_ok) then
         call stop_error('Input file does not exist: '//trim(par_time%file_field))
      end if
      write(output_unit,fmt_input) 'External field from file: '//trim(par_time%file_field)
      call pulse%Load_ElectricField(par_time%file_field)
   end if

   write(output_unit,*)
   call Timer_Toc(N=2)
!--------------------------------------------------------------------------------------
!                               ++  Equilibrium ++
!--------------------------------------------------------------------------------------

   if(par_time%restart_evolution) then
      call lattsys%Init(par_ham%Beta,par_ham%MuChem,ham,kp,par_ham%lm_gauge,&
         T1=par_time%T1_relax,T2=par_time%T2_relax,propagator=par_time%propagator,&
         Rhok_start=Rhok,tstart=tstart)
      deallocate(Rhok)
   else
      call lattsys%Init(par_ham%Beta,par_ham%MuChem,ham,kp,par_ham%lm_gauge,&
         T1=par_time%T1_relax,T2=par_time%T2_relax,propagator=par_time%propagator)
   end if

   call lattsys%SetLaserPulse(external_field)

   call print_header(output_unit,"Equilibrium","*")
   call Timer_Tic('equilibrium', 2)

   if(par_ham%FixMuChem) then
      call lattsys%SolveEquilibrium()
      write(output_unit,fmt149) "number of electrons", lattsys%nelec
   else
      call lattsys%SolveEquilibrium(par_ham%Filling)
      write(output_unit,fmt149) "chemical potential", lattsys%MuChem
   end if

   write(output_unit,*)
   call Timer_Toc(N=2)
   write(output_unit,*)

   write(output_unit,fmt_info,advance='no') "gauge: "
   select case(par_ham%lm_gauge)
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
   if(spin_current) write(output_unit,fmt_info) "spin currents will be computed"
!--------------------------------------------------------------------------------------
!                         ++  Time propagation ++
!--------------------------------------------------------------------------------------
   call print_header(output_unit,"Time propagation","*")
   call Timer_Tic('propagation', 2)

   Nsteps = ceiling((par_time%Nt+1)/dble(par_time%output_step)) - 1

   dt = par_time%Tmax/dble(par_time%Nt)

   allocate(Ekin(0:Nsteps),Etot(0:Nsteps),Jcurr(3,0:Nsteps),Dip(3,0:Nsteps))
   allocate(BandOcc(Ham%num_wann,0:Nsteps))
   if(par_ham%lm_gauge == dipole_gauge .or. par_ham%lm_gauge == dip_emp_gauge) then
      allocate(JHk(3,0:Nsteps),Jpol(3,0:Nsteps))
   else
      allocate(Jpara(3,0:Nsteps),Jdia(3,0:Nsteps),Jintra(3,0:Nsteps))
   end if

   if(spin_current) then
      allocate(Jspin(3,3,0:Nsteps))
   end if

   if(Output_Occ_KPTS) then
      allocate(Occk(lattsys%nbnd,lattsys%Nk,0:Nsteps))
   end if

   if(par_ham%lm_gauge == dipole_gauge .or. par_ham%lm_gauge == dip_emp_gauge) then
      call lattsys%CalcObservables_dip(0,dt,Ekin(0),Etot(0),Jcurr(:,0),JHk(:,0),Jpol(:,0),&
               Dip(:,0),BandOcc(:,0))
      if(spin_current) call lattsys%CalcSpinCurrent_dip(tstp,dt,Jspin(:,:,0))
   else
      call lattsys%CalcObservables_velo(0,dt,Ekin(0),Etot(0),Jcurr(:,0),Jpara(:,0),Jdia(:,0),&
         Jintra(:,0),Dip(:,0),BandOcc(:,0))
      if(spin_current) call lattsys%CalcSpinCurrent_velo(tstp,dt,Jspin(:,:,0))
   end if

   if(Output_Occ_KPTS) then
      call lattsys%GetOccupationKPTS(Occk(:,:,0))
   end if

   pulse_tmin = pulse%Tmin
   pulse_tmax = pulse%Tmax
   if(.not. ApplyField) pulse_tmax = 0.0_dp

   step = 0
   do tstp=0,par_time%Nt-1
      if(par_time%propagator == prop_unitary) then
         call lattsys%Timestep(tstp,dt,field_tmax=pulse_tmax)
      else
         call lattsys%Timestep_RelaxTime(tstp,dt)
      end if

      if(mod(tstp+1, par_time%output_step) == 0) then
         step = step + 1

         if(par_ham%lm_gauge == dipole_gauge .or. par_ham%lm_gauge == dip_emp_gauge) then
            call lattsys%CalcObservables_dip(tstp+1,dt,Ekin(step),Etot(step),Jcurr(:,step),JHk(:,step),&
                  Jpol(:,step),Dip(:,step),BandOcc(:,step))
            if(spin_current) call lattsys%CalcSpinCurrent_dip(tstp,dt,Jspin(:,:,step))
         else
            call lattsys%CalcObservables_velo(tstp+1,dt,Ekin(step),Etot(step),Jcurr(:,step),Jpara(:,step),&
                  Jdia(:,step),Jintra(:,step),Dip(:,step),BandOcc(:,step))
            if(spin_current) call lattsys%CalcSpinCurrent_velo(tstp,dt,Jspin(:,:,step))
         end if

         if(Output_Occ_KPTS) then
            call lattsys%GetOccupationKPTS(Occk(:,:,step))
         end if

         write(output_unit,fmt145) "tstp",tstp+1

      end if

   end do

   write(output_unit,fmt50)
   call Timer_Toc(N=2)
!--------------------------------------------------------------------------------------
   if(PrintToFile) then
      select case(par_ham%lm_gauge)
      case(velocity_gauge, velo_emp_gauge)
         call SaveTDObs(FlOutPref,tstart,Nsteps,par_time%output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
            Jpara,Jdia,Jintra)
      case(dipole_gauge, dip_emp_gauge)
         call SaveTDObs(FlOutPref,tstart,Nsteps,par_time%output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
            JHk,Jpol)         
      end select

      if(spin_current) then
         call SaveSpinCurrent(FlOutPref,tstart,Nsteps,par_time%output_step,dt,Jspin)
      end if

      if(Output_Occ_KPTS) then
         call SaveTDOccupation(FlOutPref,tstart,Nsteps,par_time%output_step,dt,Occk)
      end if

      if(save_dens) then
         tstop = tstart + dt * par_time%Nt
         call lattsys%SaveDensityMatrix(FlOutPref,tstop)
         ! call WriteDensityMatrix(FlOutPref,Rhok,tstop)
      end if
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
      if(ApplyField) call pulse%GetField(t,AF,EF)

   end subroutine external_field
!--------------------------------------------------------------------------------------

!======================================================================================
end program wann_evol

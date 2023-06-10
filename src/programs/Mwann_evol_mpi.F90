module Mwann_evol_mpi
!======================================================================================
   use mpi
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,one,zero,nfermi
   use scitools_utils,only: stop_error
   use scitools_linalg,only: get_large_size,util_matmul,util_rotate,util_rotate_cc
   use scitools_array1d_dist,only: dist_array1d_t,GetDisplSize1D
   use io_density,only: WriteDensityMatrix
   use wan_utils,only: Batch_Diagonalize_t
   use wan_latt_kpts,only: kpoints_t
   use wan_hamiltonian,only: wann90_tb_t
   use wan_equilibrium,only: GetChemicalPotential_mpi, Wann_GenRhok_eq
   use wan_rungekutta,only: Init_RungeKutta
   use wan_dynamics
#ifdef WITHFFTW
   use wan_fft_ham_mpi,only: wann_fft_t
   use wan_fft_propagation_mpi,only: Wann_FFT_UnitaryTimestep_dip, Wann_FFT_RelaxTimestep_dip
   use wan_fft_observables_mpi,only: Wann_FFT_Observables_dip
#endif
   implicit none
   include '../formats.h'
   include '../units_inc.f90'
!--------------------------------------------------------------------------------------
   ! .. external vector potential ..
   procedure(vecpot_efield_func),pointer :: field => null()
   ! .. parallelization ..
   integer,parameter  :: master=0,from_master=1,from_worker=2
   integer :: ntasks,taskid,ierr
   logical :: on_root
   integer :: status(MPI_STATUS_SIZE)
   type(dist_array1d_t),private :: kdist
!--------------------------------------------------------------------------------------
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   integer,parameter :: velocity_gauge=0,dipole_gauge=1,velo_emp_gauge=2,dip_emp_gauge=3
   integer,parameter :: kp_list=0, kp_path=1, kp_grid=2, kp_fft_grid_2d=3, kp_fft_grid_3d=4
   integer,parameter :: prop_unitary=0, prop_rk4=1, prop_rk5=2, prop_hybrid=3
   logical :: large_size
!--------------------------------------------------------------------------------------
   private
   public :: wann_evol_t
!--------------------------------------------------------------------------------------
   type wann_evol_t
      logical  :: fft_mode=.false.
      logical  :: free_evol=.false., spin_current=.false., restart=.false.
      integer  :: gauge
      integer  :: propagator=prop_unitary
      integer  :: Nk,Nk_loc,nbnd
      real(dp) :: Beta,MuChem,nelec
      real(dp) :: tstart
      real(dp) :: T1,T2
      real(dp),allocatable :: kcoord_loc(:,:)
      complex(dp),allocatable,dimension(:,:,:)  :: Hk,Udt,wan_rot
      complex(dp),allocatable,dimension(:,:,:,:) :: grad_Hk,Dk,velok
      ! .. Hamiltonian class ..
      type(wann90_tb_t) :: Ham
#ifdef WITHFFTW
      type(wann_fft_t)  :: ham_fft
#endif
      ! .. density matrix ..
      complex(dp),allocatable,dimension(:,:,:)   :: Rhok
   contains
      procedure,public  :: Init
      procedure,public  :: SetLaserpulse
      procedure,public  :: SolveEquilibrium
      procedure,public  :: Timestep_RelaxTime
      procedure,public  :: Timestep
      procedure,public  :: CalcObservables_velo
      procedure,public  :: CalcObservables_dip
      procedure,public  :: CalcSpinCurrent_velo
      procedure,public  :: CalcSpinCurrent_dip    
      procedure,public  :: GetOccupationKPTS
      procedure,public  :: SaveDensityMatrix
   end type wann_evol_t
!--------------------------------------------------------------------------------------
   abstract interface
      subroutine vecpot_efield_func(t,A,E)
         import :: dp
         implicit none
         real(dp),intent(in) :: t
         real(dp),intent(out) :: A(3),E(3)
      end subroutine vecpot_efield_func
   end interface
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,Beta,MuChem,ham,kp,gauge,T1,T2,spin_current,propagator,Rhok_start,&
      tstart)
      class(wann_evol_t)           :: me
      real(dp),intent(in)          :: Beta
      real(dp),intent(in)          :: MuChem
      type(wann90_tb_t),intent(in) :: ham
      type(kpoints_t),intent(in)   :: kp
      integer,intent(in)           :: gauge
      real(dp),intent(in),optional :: T1
      real(dp),intent(in),optional :: T2
      logical,intent(in),optional  :: spin_current
      integer,intent(in),optional  :: propagator
      complex(dp),intent(in),optional :: Rhok_start(:,:,:)
      real(dp),intent(in),optional :: tstart
      integer :: ik,ik_glob

      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
      on_root = taskid == master

      me%gauge = gauge

      me%Nk = size(kp%kpts,dim=1)

      call me%Ham%Set(Ham)

      me%nbnd = me%ham%num_wann
      large_size = get_large_size(me%nbnd)

      call kdist%Init(ntasks,taskid,me%Nk)
      me%Nk_loc = kdist%N_loc(taskid)

      allocate(me%kcoord_loc(me%Nk_loc,3))
      do ik=1,me%Nk_loc
         ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
         me%kcoord_loc(ik,1:3) = kp%kpts(ik_glob,1:3)
      end do

      me%Beta = Beta
      me%MuChem = MuChem
      me%T1 = 1.0e10_dp; me%T2 = 1.0e10_dp
      if(present(T1)) me%T1 = T1
      if(present(T2)) me%T2 = T2

      allocate(me%Rhok(me%nbnd,me%nbnd,me%Nk_loc))

      if(present(Rhok_start)) then
         if(any(shape(Rhok_start) /= [me%nbnd,me%nbnd,me%Nk] )) then
            call stop_error("Incompatible restart density matrix",&
               root_flag=on_root)
         end if
         do ik=1,me%Nk_loc
            ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
            me%Rhok(:,:,ik) = Rhok_start(:,:,ik_glob) 
         end do
         me%restart = .true.
      end if

      me%tstart = 0.0_dp
      if(present(tstart)) me%tstart = tstart

      if(present(spin_current)) me%spin_current = spin_current
      if(me%spin_current) then
         if(mod(me%nbnd, 2) /= 0) then
            call stop_error("Number of bands odd. Not compatible with spinor mode.",&
               root_flag=on_root)
         end if
      end if

      if(present(propagator)) me%propagator = propagator

#ifdef WITHFFTW
      if(kp%kpoints_type == kp_fft_grid_2d) then
         if(on_root) write(output_unit,fmt_info) "Fourier transform: 2D FFT"
         call me%ham_fft%InitFromW90(me%ham, [kp%nk1, kp%nk2], kdist)
         me%fft_mode = .true.
      elseif(kp%kpoints_type == kp_fft_grid_3d) then
         if(on_root) write(output_unit,fmt_info) "Fourier transform: 3D FFT"
         call me%ham_fft%InitFromW90(me%ham, [kp%nk1, kp%nk2, kp%nk3], kdist)   
         me%fft_mode = .true.   
      else
         if(on_root) write(output_unit,fmt_info) "Fourier transform: DFT"
      end if
#endif

      if(.not. me%fft_mode) then
         select case(me%propagator)
         case(prop_rk5)
            call Init_RungeKutta(me%nbnd)
         case(prop_rk4)
            call Init_RungeKutta(me%nbnd,Nk=me%Nk_loc)
         end select
      end if

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine SetLaserpulse(me,Laserfield)
      class(wann_evol_t)      :: me
      procedure(vecpot_efield_func) :: Laserfield

      field => Laserfield

   end subroutine SetLaserpulse
!--------------------------------------------------------------------------------------
   subroutine SolveEquilibrium(me,filling)
      use scitools_linalg,only: EigH,TRace
      use scitools_root,only: brent
      real(dp),parameter :: mu_tol=1.0e-8_dp
      class(wann_evol_t)      :: me
      real(dp),intent(in),optional  :: filling
      logical :: calc_mu=.false.
      integer :: ik,j
      real(dp) :: nel_loc
      real(dp) :: kpt(3)
      type(Batch_Diagonalize_t) :: batch_diag
      real(dp),allocatable :: Ek(:,:)
      complex(dp),allocatable :: Hk(:,:,:)

      calc_mu = present(filling)


      allocate(me%wan_rot(me%nbnd,me%nbnd,me%Nk_loc))

      call batch_diag%Init(me%nbnd)

      allocate(Hk(me%nbnd,me%nbnd,me%Nk_loc),Ek(me%nbnd,me%Nk_loc))

      call Wann_GenHk(me%ham,me%Nk_loc,me%kcoord_loc,Hk)
      do ik=1,me%Nk_loc
         call batch_diag%Diagonalize(Hk(:,:,ik),epsk=Ek(:,ik),vectk=me%wan_rot(:,:,ik))
      end do

      if(calc_mu) then
         me%MuChem = GetChemicalPotential_mpi(me%Nk,Ek,me%Beta,filling)
      end if

      deallocate(Ek,Hk)

      call batch_diag%Clean()

      allocate(me%Hk(me%nbnd,me%nbnd,me%Nk_loc))
      if(me%gauge == velocity_gauge .or. me%gauge == velo_emp_gauge) then
         allocate(me%velok(me%nbnd,me%nbnd,3,me%Nk_loc))
      end if

      select case(me%gauge)
      case(velocity_gauge)
         if(.not. me%restart) then
            call Wann_GenRhok_eq(me%ham,me%Nk_loc,me%kcoord_loc,me%MuChem,me%Beta,me%Rhok,&
               band_basis=.true.)
         end if
         do ik=1,me%Nk_loc
            kpt = me%kcoord_loc(ik,:)
            me%Hk(:,:,ik) = me%ham%get_ham_diag(kpt)
            me%velok(:,:,:,ik) = me%ham%get_velocity(kpt)
         end do
      case(velo_emp_gauge)
         if(.not. me%restart) then
            call Wann_GenRhok_eq(me%ham,me%Nk_loc,me%kcoord_loc,me%MuChem,me%Beta,me%Rhok,&
               band_basis=.true.)   
         end if
         do ik=1,me%Nk_loc
            kpt = me%kcoord_loc(ik,:)
            me%Hk(:,:,ik) = me%ham%get_ham_diag(kpt)
            me%velok(:,:,:,ik) = me%ham%get_emp_velocity(kpt)
         end do
      case(dipole_gauge, dip_emp_gauge)
         if(.not. me%restart) then
            call Wann_GenRhok_eq(me%ham,me%Nk_loc,me%kcoord_loc,me%MuChem,me%Beta,me%Rhok,&
               band_basis=.false.)
         end if
      end select

      nel_loc = 0.0_dp
      do ik=1,me%Nk_loc
         do j=1,me%nbnd
            nel_loc = nel_loc + dble(me%Rhok(j,j,ik)) / me%Nk
         end do
      end do

      call MPI_ALLREDUCE(nel_loc,me%nelec,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)

   end subroutine SolveEquilibrium
!--------------------------------------------------------------------------------------
   subroutine Timestep_RelaxTime(me,tstp,dt)
      class(wann_evol_t)        :: me
      integer,intent(in)        :: tstp
      real(dp),intent(in)       :: dt

      select case(me%gauge)
      case(dipole_gauge)
         ! call Wann_timestep_RelaxTime(me%ham,me%Nk,me%kcoord,tstp,dt,field,T1,T2,&
            ! me%Rhok_Eq,me%Rhok)
         if(me%fft_mode) then
#ifdef WITHFFTW
            call Wann_FFT_RelaxTimestep_dip(me%ham,me%ham_fft,kdist,tstp,dt,me%tstart,field,me%T1,me%T2,&
               me%Beta,me%MuChem,me%Rhok,method=me%propagator)
#endif
         else
            call Wann_timestep_RelaxTime_dip(me%ham,me%Nk_loc,me%kcoord_loc,tstp,dt,me%tstart,field,me%T1,me%T2,&
               me%Beta,me%MuChem,me%Rhok,method=me%propagator)
         end if
      case(dip_emp_gauge)
         if(me%fft_mode) then
#ifdef WITHFFTW
            call Wann_FFT_RelaxTimestep_dip(me%ham,me%ham_fft,kdist,tstp,dt,me%tstart,field,me%T1,me%T2,&
               me%Beta,me%MuChem,me%Rhok,Peierls_only=.true.,method=me%propagator)       
#endif  
         else
            call Wann_timestep_RelaxTime_dip(me%ham,me%Nk_loc,me%kcoord_loc,tstp,dt,me%tstart,field,me%T1,me%T2,&
               me%Beta,me%MuChem,me%Rhok,empirical=.true.,method=me%propagator)      
         end if  
      case(velocity_gauge,velo_emp_gauge)
         call Wann_timestep_RelaxTime_velo_calc(me%nbnd,me%Nk_loc,me%Hk,me%velok,tstp,dt,me%tstart,field,&
            me%T1,me%T2,me%Beta,me%MuChem,me%Rhok,method=me%propagator)
      end select

   end subroutine Timestep_RelaxTime
!--------------------------------------------------------------------------------------
   subroutine Timestep(me,tstp,dt,field_Tmax)
      use scitools_linalg,only: Eye,trace
      use scitools_evol,only: GenU_CF2,UnitaryStepFBW 
      integer,parameter :: qc=1
      class(wann_evol_t) :: me
      integer,intent(in)   :: tstp
      real(dp),intent(in)  :: dt
      real(dp),intent(in),optional :: field_Tmax
      real(dp) :: field_Tmax_,err,np
      integer :: ik,idir
      real(dp) :: AF(3),EF(3),kpt(3)
      complex(dp),dimension(me%nbnd,me%nbnd) :: Hk,Rho_old

      field_Tmax_ = HUGE(1.0_dp)
      if(present(field_Tmax)) field_Tmax_ = field_Tmax

      if((tstp * dt < field_Tmax_) .or. tstp == 0) then
         select case(me%gauge)
         case(dipole_gauge)
            if(me%fft_mode) then
#ifdef WITHFFTW
               call Wann_FFT_UnitaryTimestep_dip(me%ham,me%ham_fft,kdist,tstp,dt,me%tstart,field,me%Rhok)
#endif
            else
               call Wann_Rhok_timestep_dip(me%ham,me%Nk_loc,me%kcoord_loc,tstp,dt,me%tstart,&
                  field,me%Rhok)
            end if
         case(dip_emp_gauge)
            if(me%fft_mode) then
#ifdef WITHFFTW
               call Wann_FFT_UnitaryTimestep_dip(me%ham,me%ham_fft,kdist,tstp,dt,me%tstart,field,me%Rhok,&
                  Peierls_only=.true.)
#endif
            else
               call Wann_Rhok_timestep_dip(me%ham,me%Nk,me%kcoord_loc,tstp,dt,me%tstart,&
                  field,me%Rhok,Peierls_only=.true.)
            end if
         case(velocity_gauge,velo_emp_gauge)
            call Wann_Rhok_timestep_velo_calc(me%nbnd,me%Nk_loc,me%Hk,me%velok,tstp,dt,me%tstart,&
               field,me%Rhok)
         end select

      else
         if(.not. me%free_evol) then
            ! compute time-evolution operator

            if(on_root) then
               write(output_unit,*) "---> Switching to free evolution"
            end if

            allocate(me%Udt(me%nbnd,me%nbnd,me%Nk_loc))
            allocate(me%Dk(me%nbnd,me%nbnd,3,me%Nk_loc)); me%Dk = zero
            allocate(me%grad_Hk(me%nbnd,me%nbnd,3,me%Nk_loc))

            call field(field_Tmax_,AF,EF)
            do ik=1,me%Nk_loc
               kpt = me%kcoord_loc(ik,:)
               select case(me%gauge)
               case(dipole_gauge)
                  me%Hk(:,:,ik) = Wann_GetHk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
                  me%Dk(:,:,:,ik) = Wann_GetDk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
                  me%grad_Hk(:,:,:,ik) = Wann_GetGradHk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
                  call GenU_CF2(dt,me%Hk(:,:,ik),me%Udt(:,:,ik))
               case(dip_emp_gauge)
                  me%Hk(:,:,ik) = Wann_GetHk_dip(me%ham,AF,EF,kpt,reducedA=.false.,Peierls_only=.true.)
                  me%Dk(:,:,:,ik) = Wann_GetDk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
                  me%grad_Hk(:,:,:,ik) = Wann_GetGradHk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
                  call GenU_CF2(dt,me%Hk(:,:,ik),me%Udt(:,:,ik))
               case(velocity_gauge, velo_emp_gauge)
                  Hk = me%Hk(:,:,ik)
                  do idir=1,3
                     Hk = Hk - qc * AF(idir) * me%velok(:,:,idir,ik)
                  end do
                  me%Dk(:,:,:,ik) = me%ham%get_dipole(kpt,band_basis=.true.) 
                  call GenU_CF2(dt,Hk,me%Udt(:,:,ik))
               end select
            end do

            me%free_evol = .true.
         end if

         if(me%free_evol) then
            do ik=1,me%Nk_loc
               Rho_old = me%Rhok(:,:,ik)
               call UnitaryStepFBW(me%nbnd,me%Udt(:,:,ik),Rho_old,me%Rhok(:,:,ik),large_size=large_size)
            end do
         end if
      end if

   end subroutine Timestep
!--------------------------------------------------------------------------------------
   subroutine CalcObservables_velo(me,tstp,dt,Ekin,Etot,Jcurr,Jpara,Jdia,Jintra,Dip,BandOcc)
      class(wann_evol_t),intent(in) :: me
      integer,intent(in)   :: tstp
      real(dp),intent(in)  :: dt
      real(dp),intent(out) :: Ekin
      real(dp),intent(out) :: Etot
      real(dp),intent(out) :: Jcurr(3)
      real(dp),intent(out) :: Jpara(3)
      real(dp),intent(out) :: Jdia(3)
      real(dp),intent(out) :: Jintra(3)
      real(dp),intent(out) :: Dip(3)
      real(dp),intent(out) :: BandOcc(me%nbnd)
      integer :: ik
      real(dp) :: ck,Ekin_loc,Etot_loc,Jp_loc(3),Jd_loc(3),Ji_loc(3)
      real(dp) :: dip_loc(3)
      real(dp) :: AF(3),EF(3)
      real(dp) :: occ(me%nbnd)
      complex(dp),dimension(me%nbnd,me%nbnd) :: rhok_bnd,rot_cc

      AF = 0.0_dp; EF = 0.0_dp
      call field(tstp*dt + me%tstart,AF,EF)

      ck = me%Nk_loc / dble(me%Nk)
      Occ = 0.0_dp

      Ekin_loc = ck * Wann_KineticEn_calc(me%nbnd,me%Nk_loc,me%Hk,me%Rhok)
      Etot_loc = ck * Wann_TotalEn_velo_calc(me%nbnd,me%Nk_loc,me%Hk,me%velok,AF,me%Rhok)
      Jp_loc = ck * Wann_Current_para_velo_calc(me%nbnd,me%Nk_loc,me%velok,me%Rhok) 
      Jd_loc = ck * Wann_Current_dia_velo(me%Nk_loc,AF,me%Rhok) 
      Ji_loc = ck * Wann_Current_Intra_velo_calc(me%nbnd,me%Nk_loc,me%velok,me%Rhok) 

      call MPI_ALLREDUCE(Ekin_loc,Ekin,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(Etot_loc,Etot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(Jp_loc,Jpara,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(Jd_loc,Jdia,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(Ji_loc,Jintra,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      Jcurr = Jpara + Jdia

      if(me%free_evol) then
         Dip_loc = ck * Wann_Pol_dip_calc(me%nbnd,me%Nk_loc,me%Dk,me%Rhok)
      else
         Dip_loc = ck * Wann_Pol_velo(me%ham,me%Nk_loc,me%kcoord_loc,me%Rhok,band_basis=.true.) 
      end if
      call MPI_ALLREDUCE(Dip_loc,Dip,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)

      !$OMP PARALLEL DO REDUCTION(+:Occ)
      do ik=1,me%Nk_loc
         Occ = Occ + dble(GetDiag(me%nbnd,me%Rhok(:,:,ik))) / me%Nk
      end do
      call MPI_ALLREDUCE(Occ,BandOcc,me%nbnd,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      
   end subroutine CalcObservables_velo
!--------------------------------------------------------------------------------------
   subroutine CalcSpinCurrent_velo(me,tstp,dt,Jspin)
      class(wann_evol_t),intent(in) :: me
      integer,intent(in)   :: tstp
      real(dp),intent(in)  :: dt
      real(dp),intent(out) :: Jspin(3,3)
      real(dp) :: ck
      real(dp) :: AF(3),EF(3)
      real(dp),dimension(3,3) :: Jpara_loc,Jdia_loc,Jspin_loc

      AF = 0.0_dp; EF = 0.0_dp
      call field(tstp*dt + me%tstart,AF,EF)

      ck = me%Nk_loc / dble(me%Nk)

      Jpara_loc = Wann_SpinCurrent_para_velo(me%nbnd,me%Nk_loc,me%velok,me%wan_rot,me%Rhok)
      Jdia_loc = Wann_SpinCurrent_dia_velo(me%nbnd,me%Nk_loc,AF,me%wan_rot,me%Rhok)
      Jspin_loc = ck * (Jpara_loc + Jdia_loc)

      call MPI_ALLREDUCE(Jspin_loc,Jspin,9,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)

   end subroutine CalcSpinCurrent_velo
!--------------------------------------------------------------------------------------
   subroutine CalcObservables_dip(me,tstp,dt,Ekin,Etot,Jcurr,Jhk,Jpol,Dip,BandOcc)
      class(wann_evol_t),intent(in) :: me
      integer,intent(in)   :: tstp
      real(dp),intent(in)  :: dt
      real(dp),intent(out) :: Ekin
      real(dp),intent(out) :: Etot
      real(dp),intent(out) :: Jcurr(3)
      real(dp),intent(out) :: JHk(3)
      real(dp),intent(out) :: Jpol(3)
      real(dp),intent(out) :: Dip(3)
      real(dp),intent(out) :: BandOcc(me%nbnd)
      integer :: ik
      real(dp) :: ck,Ekin_loc,Etot_loc,J_loc(3)
      real(dp) :: dip_loc(3),JHk_loc(3)
      real(dp) :: AF(3),EF(3)
      real(dp) :: occ(me%nbnd)
      complex(dp),dimension(me%nbnd,me%nbnd) :: rhok_bnd

      AF = 0.0_dp; EF = 0.0_dp
      call field(tstp*dt + me%tstart,AF,EF)

      ck = me%Nk_loc / dble(me%Nk)
      Occ = 0.0_dp

      if(me%free_evol) then
         Etot_loc = ck * Wann_DTRAB_kpts(me%nbnd,me%Nk_loc,me%Hk,me%Rhok)
         Ekin_loc = ck * Wann_KineticEn(me%ham,me%Nk_loc,me%kcoord_loc,me%Rhok,&
            band_basis=.false.)

         if(me%gauge == dipole_gauge) then
            J_loc = ck * Wann_Current_dip_calc(me%nbnd,me%Nk_loc,me%Hk,me%grad_Hk,me%Dk,me%Rhok)
            JHk_loc = ck * Wann_Current_dip_calc(me%nbnd,me%Nk_loc,me%Hk,me%grad_Hk,me%Dk,me%Rhok,&
               dipole_current=.false.)
         else
            J_loc = ck * Wann_Current_dip_calc(me%nbnd,me%Nk_loc,me%Hk,me%grad_Hk,me%Dk,me%Rhok,&
               dipole_current=.false.)
            JHk_loc = J_loc
         end if

         Dip_loc = ck * Wann_Pol_dip_calc(me%nbnd,me%Nk_loc,me%Dk,me%Rhok)
      elseif(me%fft_mode) then

      print*, taskid, "step", tstp, ck

#ifdef WITHFFTW
         if(me%gauge == dipole_gauge) then
            call Wann_FFT_Observables_dip(me%ham,me%ham_fft,kdist,AF,EF,me%Rhok,Ekin,Etot,Jcurr,Jhk,Dip,&
               dipole_current=.true.)
         else
            call Wann_FFT_Observables_dip(me%ham,me%ham_fft,kdist,AF,EF,me%Rhok,Ekin,Etot,Jcurr,Jhk,Dip,&
               dipole_current=.false.)            
         end if


#endif
      else
         Etot_loc = ck * Wann_TotalEn_dip(me%ham,me%Nk_loc,me%kcoord_loc,AF,EF,me%Rhok)
         Ekin_loc = ck * Wann_KineticEn(me%ham,me%Nk_loc,me%kcoord_loc,me%Rhok,&
            band_basis=.false.)

         if(me%gauge == dipole_gauge) then
            J_loc = ck * Wann_Current_dip(me%ham,me%Nk_loc,me%kcoord_loc,AF,EF,me%Rhok)
            JHk_loc = ck * Wann_Current_dip(me%ham,me%Nk_loc,me%kcoord_loc,AF,EF,me%Rhok,&
               dipole_current=.false.)
         else
            J_loc = ck * Wann_Current_dip(me%ham,me%Nk_loc,me%kcoord_loc,AF,EF,me%Rhok,&
               dipole_current=.false.)
            JHk_loc = J_loc
         end if

         Dip_loc = ck * Wann_Pol_dip(me%ham,me%Nk_loc,me%kcoord_loc,AF,me%Rhok)
      end if

      !$OMP PARALLEL DO REDUCTION(+:Occ) PRIVATE(rhok_bnd)
      do ik=1,me%Nk_loc
         rhok_bnd = util_rotate(me%nbnd,me%wan_rot(:,:,ik),me%Rhok(:,:,ik),large_size=large_size)
         Occ = Occ + dble(GetDiag(me%nbnd,rhok_bnd)) / me%Nk
      end do

      if(.not. me%fft_mode) then      
         call MPI_ALLREDUCE(JHk_loc,Jhk,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(Ekin_loc,Ekin,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(Etot_loc,Etot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(J_loc,Jcurr,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(Dip_loc,Dip,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      end if
      call MPI_ALLREDUCE(Occ,BandOcc,me%nbnd,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)

      if(me%gauge == dipole_gauge) then
         Jpol = Jcurr - JHk 
      else
         Jpol = 0.0_dp
      end if

   end subroutine CalcObservables_dip
!--------------------------------------------------------------------------------------
   subroutine CalcSpinCurrent_dip(me,tstp,dt,Jspin)
      class(wann_evol_t),intent(in) :: me
      integer,intent(in)   :: tstp
      real(dp),intent(in)  :: dt
      real(dp),intent(out) :: Jspin(3,3)
      real(dp) :: AF(3),EF(3)
      real(dp),dimension(3,3) :: Jspin_loc

      AF = 0.0_dp; EF = 0.0_dp
      call field(tstp*dt + me%tstart,AF,EF)

      if(me%free_evol) then
         Jspin_loc = Wann_SpinCurrent_dip_calc(me%nbnd,me%Nk_loc,me%Hk,me%grad_Hk,me%Dk,me%Rhok,&
            dipole_current=(me%gauge == dipole_gauge))
      else
         Jspin_loc = Wann_SpinCurrent_dip(me%ham,me%Nk_loc,me%kcoord_loc,AF,EF,me%Rhok,&
            dipole_current=(me%gauge == dipole_gauge))
      end if

      Jspin_loc = me%Nk_loc / dble(me%Nk) * Jspin_loc

      call MPI_ALLREDUCE(Jspin_loc,Jspin,9,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)      

   end subroutine CalcSpinCurrent_dip
!--------------------------------------------------------------------------------------
   subroutine GetOccupationKPTS(me,Occk)
      class(wann_evol_t),intent(in) :: me
      real(dp),intent(inout)        :: Occk(:,:)
      integer :: ik
      integer :: nsize
      integer,allocatable :: displ(:),size_loc(:)
      real(dp),allocatable :: Occk_loc(:,:)
      complex(dp),allocatable :: Rhok_bnd(:,:)

      call assert_shape(Occk, [me%nbnd,me%Nk], "GetOccupationKPTS", "Occk")

      allocate(Occk_loc(me%nbnd,me%Nk_loc))
      select case(me%gauge)
         case(dipole_gauge, dip_emp_gauge)
            allocate(Rhok_bnd(me%nbnd,me%nbnd))
            do ik=1,me%Nk_loc
               rhok_bnd = util_rotate(me%nbnd,me%wan_rot(:,:,ik),me%Rhok(:,:,ik),large_size=large_size)
               Occk_loc(:,ik) = dble(GetDiag(me%nbnd,rhok_bnd))
            end do
            deallocate(Rhok_bnd)
         case(velocity_gauge,velo_emp_gauge)
            do ik=1,me%Nk_loc
               Occk_loc(:,ik) = dble(GetDiag(me%nbnd,me%Rhok(:,:,ik)))
            end do
      end select

      allocate(displ(0:ntasks-1),size_loc(0:ntasks-1))
      call GetDisplSize1D(kdist%N_loc,me%nbnd,displ,nsize,size_loc)

      call MPI_Allgatherv(Occk_loc,nsize,MPI_DOUBLE_PRECISION,Occk,size_loc,displ,&
         MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)

      deallocate(displ,size_loc)
      deallocate(Occk_loc)

   end subroutine GetOccupationKPTS
!--------------------------------------------------------------------------------------
   subroutine SaveDensityMatrix(me,prefix,tstop)
      class(wann_evol_t),intent(in) :: me
      character(len=*),intent(in)   :: prefix
      real(dp),intent(in)           :: tstop
      integer :: nsize
      integer,allocatable :: displ(:),size_loc(:)
      complex(dp),allocatable :: Rhok(:,:,:)

      if(on_root .or. debug_mode) allocate(Rhok(me%nbnd,me%nbnd,me%Nk))

      allocate(displ(0:ntasks-1),size_loc(0:ntasks-1))
      call GetDisplSize1D(kdist%N_loc,me%nbnd*me%nbnd,displ,nsize,size_loc)

      call MPI_Gatherv(me%Rhok,nsize,MPI_DOUBLE_COMPLEX,Rhok,size_loc,displ,&
         MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)

      deallocate(displ,size_loc)

      if(on_root) then
         call WriteDensityMatrix(prefix,Rhok,tstop)
      end if

      if(allocated(Rhok)) deallocate(Rhok)

   end subroutine SaveDensityMatrix
!--------------------------------------------------------------------------------------
   function GetDiag(ndim,A) result(Adiag)
      integer,intent(in)     :: ndim
      complex(dp),intent(in) :: A(:,:)
      complex(dp)            :: Adiag(ndim)
      integer :: j

      do j=1,ndim
         Adiag(j) = A(j,j)
      end do

   end function GetDiag
!--------------------------------------------------------------------------------------

!======================================================================================
end module Mwann_evol_mpi

module Mwann_evol
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use omp_lib
   use Mdebug
   use scitools_def,only: dp,iu,one,zero,nfermi
   use scitools_utils,only: stop_error
   use scitools_linalg,only: get_large_size,util_matmul,util_rotate,util_rotate_cc
   use wan_latt_kpts,only: kpoints_t
   use wan_utils,only: Batch_Diagonalize_t
   use wan_hamiltonian,only: wann90_tb_t
   use wan_equilibrium,only: GetChemicalPotential, Wann_GenRhok_eq
   use wan_rungekutta,only: Init_RungeKutta
   use wan_dynamics
#ifdef WITHFFTW
   use wan_fft_ham,only: wann_fft_t
   use wan_fft_propagation,only: Wann_FFT_UnitaryTimestep_dip, Wann_FFT_RelaxTimestep_dip
#endif
   implicit none
   include '../formats.h'
   include '../units_inc.f90'
!--------------------------------------------------------------------------------------
   ! .. external vector potential ..
   procedure(vecpot_efield_func),pointer :: field => null()
!--------------------------------------------------------------------------------------
   logical :: large_size
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   integer,parameter :: velocity_gauge=0,dipole_gauge=1,velo_emp_gauge=2,dip_emp_gauge=3
   integer,parameter :: kp_list=0, kp_path=1, kp_grid=2, kp_fft_grid_2d=3, kp_fft_grid_3d=4
   integer,parameter :: prop_unitary=0, prop_rk4=1, prop_rk5=2, prop_hybrid=3
!--------------------------------------------------------------------------------------
   private
   public :: wann_evol_t
!--------------------------------------------------------------------------------------
   type wann_evol_t
      logical  :: free_evol=.false., spin_current=.false., fft_mode=.false.
      integer  :: gauge
      integer  :: propagator=prop_unitary
      integer  :: Nk,nbnd
      real(dp) :: Beta,MuChem,nelec
      real(dp) :: T1,T2
      complex(dp),allocatable,dimension(:,:,:)  :: Hk,Udt,wan_rot
      complex(dp),allocatable,dimension(:,:,:,:) :: grad_Hk,Dk,velok
      ! .. kpoints ..
      real(dp),pointer,dimension(:,:) :: kcoord
      ! .. Hamiltonian class ..
      type(wann90_tb_t) :: Ham
#ifdef WITHFFTW
      type(wann_fft_t)  :: ham_fft
#endif
      ! .. density matrix ..
      complex(dp),allocatable,dimension(:,:,:)   :: Rhok,Rhok_eq
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
   subroutine Init(me,Beta,MuChem,ham,kp,gauge,T1,T2,spin_current,propagator)
      class(wann_evol_t)                :: me
      real(dp),intent(in)               :: Beta
      real(dp),intent(in)               :: MuChem
      type(wann90_tb_t),intent(in)      :: ham
      type(kpoints_t),target,intent(in) :: kp
      integer,intent(in)                :: gauge
      real(dp),intent(in),optional      :: T1
      real(dp),intent(in),optional      :: T2
      logical,intent(in),optional       :: spin_current
      integer,intent(in),optional       :: propagator

      me%gauge = gauge
      me%Nk = kp%nk
      me%kcoord => kp%kpts

      call me%Ham%Set(Ham)

      me%nbnd = me%ham%num_wann
      me%Beta = Beta
      me%MuChem = MuChem
      me%T1 = 1.0e10_dp; me%T2 = 1.0e10_dp
      if(present(T1)) me%T1 = T1
      if(present(T2)) me%T2 = T2

      allocate(me%Rhok(me%nbnd,me%nbnd,me%Nk))

      large_size = get_large_size(me%nbnd)

      if(present(spin_current)) me%spin_current = spin_current
      if(me%spin_current) then
         if(mod(me%nbnd, 2) /= 0) then
            call stop_error("Number of bands odd. Not compatible with spinor mode.")
         end if
      end if

      if(present(propagator)) me%propagator = propagator

      select case(me%propagator)
      case(prop_rk5)
         call Init_RungeKutta(me%nbnd)
      case(prop_rk4)
         call Init_RungeKutta(me%nbnd,Nk=me%Nk)
      end select

#ifdef WITHFFTW
      if(kp%kpoints_type == kp_fft_grid_2d) then
         write(output_unit,fmt_info) "Fourier transform: 2D FFT"
         call me%ham_fft%InitFromW90(me%ham, [kp%nk1, kp%nk2])
         me%fft_mode = .true.
      elseif(kp%kpoints_type == kp_fft_grid_3d) then
         write(output_unit,fmt_info) "Fourier transform: 3D FFT"
         call me%ham_fft%InitFromW90(me%ham, [kp%nk1, kp%nk2, kp%nk3])   
         me%fft_mode = .true.   
      else
         write(output_unit,fmt_info) "Fourier transform: DFT"
      end if
#endif

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
      real(dp) :: kpt(3)
      real(dp),allocatable :: Ek(:,:)
      type(Batch_Diagonalize_t) :: batch_diag
      integer :: nthreads,tid

      calc_mu = present(filling)

      allocate(Ek(me%nbnd,me%Nk))
      allocate(me%Hk(me%nbnd,me%nbnd,me%Nk))
      allocate(me%wan_rot(me%nbnd,me%nbnd,me%Nk))
      allocate(me%Rhok_eq(me%nbnd,me%nbnd,me%Nk))

      !$OMP PARALLEL PRIVATE(tid) DEFAULT(SHARED)
      tid = omp_get_thread_num()
      if(tid == 0) then 
         nthreads = omp_get_num_threads()
      end if
      !$OMP END PARALLEL      

      call batch_diag%Init(me%nbnd,nthreads=nthreads)      

      if(me%fft_mode) then
         call me%ham_fft%GetHam(me%Hk)
      else
         call Wann_GenHk(me%ham,me%Nk,me%kcoord,me%Hk)
      end if

      !$OMP PARALLEL PRIVATE(tid)
      tid = omp_get_thread_num()
      !$OMP DO
      do ik=1,me%Nk
         call batch_diag%Diagonalize(me%Hk(:,:,ik),epsk=Ek(:,ik),vectk=me%wan_rot(:,:,ik),tid=tid)
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      if(calc_mu) then
         me%MuChem = GetChemicalPotential(me%Nk,Ek,me%Beta,filling)
      end if

      if(me%gauge == velocity_gauge .or. me%gauge == velo_emp_gauge) then
         allocate(me%velok(me%nbnd,me%nbnd,3,me%Nk))
      end if

      select case(me%gauge)
      case(velocity_gauge)
         call Wann_GenRhok_eq(me%ham,me%Nk,me%kcoord,me%MuChem,me%Beta,me%Rhok,&
            band_basis=.true.)
         
         !$OMP PARALLEL PRIVATE(kpt)
         !$OMP DO
         do ik=1,me%Nk
            kpt = me%kcoord(ik,:)
            me%Hk(:,:,ik) = me%ham%get_ham_diag(kpt)
            me%velok(:,:,:,ik) = me%ham%get_velocity(kpt)
         end do
         !$OMP END DO
         !$OMP END PARALLEL
      case(velo_emp_gauge)
         call Wann_GenRhok_eq(me%ham,me%Nk,me%kcoord,me%MuChem,me%Beta,me%Rhok,&
            band_basis=.true.)   
         !$OMP PARALLEL PRIVATE(kpt)
         !$OMP DO
         do ik=1,me%Nk
            kpt = me%kcoord(ik,:)
            me%Hk(:,:,ik) = me%ham%get_ham_diag(kpt)
            me%velok(:,:,:,ik) = me%ham%get_emp_velocity(kpt)
         end do
         !$OMP END DO
         !$OMP END PARALLEL
      case(dipole_gauge, dip_emp_gauge)
         call Wann_GenRhok_eq(me%ham,me%Nk,me%kcoord,me%MuChem,me%Beta,me%Rhok,&
            band_basis=.false.)
      end select

      me%Rhok_eq = me%Rhok
      me%nelec = sum(nfermi(me%Beta,Ek - me%MuChem)) / me%Nk

      deallocate(Ek)

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
            call Wann_FFT_RelaxTimestep_dip(me%ham,me%ham_fft,tstp,dt,field,me%T1,me%T2,&
               me%Beta,me%MuChem,me%Rhok,method=me%propagator)
         else
            call Wann_timestep_RelaxTime_dip(me%ham,me%Nk,me%kcoord,tstp,dt,field,me%T1,me%T2,&
               me%Beta,me%MuChem,me%Rhok,method=me%propagator)
         end if
      case(dip_emp_gauge)
         if(me%fft_mode) then
            call Wann_FFT_RelaxTimestep_dip(me%ham,me%ham_fft,tstp,dt,field,me%T1,me%T2,&
               me%Beta,me%MuChem,me%Rhok,Peierls_only=.true.,method=me%propagator)         
         else
            call Wann_timestep_RelaxTime_dip(me%ham,me%Nk,me%kcoord,tstp,dt,field,me%T1,me%T2,&
               me%Beta,me%MuChem,me%Rhok,empirical=.true.,method=me%propagator)      
         end if  
      case(velocity_gauge,velo_emp_gauge)
         call Wann_timestep_RelaxTime_velo_calc(me%nbnd,me%Nk,me%Hk,me%velok,tstp,dt,field,&
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
               call Wann_FFT_UnitaryTimestep_dip(me%ham,me%ham_fft,tstp,dt,field,me%Rhok)
            else
               call Wann_Rhok_timestep_dip(me%ham,me%Nk,me%kcoord,tstp,dt,&
                  field,me%Rhok)
            end if
         case(dip_emp_gauge)
            if(me%fft_mode) then
               call Wann_FFT_UnitaryTimestep_dip(me%ham,me%ham_fft,tstp,dt,field,me%Rhok,&
                  Peierls_only=.true.)
            else
               call Wann_Rhok_timestep_dip(me%ham,me%Nk,me%kcoord,tstp,dt,&
                  field,me%Rhok,Peierls_only=.true.)
            end if
         case(velocity_gauge,velo_emp_gauge)
            call Wann_Rhok_timestep_velo_calc(me%nbnd,me%Nk,me%Hk,me%velok,tstp,dt,&
               field,me%Rhok)
         end select

      else
         if(.not. me%free_evol) then
            ! compute time-evolution operator

            write(output_unit,*) "---> Switching to free evolution"

            allocate(me%Udt(me%nbnd,me%nbnd,me%Nk))
            allocate(me%Dk(me%nbnd,me%nbnd,3,me%Nk)); me%Dk = zero
            allocate(me%grad_Hk(me%nbnd,me%nbnd,3,me%Nk))

            call field(field_Tmax_,AF,EF)
            !$OMP PARALLEL PRIVATE(kpt,Hk)
            !$OMP DO
            do ik=1,me%Nk
               kpt = me%kcoord(ik,:)
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
            !$OMP END DO
            !$OMP END PARALLEL

            me%free_evol = .true.
         end if

         if(me%free_evol) then
            !$OMP PARALLEL PRIVATE(Rho_old)
            !$OMP DO
            do ik=1,me%Nk
               Rho_old = me%Rhok(:,:,ik)
               call UnitaryStepFBW(me%nbnd,me%Udt(:,:,ik),Rho_old,me%Rhok(:,:,ik),large_size=large_size)
            end do
            !$OMP END DO
            !$OMP END PARALLEL            
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
      real(dp) :: AF(3),EF(3)

      AF = 0.0_dp; EF = 0.0_dp
      call field(tstp*dt,AF,EF)

      Ekin = Wann_KineticEn_calc(me%nbnd,me%Nk,me%Hk,me%Rhok)
      Etot = Wann_TotalEn_velo_calc(me%nbnd,me%Nk,me%Hk,me%velok,AF,me%Rhok)
      Jpara = Wann_Current_para_velo_calc(me%nbnd,me%Nk,me%velok,me%Rhok) 
      Jdia = Wann_Current_dia_velo(me%Nk,AF,me%Rhok) 
      Jintra = Wann_Current_Intra_velo_calc(me%nbnd,me%Nk,me%velok,me%Rhok) 
      Jcurr = Jpara + Jdia

      if(me%free_evol) then
         Dip = Wann_Pol_dip_calc(me%nbnd,me%Nk,me%Dk,me%Rhok)
      else
         Dip = Wann_Pol_velo(me%ham,me%Nk,me%kcoord,me%Rhok,band_basis=.true.) 
      end if
    
      BandOcc = 0.0_dp
      do ik=1,me%Nk
         BandOcc = BandOcc + dble(GetDiag(me%nbnd,me%Rhok(:,:,ik))) / me%Nk
      end do
            
   end subroutine CalcObservables_velo
!--------------------------------------------------------------------------------------
   subroutine CalcSpinCurrent_velo(me,tstp,dt,Jspin)
      class(wann_evol_t),intent(in) :: me
      integer,intent(in)   :: tstp
      real(dp),intent(in)  :: dt
      real(dp),intent(out) :: Jspin(3,3)
      real(dp) :: AF(3),EF(3)
      real(dp) :: Jpara(3,3),Jdia(3,3)

      AF = 0.0_dp; EF = 0.0_dp
      call field(tstp*dt,AF,EF)

      Jpara = Wann_SpinCurrent_para_velo(me%nbnd,me%Nk,me%velok,me%wan_rot,me%Rhok)
      Jdia = Wann_SpinCurrent_dia_velo(me%nbnd,me%Nk,AF,me%wan_rot,me%Rhok)
      Jspin = Jpara + Jdia

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
      real(dp) :: AF(3),EF(3)
      complex(dp),dimension(me%nbnd,me%nbnd) :: rhok_bnd

      AF = 0.0_dp; EF = 0.0_dp
      call field(tstp*dt,AF,EF)

      if(me%free_evol) then
         Etot = Wann_DTRAB_kpts(me%nbnd,me%Nk,me%Hk,me%Rhok)

         if(me%gauge == dipole_gauge) then
            Jcurr = Wann_Current_dip_calc(me%nbnd,me%Nk,me%Hk,me%grad_Hk,me%Dk,me%Rhok)
            JHk = Wann_Current_dip_calc(me%nbnd,me%Nk,me%Hk,me%grad_Hk,me%Dk,me%Rhok,&
               dipole_current=.false.)
         else
            Jcurr = Wann_Current_dip_calc(me%nbnd,me%Nk,me%Hk,me%grad_Hk,me%Dk,me%Rhok,&
               dipole_current=.false.)
            JHk = Jcurr
         end if

         Dip = Wann_Pol_dip_calc(me%nbnd,me%Nk,me%Dk,me%Rhok)
      else
         Etot = Wann_TotalEn_dip(me%ham,me%Nk,me%kcoord,AF,EF,me%Rhok)

         if(me%gauge == dipole_gauge) then
            Jcurr = Wann_Current_dip(me%ham,me%Nk,me%kcoord,AF,EF,me%Rhok)
            JHk = Wann_Current_dip(me%ham,me%Nk,me%kcoord,AF,EF,me%Rhok,&
               dipole_current=.false.)
         else
            Jcurr = Wann_Current_dip(me%ham,me%Nk,me%kcoord,AF,EF,me%Rhok,&
               dipole_current=.false.)
            JHk = Jcurr
         end if

         Dip = Wann_Pol_dip(me%ham,me%Nk,me%kcoord,AF,me%Rhok)
      end if

      Ekin = Wann_KineticEn(me%ham,me%Nk,me%kcoord,me%Rhok,&
         band_basis=.false.)
      
      BandOcc = 0.0_dp
      do ik=1,me%Nk
         rhok_bnd = util_rotate(me%nbnd,me%wan_rot(:,:,ik),me%Rhok(:,:,ik),large_size=large_size)
         BandOcc = BandOcc + dble(GetDiag(me%nbnd,rhok_bnd)) / me%Nk
      end do
      
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

      AF = 0.0_dp; EF = 0.0_dp
      call field(tstp*dt,AF,EF)

      if(me%free_evol) then
         Jspin = Wann_SpinCurrent_dip_calc(me%nbnd,me%Nk,me%Hk,me%grad_Hk,me%Dk,me%Rhok,&
            dipole_current=(me%gauge == dipole_gauge))
      else
         Jspin = Wann_SpinCurrent_dip(me%ham,me%Nk,me%kcoord,AF,EF,me%Rhok,&
            dipole_current=(me%gauge == dipole_gauge))
      end if

   end subroutine CalcSpinCurrent_dip
!--------------------------------------------------------------------------------------
   subroutine GetOccupationKPTS(me,Occk)
      class(wann_evol_t),intent(in) :: me
      real(dp),intent(inout)        :: Occk(:,:)
      integer :: ik
      complex(dp),allocatable :: Rhok_bnd(:,:)

      call assert_shape(Occk, [me%nbnd,me%Nk], "GetOccupationKPTS", "Occk")

      select case(me%gauge)
         case(dipole_gauge, dip_emp_gauge)
            !$OMP PARALLEL PRIVATE(Rhok_bnd)
            allocate(Rhok_bnd(me%nbnd,me%nbnd))
            !$OMP DO
            do ik=1,me%Nk
               rhok_bnd = util_rotate(me%nbnd,me%wan_rot(:,:,ik),me%Rhok(:,:,ik),large_size=large_size)
               Occk(:,ik) = dble(GetDiag(me%nbnd,rhok_bnd))
            end do
            !$OMP END DO
            deallocate(Rhok_bnd)
            !$OMP END PARALLEL
         case(velocity_gauge,velo_emp_gauge)
            !$OMP PARALLEL DO
            do ik=1,me%Nk
               Occk(:,ik) = dble(GetDiag(me%nbnd,me%Rhok(:,:,ik)))
            end do
      end select

   end subroutine GetOccupationKPTS
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
end module Mwann_evol

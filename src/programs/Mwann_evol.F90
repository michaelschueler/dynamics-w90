module Mwann_evol
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,iu,one,zero,nfermi
   use Munits,only: DPi
   use Mlinalg,only: util_rotate
   use Mrungekutta,only: ODE_step_rk5
   use Mlatt,only: latt3d_t
   use Mham_w90,only: wann90_tb_t
   use Mwann_dyn
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! .. external vector potential ..
   procedure(vecpot_efield_func),pointer :: field => null()
!--------------------------------------------------------------------------------------
   logical :: large_size
   integer,parameter :: velocity_gauge=0,dipole_gauge=1,velo_emp_gauge=2,dip_emp_gauge=3
!--------------------------------------------------------------------------------------
   private
   public :: wann_evol_t
!--------------------------------------------------------------------------------------
   type wann_evol_t
      logical  :: free_evol=.false.
      integer  :: gauge
      integer  :: Nk,nbnd
      real(dp) :: Beta,MuChem
      complex(dp),allocatable,dimension(:,:,:)  :: Hk,Udt,wan_rot
      complex(dp),allocatable,dimension(:,:,:,:) :: grad_Hk,Dk,velok
      ! .. lattice class ..
      type(latt3d_t) :: latt
      ! .. Hamiltonian class ..
      type(wann90_tb_t) :: Ham
      ! .. density matrix ..
      complex(dp),allocatable,dimension(:,:,:)   :: Rhok,Rhok_eq
   contains
      procedure,public  :: Init
      procedure,public  :: SetLaserpulse
      procedure,public  :: SolveEquilibrium
      procedure,public  :: SaveDensM
      procedure,public  :: ReadDensM
      procedure,public  :: Timestep_RelaxTime
      procedure,public  :: Timestep
      procedure,public  :: CalcObservables_velo
      procedure,public  :: CalcObservables_dip
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
   subroutine Init(me,Beta,MuChem,ham,Nk1,Nk2,Nk3,gauge)
      class(wann_evol_t)           :: me
      real(dp),intent(in)          :: Beta
      real(dp),intent(in)          :: MuChem
      type(wann90_tb_t),intent(in) :: ham
      integer,intent(in)           :: Nk1,Nk2,Nk3
      integer,intent(in)           :: gauge
      character(len=32) :: sampling
      integer :: ik,ik_glob

      me%gauge = gauge

      call me%latt%Init(Nk1,Nk2,Nk3,ham%recip_lattice)
      me%Nk = me%latt%Nk

      call me%Ham%Set(Ham)

      me%nbnd = me%ham%num_wann
      me%Beta = Beta
      me%MuChem = MuChem

      allocate(me%Rhok(me%nbnd,me%nbnd,me%Nk))

      large_size = me%nbnd >= 8

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine SetLaserpulse(me,Laserfield)
      class(wann_evol_t)      :: me
      procedure(vecpot_efield_func) :: Laserfield

      field => Laserfield

   end subroutine SetLaserpulse
!--------------------------------------------------------------------------------------
   subroutine SaveDensM(me,fname)
      use Mhdf5_utils
      class(wann_evol_t)        :: me
      character(len=*),intent(in) :: fname
      real(dp),allocatable :: rdata(:,:,:)
      integer(HID_t) :: file_id

      call hdf_open_file(file_id, trim(fname), STATUS='NEW')
      call hdf_write_attribute(file_id,'','nk1', me%latt%Nk1)
      call hdf_write_attribute(file_id,'','nk2', me%latt%Nk2)
      call hdf_write_attribute(file_id,'','nk3', me%latt%Nk3)
      call hdf_write_attribute(file_id,'','nbnd', me%nbnd)

      allocate(rdata(me%nbnd,me%nbnd,me%Nk))
      rdata = dble(me%Rhok)
      call hdf_write_dataset(file_id,'real',rdata)

      rdata = aimag(me%Rhok)
      call hdf_write_dataset(file_id,'imag',rdata)

      call hdf_close_file(file_id)
      deallocate(rdata)

   end subroutine SaveDensM
!--------------------------------------------------------------------------------------
   subroutine ReadDensM(me,fname)
      use Mhdf5_utils
      class(wann_evol_t)        :: me
      character(len=*),intent(in) :: fname
      integer :: nk1,nk2,nk3,nbnd
      real(dp),allocatable :: rdata(:,:,:)
      integer(HID_t) :: file_id

      call hdf_open_file(file_id, trim(fname), STATUS='OLD', ACTION='READ')
      call hdf_read_attribute(file_id,'','nk1', nk1)
      call hdf_read_attribute(file_id,'','nk2', nk2)
      call hdf_read_attribute(file_id,'','nk3', nk3)
      call hdf_read_attribute(file_id,'','nbnd', nbnd)

      call assert(nk1 == me%latt%Nk1)
      call assert(nk2 == me%latt%Nk2)
      call assert(nk3 == me%latt%Nk3)
      call assert(nbnd == me%nbnd)

      allocate(rdata(me%nbnd,me%nbnd,me%Nk))
      call hdf_read_dataset(file_id,'real',rdata)
      me%Rhok = rdata
      call hdf_read_dataset(file_id,'imag',rdata)
      me%Rhok = me%Rhok + iu * rdata

      call hdf_close_file(file_id)
      deallocate(rdata)

   end subroutine ReadDensM
!--------------------------------------------------------------------------------------
   subroutine SolveEquilibrium(me,filling)
      use Mlinalg,only: EigHE,TRace
      use Mroot,only: brent
      real(dp),parameter :: mu_tol=1.0e-8_dp
      class(wann_evol_t)      :: me
      real(dp),intent(in),optional  :: filling
      logical :: calc_mu=.false.
      integer :: ik,j
      real(dp) :: npart,npart_target,Emin,Emax
      real(dp) :: kpt(3)
      real(dp),allocatable :: Ek(:,:)
      complex(dp),allocatable :: Hk(:,:,:)

      calc_mu = present(filling)

      allocate(Hk(me%nbnd,me%nbnd,me%Nk),Ek(me%nbnd,me%Nk))
      allocate(me%wan_rot(me%nbnd,me%nbnd,me%Nk))
      allocate(me%Rhok_eq(me%nbnd,me%nbnd,me%Nk))

      call Wann_GenHk(me%ham,me%Nk,me%latt%kcoord,Hk)
      do ik=1,me%Nk
         call EigHE(Hk(:,:,ik),Ek(:,ik),me%wan_rot(:,:,ik))
         Ek(:,ik) = me%ham%get_eig(me%latt%kcoord(ik,:))
      end do

      if(calc_mu) then
         npart_target = filling
         Emin = minval(Ek); Emax = maxval(Ek)
         me%MuChem = brent(part_func,Emin,Emax,mu_tol)
      end if

      allocate(me%Hk(me%nbnd,me%nbnd,me%Nk))
      if(me%gauge == velocity_gauge .or. me%gauge == velo_emp_gauge) then
         allocate(me%velok(me%nbnd,me%nbnd,me%Nk,3))
      end if

      select case(me%gauge)
      case(velocity_gauge)
         call Wann_GenRhok_eq(me%ham,me%Nk,me%latt%kcoord,me%MuChem,me%Beta,me%Rhok,&
            band_basis=.true.)
         
         do ik=1,me%Nk
            kpt = me%latt%kcoord(ik,:)
            me%Hk(:,:,ik) = me%ham%get_ham_diag(kpt)
            me%velok(:,:,ik,1:3) = me%ham%get_velocity(kpt)
         end do
      case(velo_emp_gauge)
         call Wann_GenRhok_eq(me%ham,me%Nk,me%latt%kcoord,me%MuChem,me%Beta,me%Rhok,&
            band_basis=.true.)   
         do ik=1,me%Nk
            kpt = me%latt%kcoord(ik,:)
            me%Hk(:,:,ik) = me%ham%get_ham_diag(kpt)
            me%velok(:,:,ik,1:3) = me%ham%get_emp_velocity(kpt)
         end do
      case(dipole_gauge)
         call Wann_GenRhok_eq(me%ham,me%Nk,me%latt%kcoord,me%MuChem,me%Beta,me%Rhok,&
            band_basis=.false.)
      case(dip_emp_gauge)
         call Wann_GenRhok_eq(me%ham,me%Nk,me%latt%kcoord,me%MuChem,me%Beta,me%Rhok,&
            band_basis=.false.)
      end select

      me%Rhok_eq = me%Rhok

      deallocate(Ek,Hk)
   !.............................................
   contains
   !............................................
      real(dp) function num_part(mu)
         real(dp),intent(in) :: mu
         integer :: ik,al
         real(dp) :: np_loc

         num_part = sum(nfermi(me%Beta,Ek - mu))/me%Nk

      end function num_part
   !............................................
      real(dp) function part_func(mu)
         real(dp),intent(in) :: mu
         real(dp) :: npart_mu

         part_func = num_part(mu) - npart_target

      end function part_func
   !............................................
   end subroutine SolveEquilibrium
!--------------------------------------------------------------------------------------
   subroutine Timestep_RelaxTime(me,tau,tstp,dt)
      integer,parameter :: qc=1
      class(wann_evol_t)      :: me
      real(dp),intent(in)       :: tau
      integer,intent(in)        :: tstp
      real(dp),intent(in)       :: dt
      integer :: ik
      real(dp) :: kpt(3)
      complex(dp),allocatable :: Rhok_eq(:,:),Rhok_old(:,:)

      allocate(Rhok_eq(me%nbnd,me%nbnd),Rhok_old(me%nbnd,me%nbnd))

      do ik=1,me%Nk
         kpt = me%latt%kcoord(ik,:)
         Rhok_eq = me%Rhok_eq(:,:,ik)
         Rhok_old = me%Rhok(:,:,ik)

         select case(me%gauge)
         case(dipole_gauge)
            me%Rhok(:,:,ik) = ODE_step_RK5(me%nbnd,tstp,dt,deriv_dipole_gauge,Rhok_old)
         case(dip_emp_gauge)
            me%Rhok(:,:,ik) = ODE_step_RK5(me%nbnd,tstp,dt,deriv_dip_emp_gauge,Rhok_old)
         case(velocity_gauge,velo_emp_gauge)
            me%Rhok(:,:,ik) = ODE_step_RK5(me%nbnd,tstp,dt,deriv_velocity_gauge,Rhok_old)
         end select

      end do

      deallocate(Rhok_eq,Rhok_old)
      !......................................................
      contains
      !......................................................
      function deriv_dipole_gauge(nst,t,yt) result(dydt)
         integer,intent(in) :: nst
         real(dp),intent(in) :: t 
         complex(dp),intent(in) :: yt(:,:)
         complex(dp) :: dydt(nst,nst)
         complex(dp) :: ht(nst,nst)
         real(dp) :: AF(3),EF(3)

         AF = 0.0_dp; EF = 0.0_dp
         if(associated(field)) call field(t,AF,EF)

         ht = Wann_GetHk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
         dydt = -iu*(matmul(ht,yt) - matmul(yt,ht))
         dydt = dydt - (yt - Rhok_eq) / tau

      end function deriv_dipole_gauge
      !....................................................
      function deriv_dip_emp_gauge(nst,t,yt) result(dydt)
         integer,intent(in) :: nst
         real(dp),intent(in) :: t 
         complex(dp),intent(in) :: yt(:,:)
         complex(dp) :: dydt(nst,nst)
         complex(dp) :: ht(nst,nst)
         real(dp) :: AF(3),EF(3)

         AF = 0.0_dp; EF = 0.0_dp
         if(associated(field)) call field(t,AF,EF)

         ht = Wann_GetHk_dip(me%ham,AF,EF,kpt,reducedA=.false.,&
            Peierls_only=.true.)
         dydt = -iu*(matmul(ht,yt) - matmul(yt,ht))
         dydt = dydt - (yt - Rhok_eq) / tau

      end function deriv_dip_emp_gauge
      !....................................................
      function deriv_velocity_gauge(nst,t,yt) result(dydt)
         integer,intent(in) :: nst
         real(dp),intent(in) :: t 
         complex(dp),intent(in) :: yt(:,:)
         complex(dp) :: dydt(nst,nst)
         complex(dp) :: ht(nst,nst)
         real(dp) :: AF(3),EF(3)
         integer :: idir

         AF = 0.0_dp; EF = 0.0_dp
         if(associated(field)) call field(t,AF,EF)

         ht = me%Hk(:,:,ik)
         do idir=1,3
            ht = ht - qc * AF(idir) * me%velok(:,:,ik,idir)
         end do
         dydt = -iu*(matmul(ht,yt) - matmul(yt,ht))
         dydt = dydt - (yt - Rhok_eq) / tau

      end function deriv_velocity_gauge
      !....................................................

   end subroutine Timestep_RelaxTime
!--------------------------------------------------------------------------------------
   subroutine Timestep(me,tstp,dt,field_Tmax)
      use Mlinalg,only: Eye,trace
      use Mevol,only: GenU_CF2,UnitaryStepFBW 
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
            call Wann_Rhok_timestep_dip(me%ham,me%Nk,me%latt%kcoord,tstp,dt,&
               field,me%Rhok)
         case(dip_emp_gauge)
            call Wann_Rhok_timestep_dip(me%ham,me%Nk,me%latt%kcoord,tstp,dt,&
               field,me%Rhok,Peierls_only=.true.)
         case(velocity_gauge,velo_emp_gauge)
            call Wann_Rhok_timestep_velo_calc(me%nbnd,me%Nk,me%Hk,me%velok,tstp,dt,&
               field,me%Rhok)
         end select

      else
         if(.not. me%free_evol) then
            ! compute time-evolution operator

            write(output_unit,*) "---> Switching to free evolution"

            allocate(me%Udt(me%nbnd,me%nbnd,me%Nk))
            allocate(me%Dk(me%nbnd,me%nbnd,me%Nk,3)); me%Dk = zero
            allocate(me%grad_Hk(me%nbnd,me%nbnd,me%Nk,3))

            call field(field_Tmax_,AF,EF)
            do ik=1,me%Nk
               kpt = me%latt%kcoord(ik,:)
               select case(me%gauge)
               case(dipole_gauge)
                  me%Hk(:,:,ik) = Wann_GetHk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
                  me%Dk(:,:,ik,:) = Wann_GetDk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
                  me%grad_Hk(:,:,ik,:) = Wann_GetGradHk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
                  call GenU_CF2(dt,me%Hk(:,:,ik),me%Udt(:,:,ik))
               case(dip_emp_gauge)
                  me%Hk(:,:,ik) = Wann_GetHk_dip(me%ham,AF,EF,kpt,reducedA=.false.,Peierls_only=.true.)
                  me%Dk(:,:,ik,:) = Wann_GetDk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
                  me%grad_Hk(:,:,ik,:) = Wann_GetGradHk_dip(me%ham,AF,EF,kpt,reducedA=.false.)
                  call GenU_CF2(dt,me%Hk(:,:,ik),me%Udt(:,:,ik))
               case(velocity_gauge, velo_emp_gauge)
                  Hk = me%Hk(:,:,ik)
                  do idir=1,3
                     Hk = Hk - qc * AF(idir) * me%velok(:,:,ik,idir)
                  end do
                  me%Dk(:,:,ik,:) = me%ham%get_dipole(kpt,band_basis=.true.) 
                  call GenU_CF2(dt,Hk,me%Udt(:,:,ik))
               end select
            end do

            me%free_evol = .true.
         end if

         if(me%free_evol) then
            do ik=1,me%Nk
               Rho_old = me%Rhok(:,:,ik)
               call UnitaryStepFBW(me%nbnd,me%Udt(:,:,ik),Rho_old,me%Rhok(:,:,ik))
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
      real(dp) :: AF(3),EF(3)

      AF = 0.0_dp; EF = 0.0_dp
      if(associated(field)) call field(tstp*dt,AF,EF)

      Ekin = Wann_KineticEn_calc(me%nbnd,me%Nk,me%Hk,me%Rhok)
      Etot = Wann_TotalEn_velo_calc(me%nbnd,me%Nk,me%Hk,me%velok,AF,me%Rhok)
      Jpara = Wann_Current_para_velo_calc(me%nbnd,me%Nk,me%velok,me%Rhok) 
      Jdia = Wann_Current_dia_velo(me%Nk,AF,me%Rhok) 
      Jintra = Wann_Current_Intra_velo_calc(me%nbnd,me%Nk,me%velok,me%Rhok) 
      Jcurr = Jpara + Jdia

      if(me%free_evol) then
         Dip = Wann_Pol_dip_calc(me%nbnd,me%Nk,me%Dk,me%Rhok)
      else
         Dip = Wann_Pol_velo(me%ham,me%Nk,me%latt%kcoord,me%Rhok,band_basis=.true.) 
      end if
    
      BandOcc = 0.0_dp
      do ik=1,me%Nk
         BandOcc = BandOcc + dble(GetDiag(me%nbnd,me%Rhok(:,:,ik))) / me%Nk
      end do
            
   end subroutine CalcObservables_velo
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
      if(associated(field)) call field(tstp*dt,AF,EF)

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
         Etot = Wann_TotalEn_dip(me%ham,me%Nk,me%latt%kcoord,AF,EF,me%Rhok)

         if(me%gauge == dipole_gauge) then
            Jcurr = Wann_Current_dip(me%ham,me%Nk,me%latt%kcoord,AF,EF,me%Rhok)
            JHk = Wann_Current_dip(me%ham,me%Nk,me%latt%kcoord,AF,EF,me%Rhok,&
               dipole_current=.false.)
         else
            Jcurr = Wann_Current_dip(me%ham,me%Nk,me%latt%kcoord,AF,EF,me%Rhok,&
               dipole_current=.false.)
            JHk = Jcurr
         end if

         Dip = Wann_Pol_dip(me%ham,me%Nk,me%latt%kcoord,AF,me%Rhok)
      end if

      Ekin = Wann_KineticEn(me%ham,me%Nk,me%latt%kcoord,me%Rhok,&
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

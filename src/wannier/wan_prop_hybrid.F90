module wan_prop_hybrid
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp, zero, one, iu, nfermi
   use scitools_linalg,only: util_matmul,util_rotate,util_rotate_cc,get_large_size
   use scitools_evol,only: GenU_CF2,UnitaryStepFBW
   use wan_hamiltonian,only: wann90_tb_t
   use wan_utils,only: Batch_Diagonalize_t
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: hybrid_propagator_t
!--------------------------------------------------------------------------------------
   type :: hybrid_propagator_t
      logical :: large_size
      integer :: nthreads
      integer :: nbnd
      real(dp) :: Beta,Mu
      real(dp) :: Gmm(3),AF(3),Ared(3),EF(3)
      type(Batch_Diagonalize_t) :: batch_diag
      ! .. temporary work space ..
      complex(dp),allocatable,dimension(:,:,:)   :: Hkt,Udt,Udt2,Dscatt
      complex(dp),allocatable,dimension(:,:,:)   :: Rho_eq,Rho_off,Rho_tmp
      complex(dp),allocatable,dimension(:,:,:,:) :: Dipk,Ckt
   contains
      procedure, public :: Init
      procedure, public :: Clean
      procedure, public :: Prepare_Timestep
      procedure, public :: TimeStep_dip
      procedure, public :: TimeStep_velo
   end type hybrid_propagator_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,nbnd,Beta,Mu,T1,T2,nthreads)
      class(hybrid_propagator_t) :: me
      integer,intent(in) :: nbnd
      real(dp),intent(in) :: Beta
      real(dp),intent(in) :: Mu
      real(dp),intent(in) :: T1
      real(dp),intent(in) :: T2
      integer,intent(in),optional :: nthreads

      me%nbnd = nbnd
      me%large_size = get_large_size(me%nbnd)

      me%nthreads = 1
      if(present(nthreads)) me%nthreads = nthreads

      me%Beta = Beta
      me%Mu = Mu
      me%Gmm(1) = 1.0_dp / T1
      me%Gmm(2) = 1.0_dp / T2
      me%Gmm(3) = me%Gmm(1) - me%Gmm(2)

      call me%batch_diag%Init(me%nbnd,nthreads=me%nthreads)

      allocate(me%Hkt(me%nbnd,me%nbnd,0:me%nthreads-1))
      allocate(me%Udt(me%nbnd,me%nbnd,0:me%nthreads-1))
      allocate(me%Udt2(me%nbnd,me%nbnd,0:me%nthreads-1))
      allocate(me%Dscatt(me%nbnd,me%nbnd,0:me%nthreads-1))
      allocate(me%Ckt(me%nbnd,me%nbnd,3,0:me%nthreads-1))
      allocate(me%Dipk(me%nbnd,me%nbnd,3,0:me%nthreads-1))

      allocate(me%Rho_eq(me%nbnd,me%nbnd,0:me%nthreads-1))
      allocate(me%Rho_off(me%nbnd,me%nbnd,0:me%nthreads-1))
      allocate(me%Rho_tmp(me%nbnd,me%nbnd,0:me%nthreads-1))

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine Clean(me)
      class(hybrid_propagator_t) :: me

      call me%batch_diag%Clean()
      if(allocated(me%Hkt)) deallocate(me%Hkt)
      if(allocated(me%Udt)) deallocate(me%Udt)
      if(allocated(me%Udt2)) deallocate(me%Udt2)
      if(allocated(me%Ckt)) deallocate(me%Ckt)
      if(allocated(me%Rho_eq)) deallocate(me%Rho_eq)
      if(allocated(me%Rho_off)) deallocate(me%Rho_off)
      if(allocated(me%Rho_tmp)) deallocate(me%Rho_tmp)

   end subroutine
!--------------------------------------------------------------------------------------
   subroutine Prepare_Timestep(me,tstp,dt,tstart,field,w90)
      class(hybrid_propagator_t)            :: me
      integer,intent(in)                    :: tstp
      real(dp),intent(in)                   :: dt
      real(dp),intent(in)                   :: tstart !! starting time
      type(wann90_tb_t),intent(in),optional :: w90      
      interface
         subroutine field(t,A,E)
            import :: dp
            implicit none
            real(dp),intent(in) :: t
            real(dp),intent(out) :: A(3),E(3)
         end subroutine field
      end interface
      real(dp) :: tn

      tn  = (tstp + 0.5_dp) * dt + tstart
      call field(tn, me%AF, me%EF)

      if(present(w90)) then
         me%Ared = w90%get_kreduced(me%AF)
      end if

   end subroutine Prepare_Timestep
!--------------------------------------------------------------------------------------
   subroutine TimeStep_dip(me,w90,kpt,tstp,dt,Rhok,empirical,tid)
      class(hybrid_propagator_t)   :: me
      type(wann90_tb_t),intent(in) :: w90
      real(dp),intent(in)          :: kpt(3)      
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      complex(dp),intent(inout)    :: Rhok(:,:)
      logical,intent(in),optional  :: empirical
      integer,intent(in),optional  :: tid
      logical :: empirical_
      integer :: tid_
      integer :: i
      real(dp) :: dt2,dt23,dt6
      real(dp) :: kA(3)

      empirical_ = .false.
      if(present(empirical)) empirical_ = empirical

      tid_ = 0
      if(present(tid)) tid_ = tid

      kA = kpt - me%Ared
      me%Hkt(:,:,tid_) = w90%get_ham(kA)
      if(.not.empirical_) then
         me%Dipk(:,:,:,tid_) = w90%get_dipole(kA)
         me%Hkt(:,:,tid_) = me%Hkt(:,:,tid_) - me%EF(1) * me%Dipk(:,:,1,tid_) &
            - me%EF(2) * me%Dipk(:,:,2,tid_) - me%EF(3) * me%Dipk(:,:,3,tid_) 
      end if

      ! call GenU_CF2(dt, me%Hkt, me%Udt)

      call me%batch_diag%Diagonalize(me%Hkt(:,:,tid_), tid=tid_)

      dt2 = 0.5_dp * dt
      dt23 = 2.0_dp / 3.0_dp * dt
      dt6 = dt / 6.0_dp
      call GenU_Eig(dt,me%batch_diag%eps(:,tid_),me%batch_diag%vect(:,:,tid_),me%large_size,&
         me%Ckt(:,:,1,tid_),me%Udt(:,:,tid_))
      call GenU_Eig(dt2,me%batch_diag%eps(:,tid_),me%batch_diag%vect(:,:,tid_),me%large_size,&
         me%Ckt(:,:,1,tid_),me%Udt2(:,:,tid_))

      me%Rho_tmp(:,:,tid_) = zero
      do i=1,me%nbnd
         me%Rho_tmp(i,i,tid_) = nfermi(me%Beta, me%batch_diag%eps(i,tid_) - me%Mu)
      end do
      me%Rho_eq(:,:,tid_) = util_rotate_cc(me%nbnd,me%batch_diag%vect(:,:,tid_),me%Rho_tmp(:,:,tid_),&
         large_size=me%large_size)

      me%Rho_tmp(:,:,tid_) = util_rotate(me%nbnd,me%batch_diag%vect(:,:,tid_),Rhok,large_size=me%large_size)
      do i=1,me%nbnd
         me%Rho_tmp(i,i,tid_) = zero
      end do
      me%Rho_off(:,:,tid_) = util_rotate_cc(me%nbnd,me%batch_diag%vect(:,:,tid_),me%Rho_tmp(:,:,tid_),&
         large_size=me%large_size)      

      me%Dscatt(:,:,tid_) = -me%Gmm(1) * (Rhok - me%Rho_eq(:,:,tid_)) + me%Gmm(3) * me%Rho_off(:,:,tid_)

      me%Rho_tmp(:,:,tid_) = Rhok + dt6 * me%Dscatt(:,:,tid_)
      call UnitaryStepFBW(me%nbnd,me%Udt(:,:,tid_),me%Rho_tmp(:,:,tid_),Rhok,large_size=me%large_size)

      me%Ckt(:,:,1,tid_) = dt23 * me%Dscatt(:,:,tid_)
      call UnitaryStepFBW(me%nbnd,me%Udt2(:,:,tid_),me%Ckt(:,:,1,tid_),me%Ckt(:,:,2,tid_),&
         large_size=me%large_size)

      Rhok = Rhok + me%Ckt(:,:,2,tid_) + dt6 * me%Dscatt(:,:,tid_)

   end subroutine TimeStep_dip
!--------------------------------------------------------------------------------------
   subroutine TimeStep_velo(me,ik,tstp,dt,Hk,vk,Rhok,tid)
      class(hybrid_propagator_t)   :: me
      integer,intent(in)           :: ik
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      complex(dp),intent(in)       :: Hk(:,:,:)
      complex(dp),intent(in)       :: vk(:,:,:,:)      
      complex(dp),intent(inout)    :: Rhok(:,:)
      integer,intent(in),optional  :: tid
      integer :: tid_
      integer :: i
      real(dp) :: dt2,dt23,dt6

      tid_ = 0
      if(present(tid)) tid_ = tid

      me%Hkt(:,:,tid_) = Hk(:,:,ik) - me%AF(1) * vk(:,:,1,ik) &
         - me%AF(2) * vk(:,:,2,ik) - me%AF(3) * vk(:,:,3,ik) 
      
      ! call GenU_CF2(dt, me%Hkt, me%Udt)

      call me%batch_diag%Diagonalize(me%Hkt(:,:,tid_), tid=tid_)

      dt2 = 0.5_dp * dt
      dt23 = 2.0_dp / 3.0_dp * dt
      dt6 = dt / 6.0_dp
      call GenU_Eig(dt,me%batch_diag%eps(:,tid_),me%batch_diag%vect(:,:,tid_),me%large_size,&
         me%Ckt(:,:,1,tid_),me%Udt(:,:,tid_))
      call GenU_Eig(dt2,me%batch_diag%eps(:,tid_),me%batch_diag%vect(:,:,tid_),me%large_size,&
         me%Ckt(:,:,1,tid_),me%Udt2(:,:,tid_))

      me%Rho_tmp(:,:,tid_) = zero
      do i=1,me%nbnd
         me%Rho_tmp(i,i,tid_) = nfermi(me%Beta, me%batch_diag%eps(i,tid_) - me%Mu)
      end do
      me%Rho_eq(:,:,tid_) = util_rotate_cc(me%nbnd,me%batch_diag%vect(:,:,tid_),me%Rho_tmp(:,:,tid_),large_size=me%large_size)

      me%Rho_tmp(:,:,tid_) = util_rotate(me%nbnd,me%batch_diag%vect(:,:,tid_),Rhok,large_size=me%large_size)
      do i=1,me%nbnd
         me%Rho_tmp(i,i,tid_) = zero
      end do
      me%Rho_off(:,:,tid_) = util_rotate_cc(me%nbnd,me%batch_diag%vect(:,:,tid_),me%Rho_tmp(:,:,tid_),large_size=me%large_size)      

      me%Dscatt(:,:,tid_) = -me%Gmm(1) * (Rhok - me%Rho_eq(:,:,tid_)) + me%Gmm(3) * me%Rho_off(:,:,tid_)
      
      me%Rho_tmp(:,:,tid_) = Rhok + dt6 * me%Dscatt(:,:,tid_)
      call UnitaryStepFBW(me%nbnd,me%Udt(:,:,tid_),me%Rho_tmp(:,:,tid_),Rhok,large_size=me%large_size)

      me%Ckt(:,:,1,tid_) = dt23 * me%Dscatt(:,:,tid_)
      call UnitaryStepFBW(me%nbnd,me%Udt2(:,:,tid_),me%Ckt(:,:,1,tid_),me%Ckt(:,:,2,tid_),large_size=me%large_size)

      Rhok = Rhok + me%Ckt(:,:,2,tid_) + dt6 * me%Dscatt(:,:,tid_)

   end subroutine TimeStep_velo
!--------------------------------------------------------------------------------------
   subroutine GenU_Eig(dt,Ek,rotk,large_size,work,Udt)
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Ek(:)
      complex(dp),intent(in) :: rotk(:,:)
      logical,intent(in) :: large_size
      complex(dp),intent(inout) :: work(:,:)
      complex(dp),intent(inout) :: Udt(:,:)
      integer :: n,i

      n = size(Ek)

      work = zero
      do i=1,n
         work(i,i) = exp(-iu*dt*Ek(i))
      end do

      Udt = util_rotate_cc(n,rotk,work,large_size=large_size)

   end subroutine GenU_Eig
!--------------------------------------------------------------------------------------



!======================================================================================
end module wan_prop_hybrid

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
      integer :: nbnd
      integer :: comm_order=3
      real(dp) :: Beta,Mu
      real(dp) :: Gmm(3),AF(3),Ared(3),EF(3)
      type(Batch_Diagonalize_t) :: batch_diag
      ! .. temporary work space ..
      complex(dp),allocatable,dimension(:,:)   :: Hkt,Udt,Dscatt
      complex(dp),allocatable,dimension(:,:)   :: Rho_eq,Rho_off,Rho_tmp
      complex(dp),allocatable,dimension(:,:,:) :: Dipk,Ckt
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
   subroutine Init(me,nbnd,Beta,Mu,T1,T2,comm_order)
      class(hybrid_propagator_t) :: me
      integer,intent(in) :: nbnd
      real(dp),intent(in) :: Beta
      real(dp),intent(in) :: Mu
      real(dp),intent(in) :: T1
      real(dp),intent(in) :: T2
      integer,intent(in),optional :: comm_order

      me%nbnd = nbnd
      me%large_size = get_large_size(me%nbnd)

      me%Beta = Beta
      me%Mu = Mu
      me%Gmm(1) = 1.0_dp / T1
      me%Gmm(2) = 1.0_dp / T2
      me%Gmm(3) = me%Gmm(1) - me%Gmm(2)

      me%comm_order = 3
      if(present(comm_order)) me%comm_order = comm_order

      call me%batch_diag%Init(me%nbnd)

      allocate(me%Hkt(me%nbnd,me%nbnd))
      allocate(me%Udt(me%nbnd,me%nbnd))
      allocate(me%Dscatt(me%nbnd,me%nbnd))
      allocate(me%Ckt(me%nbnd,me%nbnd,me%comm_order))
      allocate(me%Dipk(me%nbnd,me%nbnd,3))

      allocate(me%Rho_eq(me%nbnd,me%nbnd))
      allocate(me%Rho_off(me%nbnd,me%nbnd))
      allocate(me%Rho_tmp(me%nbnd,me%nbnd))

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine Clean(me)
      class(hybrid_propagator_t) :: me

      call me%batch_diag%Clean()
      if(allocated(me%Hkt)) deallocate(me%Hkt)
      if(allocated(me%Udt)) deallocate(me%Udt)
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
   subroutine TimeStep_dip(me,w90,kpt,tstp,dt,Rhok,empirical)
      class(hybrid_propagator_t)   :: me
      type(wann90_tb_t),intent(in) :: w90
      real(dp),intent(in)          :: kpt(3)      
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      complex(dp),intent(inout)    :: Rhok(:,:)
      logical,intent(in),optional  :: empirical
      logical :: empirical_
      integer :: i,icomm
      real(dp) :: kA(3)
      complex(dp) :: zh,mzh

      empirical_ = .false.
      if(present(empirical)) empirical_ = empirical

      kA = kpt - me%Ared
      me%Hkt = w90%get_ham(kA)
      if(.not.empirical_) then
         me%Dipk = w90%get_dipole(kA)
         me%Hkt = me%Hkt - me%EF(1) * me%Dipk(:,:,1) - me%EF(2) * me%Dipk(:,:,2) - me%EF(3) * me%Dipk(:,:,3) 
      end if

      call GenU_CF2(dt, me%Hkt, me%Udt)

      call me%batch_diag%Diagonalize(me%Hkt)
      me%Rho_tmp = zero
      do i=1,me%nbnd
         me%Rho_tmp(i,i) = nfermi(me%Beta, me%batch_diag%eps(i,0) - me%Mu)
      end do
      me%Rho_eq = util_rotate_cc(me%nbnd,me%batch_diag%vect(:,:,0),me%Rho_tmp,large_size=me%large_size)

      me%Rho_tmp = util_rotate(me%nbnd,me%batch_diag%vect(:,:,0),Rhok,large_size=me%large_size)
      do i=1,me%nbnd
         me%Rho_tmp(i,i) = zero
      end do
      me%Rho_off = util_rotate_cc(me%nbnd,me%batch_diag%vect(:,:,0),me%Rho_tmp,large_size=me%large_size)      

      me%Dscatt = -me%Gmm(1) * (Rhok - me%Rho_eq) + me%Gmm(3) * me%Rho_off

      me%Ckt = zero
      me%Ckt(:,:,1) = -iu * dt * me%Dscatt
      do icomm=1,me%comm_order - 1
         zh = iu * dt / (icomm + 1.0_dp)
         mzh = -zh
         ! me%Ckt(:,:,icomm+1) = zh * (matmul(me%Hkt, me%Ckt(:,:,icomm)) - matmul(me%Ckt(:,:,icomm), me%Hkt))

         call ZGEMM('N','N',me%nbnd,me%nbnd,me%nbnd,zh,me%Hkt(1,1),me%nbnd,me%Ckt(1,1,icomm),&
            me%nbnd,zero,me%Ckt(1,1,icomm+1),me%nbnd)
         call ZGEMM('N','N',me%nbnd,me%nbnd,me%nbnd,mzh,me%Ckt(1,1,icomm),me%nbnd,me%Hkt(1,1),&
            me%nbnd,one,me%Ckt(1,1,icomm+1),me%nbnd)        

         ! call util_matmul(me%Hkt, me%Ckt(:,:,icomm), me%Ckt(:,:,icomm+1), alpha=zh, &
         !    large_size=me%large_size)
         ! call util_matmul(me%Ckt(:,:,icomm), me%Hkt, me%Ckt(:,:,icomm+1), alpha=-zh, beta=one, &
         !    large_size=me%large_size)
      end do

      me%Rho_tmp = Rhok
      do icomm=1,me%comm_order
         me%Rho_tmp = me%Rho_tmp + iu * me%Ckt(:,:,icomm)
      end do

      call UnitaryStepFBW(me%nbnd,me%Udt,me%Rho_tmp,Rhok,large_size=me%large_size)

   end subroutine TimeStep_dip
!--------------------------------------------------------------------------------------
   subroutine TimeStep_velo(me,ik,tstp,dt,Hk,vk,Rhok)
      class(hybrid_propagator_t)   :: me
      integer,intent(in)           :: ik
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      complex(dp),intent(in)       :: Hk(:,:,:)
      complex(dp),intent(in)       :: vk(:,:,:,:)      
      complex(dp),intent(inout)    :: Rhok(:,:)
      integer :: i,icomm
      complex(dp) :: zh

      me%Hkt = Hk(:,:,ik)
      me%Hkt = me%Hkt - me%AF(1) * vk(:,:,ik,1) - me%AF(2) * vk(:,:,ik,2) - me%AF(3) * vk(:,:,ik,3) 

      call GenU_CF2(dt, me%Hkt, me%Udt)

      call me%batch_diag%Diagonalize(me%Hkt)
      me%Rho_tmp = zero
      do i=1,me%nbnd
         me%Rho_tmp(i,i) = nfermi(me%Beta, me%batch_diag%eps(i,0) - me%Mu)
      end do
      me%Rho_eq = util_rotate_cc(me%nbnd,me%batch_diag%vect(:,:,0),me%Rho_tmp,large_size=me%large_size)

      me%Rho_tmp = util_rotate(me%nbnd,me%batch_diag%vect(:,:,0),Rhok,large_size=me%large_size)
      do i=1,me%nbnd
         me%Rho_tmp(i,i) = zero
      end do
      me%Rho_off = util_rotate_cc(me%nbnd,me%batch_diag%vect(:,:,0),me%Rho_tmp,large_size=me%large_size)      

      me%Dscatt = -me%Gmm(1) * (Rhok - me%Rho_eq) + me%Gmm(3) * me%Rho_off

      me%Ckt = zero
      me%Ckt(:,:,1) = -iu * dt * me%Dscatt
      do icomm=1,me%comm_order - 1
         zh = iu * dt / (me%comm_order + 1.0_dp)
         call util_matmul(me%Hkt, me%Ckt(:,:,icomm), me%Ckt(:,:,icomm+1), alpha=zh, &
            large_size=me%large_size)
         call util_matmul(me%Ckt(:,:,icomm), me%Hkt, me%Ckt(:,:,icomm+1), alpha=-zh, beta=one, &
            large_size=me%large_size)
      end do

      me%Rho_tmp = Rhok
      do icomm=1,me%comm_order
         me%Rho_tmp = me%Rho_tmp + me%Ckt(:,:,icomm)
      end do

      call UnitaryStepFBW(me%nbnd,me%Udt,me%Rho_tmp,Rhok,large_size=me%large_size)

   end subroutine TimeStep_velo
!--------------------------------------------------------------------------------------




!======================================================================================
end module wan_prop_hybrid
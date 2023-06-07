module wan_rungekutta
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use omp_lib
   use Mdebug
   use scitools_def,only: dp, zero, one, iu, nfermi
   use scitools_linalg,only: util_matmul,util_rotate,util_rotate_cc,get_large_size
   use scitools_evol,only: GenU_CF2,GenU_CF4,UnitaryStepFBW
   use wan_hamiltonian,only: wann90_tb_t
   use wan_utils,only: Batch_Diagonalize_t
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: Init_RungeKutta, Clean_RungeKutta
   public :: RK5_TimeStep_Dip_k, RK5_TimeStep_velo_k
   public :: RK4_TimeStep_Dip_k, RK4_TimeStep_velo_k

   integer,parameter :: qc = 1

   logical :: large_size
   integer :: nthreads_
   integer :: nbnd_
   real(dp),allocatable,dimension(:,:)      :: epsk
   complex(dp),allocatable,dimension(:,:,:) :: Hkt, rotk, Rhok_old, Rho_eq, Rho_off, Rho_tmp
   complex(dp),allocatable,dimension(:,:,:,:) :: DRhok_step
   complex(dp),allocatable,dimension(:,:,:,:) :: Hkt_step
   type(Batch_Diagonalize_t) :: batch_diag
!--------------------------------------------------------------------------------------
   real(dp) :: ARK5(5),BRK5(6),CRK5(5,5) 

  data ARK5(1),ARK5(2),ARK5(3),ARK5(4),ARK5(5) / &
       0.25D0,0.375D0,0.92307692307692307692D0,1.0D0,0.5D0 /

  data BRK5(1),BRK5(2),BRK5(3),BRK5(4),BRK5(5),BRK5(6) / &
       0.11851851851851851851D0,&
       0.0D0,&
       0.51898635477582846003D0,&
       0.50613149034201665780D0,&
       -0.18D0,&
       0.03636363636363636363D0 /
       
  data CRK5(1,1),CRK5(1,2),CRK5(1,3),CRK5(1,4),CRK5(1,5) / &
       0.25D0,0.0D0,0.0D0,0.0D0,0.0D0 /
  data CRK5(2,1),CRK5(2,2),CRK5(2,3),CRK5(2,4),CRK5(2,5) / &
       0.09375D0,0.28125D0,0.0D0,0.0D0,0.0D0 /
  data CRK5(3,1),CRK5(3,2),CRK5(3,3),CRK5(3,4),CRK5(3,5) / &
       0.87938097405553026854D0,&
       -3.27719617660446062812D0,&
       3.32089212562585343650D0,&
       0.0D0,0.0D0 /
  data CRK5(4,1),CRK5(4,2),CRK5(4,3),CRK5(4,4),CRK5(4,5) / &
       2.03240740740740740740D0,&
       -8.0D0,&
       7.17348927875243664717D0,&
       -0.20589668615984405458D0,&
       0.0D0 /
  data CRK5(5,1),CRK5(5,2),CRK5(5,3),CRK5(5,4),CRK5(5,5) / &
       -0.29629629629629629629D0,&
       2.0D0,&
       -1.38167641325536062378D0,&
       0.45297270955165692007D0,&
       -0.275D0 /
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init_RungeKutta(nbnd,Nk,nthreads)
      integer,intent(in) :: nbnd
      integer,intent(in),optional :: Nk
      integer,intent(in),optional :: nthreads

      nbnd_ = nbnd
      large_size = get_large_size(nbnd)

      nthreads_ = 1
      if(present(nthreads)) nthreads_ = nthreads

      call batch_diag%Init(nbnd,nthreads=nthreads_)

      allocate(epsk(nbnd,0:nthreads_-1))
      allocate(Hkt(nbnd,nbnd,0:nthreads_-1))
      allocate(rotk(nbnd,nbnd,0:nthreads_-1))
      allocate(Rhok_old(nbnd,nbnd,0:nthreads_-1))
      allocate(Rho_eq(nbnd,nbnd,0:nthreads_-1))
      allocate(Rho_tmp(nbnd,nbnd,0:nthreads_-1))
      allocate(Rho_off(nbnd,nbnd,0:nthreads_-1))
      allocate(DRhok_step(nbnd,nbnd,6,0:nthreads_-1))

      if(present(Nk)) then
         allocate(Hkt_step(nbnd,nbnd,Nk,3))
      end if

   end subroutine Init_RungeKutta
!--------------------------------------------------------------------------------------
   subroutine Clean_RungeKutta()

      call batch_diag%Clean()

      deallocate(epsk)
      deallocate(Hkt)
      deallocate(rotk)
      deallocate(Rhok_old)
      deallocate(Rho_eq)
      deallocate(Rho_tmp)
      deallocate(Rho_off)
      deallocate(DRhok_step)

      if(allocated(Hkt_step)) deallocate(Hkt_step)

   end subroutine Clean_RungeKutta
!--------------------------------------------------------------------------------------
   subroutine RK5_TimeStep_Dip_k(w90,kpt,tstp,dt,field,T1,T2,beta,mu,Rhok,empirical,tid)
      type(wann90_tb_t),intent(in) :: w90
      real(dp),intent(in)          :: kpt(3)
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      real(dp),intent(in)          :: T1
      real(dp),intent(in)          :: T2
      real(dp),intent(in)          :: beta
      real(dp),intent(in)          :: mu
      complex(dp),intent(inout)    :: Rhok(:,:)
      logical,intent(in),optional  :: empirical
      integer,intent(in),optional  :: tid
      interface
         subroutine field(t,A,E)
            import :: dp
            implicit none
            real(dp),intent(in) :: t
            real(dp),intent(out) :: A(3),E(3)
         end subroutine field
      end interface
      logical :: empirical_
      integer :: tid_
      integer :: k,j
      real(dp) :: tk
      real(dp) :: gm(3),AF(3),EF(3)

      empirical_ = .false.
      if(present(empirical)) empirical_ = empirical

      tid_ = 0
      if(present(tid)) tid_ = tid

      gm(1) = 1.0_dp / T1
      gm(2) = 1.0_dp / T2
      gm(3) = gm(1) - gm(2)

      Rhok_old(:,:,tid_) = Rhok

      tk = tstp * dt
      call field(tk, AF, EF)
      DRhok_step(:,:,1,tid_) = deriv_dipole_gauge(nbnd_,w90,kpt,AF,EF,beta,mu,gm,Rhok_old(:,:,tid_),empirical_)
      do k=1,5
         Rhok = Rhok_old(:,:,tid_)
         do j=1,k
            Rhok = Rhok + dt * CRK5(k,j)*DRhok_step(:,:,j,tid_)
         end do

         tk = (tstp + ARK5(k))*dt
         call field(tk, AF, EF)
         DRhok_step(:,:,k+1,tid_) = deriv_dipole_gauge(nbnd_,w90,kpt,AF,EF,beta,mu,gm,Rhok,empirical_)
      end do

      Rhok = Rhok_old(:,:,tid_)
      do j=1,6
         Rhok = Rhok + dt*BRK5(j)*DRhok_step(:,:,j,tid_)
      end do

   end subroutine RK5_TimeStep_Dip_k
!--------------------------------------------------------------------------------------
   subroutine RK4_TimeStep_dip_k(w90,kpt,ik,tstp,dt,field,T1,T2,beta,mu,Rhok,empirical,tid)
      type(wann90_tb_t),intent(in) :: w90
      real(dp),intent(in)          :: kpt(3)
      integer,intent(in)           :: ik
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      real(dp),intent(in)          :: T1
      real(dp),intent(in)          :: T2
      real(dp),intent(in)          :: beta
      real(dp),intent(in)          :: mu
      complex(dp),intent(inout)    :: Rhok(:,:)
      logical,intent(in),optional  :: empirical
      integer,intent(in),optional  :: tid
      interface
         subroutine field(t,A,E)
            import :: dp
            implicit none
            real(dp),intent(in) :: t
            real(dp),intent(out) :: A(3),E(3)
         end subroutine field
      end interface
      logical :: empirical_
      integer :: tid_
      real(dp) :: tk
      real(dp) :: gm(3),AF(3),EF(3)

      empirical_ = .false.
      if(present(empirical)) empirical_ = empirical

      tid_ = 0
      if(present(tid)) tid_ = tid

      gm(1) = 1.0_dp / T1
      gm(2) = 1.0_dp / T2
      gm(3) = gm(1) - gm(2)

      Rhok_old(:,:,tid_) = Rhok

      ! Hk_step(:,:,1) ---> H(t_n)
      ! Hk_step(:,:,2) ---> H(t_n + dt/2)
      ! Hk_step(:,:,3) ---> H(t_n + dt)

      if(tstp == 0) then
         tk = tstp * dt
         call field(tk, AF, EF)     
         Hkt_step(:,:,ik,1) = Wann_GetHk_dip(w90,AF,EF,kpt,reducedA=.false.,&
            Peierls_only=empirical_)            
      end if

      tk =(tstp + 0.5_dp) * dt
      call field(tk, AF, EF)   
      Hkt_step(:,:,ik,2) = Wann_GetHk_dip(w90,AF,EF,kpt,reducedA=.false.,&
         Peierls_only=empirical_)   

      tk =(tstp + 1.0_dp) * dt
      call field(tk, AF, EF)   
      Hkt_step(:,:,ik,3) = Wann_GetHk_dip(w90,AF,EF,kpt,reducedA=.false.,&
            Peierls_only=empirical_)   

      ! Runge-Kutta steps
      DRhok_step(:,:,1,tid_) = GetDeriv(nbnd_,beta,mu,gm,Hkt_step(:,:,ik,1),Rhok_old(:,:,tid_),tid=tid_)
      Rhok = Rhok_old(:,:,tid_) + 0.5_dp * dt * DRhok_step(:,:,1,tid_)
      DRhok_step(:,:,2,tid_) = GetDeriv(nbnd_,beta,mu,gm,Hkt_step(:,:,ik,2),Rhok,tid=tid_)
      Rhok = Rhok_old(:,:,tid_) + 0.5_dp * dt * DRhok_step(:,:,2,tid_)
      DRhok_step(:,:,3,tid_) = GetDeriv(nbnd_,beta,mu,gm,Hkt_step(:,:,ik,2),Rhok,tid=tid_)
      Rhok = Rhok_old(:,:,tid_) + dt * DRhok_step(:,:,3,tid_)
      DRhok_step(:,:,4,tid_) = GetDeriv(nbnd_,beta,mu,gm,Hkt_step(:,:,ik,3),Rhok,tid=tid_)

      Rhok = Rhok_old(:,:,tid_) + dt/6.0_dp * (DRhok_step(:,:,1,tid_) + 2.0_dp * DRhok_step(:,:,2,tid_) &
         + 2.0_dp * DRhok_step(:,:,3,tid_) + DRhok_step(:,:,4,tid_))

      ! update next step
      Hkt_step(:,:,ik,1) = Hkt_step(:,:,ik,3)

   end subroutine RK4_TimeStep_dip_k
!--------------------------------------------------------------------------------------
   subroutine RK5_TimeStep_velo_k(ik,tstp,dt,field,T1,T2,beta,mu,Hk,vk,Rhok,tid)
      integer,intent(in)           :: ik
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      real(dp),intent(in)          :: T1
      real(dp),intent(in)          :: T2
      real(dp),intent(in)          :: beta
      real(dp),intent(in)          :: mu
      complex(dp),intent(in)       :: Hk(:,:,:)
      complex(dp),intent(in)       :: vk(:,:,:,:)
      complex(dp),intent(inout)    :: Rhok(:,:)
      integer,intent(in),optional  :: tid
      interface
         subroutine field(t,A,E)
            import :: dp
            implicit none
            real(dp),intent(in) :: t
            real(dp),intent(out) :: A(3),E(3)
         end subroutine field
      end interface
      integer :: tid_
      integer :: k,j
      real(dp) :: tk
      real(dp) :: gm(3),AF(3),EF(3)

      tid_ = 0
      if(present(tid)) tid_ = tid

      gm(1) = 1.0_dp / T1
      gm(2) = 1.0_dp / T2
      gm(3) = gm(1) - gm(2)

      Rhok_old(:,:,tid_) = Rhok

      tk = tstp * dt
      call field(tk, AF, EF)
      DRhok_step(:,:,1,tid_) = deriv_velo_gauge(nbnd_,ik,Hk,vk,AF,EF,beta,mu,gm,Rhok_old(:,:,tid_))
      do k=1,5
         Rhok = Rhok_old(:,:,tid_)
         do j=1,k
            Rhok = Rhok + dt * CRK5(k,j)*DRhok_step(:,:,j,tid_)
         end do

         tk = (tstp + ARK5(k))*dt
         call field(tk, AF, EF)
         DRhok_step(:,:,k+1,tid_) = deriv_velo_gauge(nbnd_,ik,Hk,vk,AF,EF,beta,mu,gm,Rhok_old(:,:,tid_))
      end do

      Rhok = Rhok_old(:,:,tid_)
      do j=1,6
         Rhok = Rhok + dt*BRK5(j)*DRhok_step(:,:,j,tid_)
      end do

   end subroutine RK5_TimeStep_velo_k
!--------------------------------------------------------------------------------------
   subroutine RK4_TimeStep_velo_k(ik,tstp,dt,field,T1,T2,beta,mu,Hk,vk,Rhok,tid)
     integer,intent(in)           :: ik
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      real(dp),intent(in)          :: T1
      real(dp),intent(in)          :: T2
      real(dp),intent(in)          :: beta
      real(dp),intent(in)          :: mu
      complex(dp),intent(in)       :: Hk(:,:,:)
      complex(dp),intent(in)       :: vk(:,:,:,:)
      complex(dp),intent(inout)    :: Rhok(:,:)
      integer,intent(in),optional  :: tid
      interface
         subroutine field(t,A,E)
            import :: dp
            implicit none
            real(dp),intent(in) :: t
            real(dp),intent(out) :: A(3),E(3)
         end subroutine field
      end interface
      integer :: tid_
      real(dp) :: tk
      real(dp) :: gm(3),AF(3),EF(3)

      tid_ = 0
      if(present(tid)) tid_ = tid

      gm(1) = 1.0_dp / T1
      gm(2) = 1.0_dp / T2
      gm(3) = gm(1) - gm(2)

      Rhok_old(:,:,tid_) = Rhok

      ! Hk_step(:,:,1) ---> H(t_n)
      ! Hk_step(:,:,2) ---> H(t_n + dt/2)
      ! Hk_step(:,:,3) ---> H(t_n + dt)

      if(tstp == 0) then
         tk = tstp * dt
         call field(tk, AF, EF)        
         call UpdateHk_velo(ik, Hk, vk, AF, Hkt_step(:,:,ik,1))
      end if

      tk =(tstp + 0.5_dp) * dt
      call field(tk, AF, EF)   
      call UpdateHk_velo(ik, Hk, vk, AF, Hkt_step(:,:,ik,2))

      tk =(tstp + 1.0_dp) * dt
      call field(tk, AF, EF)   
      call UpdateHk_velo(ik, Hk, vk, AF, Hkt_step(:,:,ik,3))

      ! Runge-Kutta steps
      DRhok_step(:,:,1,tid_) = GetDeriv(nbnd_,beta,mu,gm,Hkt_step(:,:,ik,1),Rhok_old(:,:,tid_),tid=tid_)
      Rhok = Rhok_old(:,:,tid_) + 0.5_dp * dt * DRhok_step(:,:,1,tid_)
      DRhok_step(:,:,2,tid_) = GetDeriv(nbnd_,beta,mu,gm,Hkt_step(:,:,ik,2),Rhok,tid=tid_)
      Rhok = Rhok_old(:,:,tid_) + 0.5_dp * dt * DRhok_step(:,:,2,tid_)
      DRhok_step(:,:,3,tid_) = GetDeriv(nbnd_,beta,mu,gm,Hkt_step(:,:,ik,2),Rhok,tid=tid_)
      Rhok = Rhok_old(:,:,tid_) + dt * DRhok_step(:,:,3,tid_)
      DRhok_step(:,:,4,tid_) = GetDeriv(nbnd_,beta,mu,gm,Hkt_step(:,:,ik,3),Rhok,tid=tid_)

      Rhok = Rhok_old(:,:,tid_) + dt/6.0_dp * (DRhok_step(:,:,1,tid_) + 2.0_dp * DRhok_step(:,:,2,tid_) &
         + 2.0_dp * DRhok_step(:,:,3,tid_) + DRhok_step(:,:,4,tid_))

      ! update next step
      Hkt_step(:,:,ik,1) = Hkt_step(:,:,ik,3)

   end subroutine RK4_TimeStep_velo_k
!--------------------------------------------------------------------------------------
   function GetDeriv(nst,beta,mu,gm,Hk,yt,commutator,tid) result(dydt)
      integer,intent(in)           :: nst
      real(dp),intent(in)          :: beta
      real(dp),intent(in)          :: mu
      real(dp),intent(in)          :: gm(3)
      complex(dp),intent(in)       :: Hk(:,:)
      complex(dp),intent(in)       :: yt(:,:)    
      logical,intent(in),optional  :: commutator
      integer,intent(in),optional  :: tid
      complex(dp) :: dydt(nst,nst) 
      integer :: tid_
      logical :: commutator_
      integer :: i

      commutator_ = .true.
      if(present(commutator)) commutator_ = commutator

      tid_ = 0
      if(present(tid)) tid_ = tid

      call batch_diag%Diagonalize(Hk, epsk=epsk(:,tid_), vectk=rotk(:,:,tid_), tid=tid_)

      Rho_tmp(:,:,tid_) = zero
      do i=1,nst
         Rho_tmp(i,i,tid_) = nfermi(Beta, epsk(i,tid_) - Mu)
      end do
      rho_eq(:,:,tid_) = util_rotate_cc(nst,rotk(:,:,tid_),Rho_tmp(:,:,tid_),large_size=large_size)

      Rho_tmp(:,:,tid_) = util_rotate(nst,rotk(:,:,tid_),yt,large_size=large_size)
      do i=1,nst
         Rho_tmp(i,i,tid_) = zero
      end do
      rho_off(:,:,tid_) = util_rotate_cc(nst,rotk(:,:,tid_),Rho_tmp(:,:,tid_),large_size=large_size)

      dydt = zero
      if(commutator_) then
         call util_matmul(Hk, yt, dydt, alpha=-iu, large_size=large_size)
         call util_matmul(yt, Hk, dydt, alpha=iu, beta=one, large_size=large_size)
      end if
      dydt = dydt - gm(1) * (yt - Rho_eq(:,:,tid_))
      dydt = dydt + gm(3) * rho_off(:,:,tid_)

   end function GetDeriv
!--------------------------------------------------------------------------------------
   function deriv_dipole_gauge(nst,w90,kpt,AF,EF,beta,mu,gm,yt,empirical,commutator,tid) result(dydt)
      integer,intent(in)           :: nst
      type(wann90_tb_t),intent(in) :: w90
      real(dp),intent(in)          :: kpt(3)
      real(dp),intent(in)          :: AF(3),EF(3)
      real(dp),intent(in)          :: beta
      real(dp),intent(in)          :: mu
      real(dp),intent(in)          :: gm(3)
      complex(dp),intent(in)       :: yt(:,:)
      logical,intent(in)           :: empirical
      logical,intent(in),optional  :: commutator
      integer,intent(in),optional  :: tid
      complex(dp) :: dydt(nst,nst)
      logical :: commutator_
      integer :: tid_

      commutator_ = .true.
      if(present(commutator)) commutator_ = commutator

      tid_ = 0
      if(present(tid)) tid_ = tid

      Hkt(:,:,tid_) = Wann_GetHk_dip(w90,AF,EF,kpt,reducedA=.false.,&
         Peierls_only=empirical)

      ! Hkt = w90%get_ham_Peierls_Dipole(kpt,AF,EF,Peierls_only=empirical)

      dydt = GetDeriv(nst,beta,mu,gm,Hkt(:,:,tid_),yt,commutator=commutator_,tid=tid_)

   end function deriv_dipole_gauge
!--------------------------------------------------------------------------------------
   subroutine UpdateHk_velo(ik,Hk,vk,AF,Hkt)
      integer,intent(in)           :: ik
      complex(dp),intent(in)       :: Hk(:,:,:)
      complex(dp),intent(in)       :: vk(:,:,:,:)
      real(dp),intent(in)          :: AF(3)
      complex(dp),intent(inout)    :: Hkt(:,:)
      integer :: idir

      Hkt = Hk(:,:,ik)
      do idir=1,3
         Hkt = Hkt - AF(idir) * vk(:,:,idir,ik)
      end do

   end subroutine UpdateHk_velo
!--------------------------------------------------------------------------------------
   function deriv_velo_gauge(nst,ik,Hk,vk,AF,EF,beta,mu,gm,yt,commutator,tid) result(dydt)
      integer,intent(in)           :: nst
      integer,intent(in)           :: ik
      complex(dp),intent(in)       :: Hk(:,:,:)
      complex(dp),intent(in)       :: vk(:,:,:,:)
      real(dp),intent(in)          :: AF(3),EF(3)
      real(dp),intent(in)          :: beta
      real(dp),intent(in)          :: mu
      real(dp),intent(in)          :: gm(3)
      complex(dp),intent(in)       :: yt(:,:)
      logical,intent(in),optional  :: commutator
      integer,intent(in),optional  :: tid
      complex(dp) :: dydt(nst,nst)
      logical :: commutator_
      integer :: tid_
      integer :: idir

      commutator_ = .true.
      if(present(commutator)) commutator_ = commutator

      tid_ = 0
      if(present(tid)) tid_ = tid

      call UpdateHk_velo(ik, Hk, vk, AF, Hkt(:,:,tid_))

      ! Hkt = Hk(:,:,ik)
      ! do idir=1,3
      !    Hkt = Hkt - AF(idir) * vk(:,:,ik,idir)
      ! end do

      dydt = GetDeriv(nst,beta,mu,gm,Hkt(:,:,tid_),yt,commutator=commutator_,tid=tid_)

   end function deriv_velo_gauge
!--------------------------------------------------------------------------------------
   function Wann_GetHk_dip(w90,Avec,Efield,kpt,reducedA,Peierls_only) result(Hk)
      type(wann90_tb_t),intent(in) :: w90
      real(dp),intent(in)          :: Avec(3),Efield(3)
      real(dp),intent(in)          :: kpt(3)
      logical,intent(in),optional  :: reducedA
      logical,intent(in),optional  :: Peierls_only
      logical :: reducedA_,peierls_
      real(dp) :: Ared(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann) :: Hk
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: Dk

      reducedA_ = .false.
      if(present(reducedA)) reducedA_ = reducedA

      peierls_ = .false.
      if(present(Peierls_only)) peierls_ = Peierls_only

      if(.not.reducedA) then
         Ared = matmul(w90%recip_reduced(1:3,1:3), Avec)
      else
         Ared = Avec
      end if

      Hk = w90%get_ham([kpt(1)-Ared(1),kpt(2)-Ared(2),kpt(3)-Ared(3)]) 
      if(peierls_) return
      Dk = w90%get_dipole([kpt(1)-Ared(1),kpt(2)-Ared(2),kpt(3)-Ared(3)])
      Hk = Hk - qc*(Dk(:,:,1) * EField(1) + Dk(:,:,2) * EField(2) + Dk(:,:,3) * EField(3))

   end function Wann_GetHk_dip
!--------------------------------------------------------------------------------------

!======================================================================================
end module wan_rungekutta
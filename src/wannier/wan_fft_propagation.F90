module wan_fft_propagation
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use omp_lib
   use Mdebug
   use scitools_def,only: dp,iu,one,zero,nfermi
   use scitools_linalg,only: util_matmul,util_rotate,util_rotate_cc,get_large_size
   use scitools_evol,only: GenU_CF2, GenU_CF4, UnitaryStepFBW
   use wan_utils,only: Batch_Diagonalize_t
   use wan_hamiltonian,only: wann90_tb_t
   use wan_fft_ham,only: wann_fft_t
   implicit none
   include '../formats.h'
   include '../units_inc.f90'
!--------------------------------------------------------------------------------------
   integer,parameter :: prop_unitary=0, prop_rk4=1, prop_rk5=2, prop_hybrid=3
   real(dp),parameter :: c1=0.5_dp-sqrt(3.0_dp)/6.0_dp
   real(dp),parameter :: c2=0.5_dp+sqrt(3.0_dp)/6.0_dp

   abstract interface
      subroutine vecpot_efield_func(t,A,E)
         import :: dp
         implicit none
         real(dp),intent(in) :: t
         real(dp),intent(out) :: A(3),E(3)
      end subroutine vecpot_efield_func
   end interface
!--------------------------------------------------------------------------------------
   private
   public :: Wann_FFT_UnitaryTimestep_dip, Wann_FFT_RelaxTimestep_dip
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Wann_FFT_UnitaryTimestep_dip(ham,ham_fft,tstp,dt,field,Rhok,Peierls_only)
      type(wann90_tb_t),intent(in) :: ham
      type(wann_fft_t),intent(in)  :: ham_fft
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      procedure(vecpot_efield_func),pointer :: field
      complex(dp),intent(inout)    :: Rhok(:,:,:)
      logical,intent(in),optional  :: Peierls_only
      logical :: peierls_
      logical :: large_size
      integer :: nbnd,Nk,ik,idir
      real(dp) :: tn
      real(dp),dimension(3,2) :: AF,EF,Ared
      complex(dp),allocatable,dimension(:,:) :: Rho_old,Udt
      complex(dp),allocatable,dimension(:,:,:) :: Hk_1,Hk_2
      complex(dp),allocatable,dimension(:,:,:,:) :: Dk_1,Dk_2
      integer :: tid

      peierls_ = .false.
      if(present(Peierls_only)) peierls_ = Peierls_only

      nbnd = ham_fft%nwan
      Nk = ham_fft%nkpts

      large_size = get_large_size(nbnd)

      call assert_shape(Rhok, [nbnd,nbnd,Nk], "Wann_FFT_UnitaryTimestep_dip", "Rhok")

      tn = (tstp + c1) * dt
      call field(tn, AF(:,1), EF(:,1))
      tn = (tstp + c2) * dt
      call field(tn, AF(:,2), EF(:,2))

      Ared(:,1) = ham%get_kreduced(AF(:,1))
      Ared(:,2) = ham%get_kreduced(AF(:,2))

      allocate(Hk_1(nbnd,nbnd,Nk),Hk_2(nbnd,nbnd,Nk))

      call ham_fft%GetHam_Dressed(Ared(:,1), Hk_1)
      call ham_fft%GetHam_Dressed(Ared(:,2), Hk_2)

      if(.not.peierls_) then
         allocate(Dk_1(nbnd,nbnd,3,Nk),Dk_2(nbnd,nbnd,3,Nk))
         call ham_fft%GetDipole_Dressed(Ared(:,1), Dk_1)
         call ham_fft%GetDipole_Dressed(Ared(:,2), Dk_2)
         do ik=1,Nk
            do idir=1,3
               Hk_1(:,:,ik) = Hk_1(:,:,ik) - EF(idir,1) * Dk_1(:,:,idir,ik)
               Hk_2(:,:,ik) = Hk_2(:,:,ik) - EF(idir,2) * Dk_2(:,:,idir,ik)
            end do
         end do
      end if

      call omp_set_num_threads(ham_fft%nthreads)
      !$OMP PARALLEL PRIVATE(Rho_old,Udt)
      allocate(Rho_old(nbnd,nbnd))
      allocate(Udt(nbnd,nbnd))

      !$OMP DO
      do ik=1,Nk
         Rho_old = Rhok(:,:,ik)
         call GenU_CF4(dt,Hk_1(:,:,ik),Hk_2(:,:,ik),Udt)
         call UnitaryStepFBW(nbnd,Udt,Rho_old,Rhok(:,:,ik),large_size=large_size)
      end do
      !$OMP END DO

      deallocate(Rho_old)
      deallocate(Udt)
      !$OMP END PARALLEL

      deallocate(Hk_1,Hk_2)
      if(allocated(Dk_1)) deallocate(Dk_1)
      if(allocated(Dk_2)) deallocate(Dk_2)


   end subroutine Wann_FFT_UnitaryTimestep_dip
!--------------------------------------------------------------------------------------
   subroutine Wann_FFT_RelaxTimestep_dip(ham,ham_fft,tstp,dt,field,T1,T2,Beta,Mu,Rhok,&
      Peierls_only,method,comm_order)
      type(wann90_tb_t),intent(in) :: ham
      type(wann_fft_t),intent(in)  :: ham_fft
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      procedure(vecpot_efield_func),pointer :: field
      real(dp),intent(in)          :: T1
      real(dp),intent(in)          :: T2
      real(dp),intent(in)          :: Beta
      real(dp),intent(in)          :: Mu
      complex(dp),intent(inout)    :: Rhok(:,:,:)
      logical,intent(in),optional  :: Peierls_only
      integer,intent(in),optional  :: method
      integer,intent(in),optional  :: comm_order
      logical :: peierls_
      integer :: method_, comm_order_

      peierls_ = .false.
      if(present(Peierls_only)) peierls_ = Peierls_only

      method_ = prop_rk4
      if(present(method)) method_ = method

      comm_order_ = 3
      if(present(comm_order)) comm_order_ = comm_order

      select case(method_)
      case(prop_rk4)
         call RelaxTimestep_RK4(ham,ham_fft,tstp,dt,field,T1,T2,Beta,Mu,Rhok,peierls_)
      case(prop_hybrid)
         call RelaxTimestep_Hyb(ham,ham_fft,tstp,dt,field,T1,T2,Beta,Mu,Rhok,peierls_,&
            comm_order_)
      case default
         write(output_unit,fmt700) "Other propagators not implemented. Switching back to Runge-Kutta-4"
         call RelaxTimestep_RK4(ham,ham_fft,tstp,dt,field,T1,T2,Beta,Mu,Rhok,peierls_)
      end select

   end subroutine Wann_FFT_RelaxTimestep_dip
!--------------------------------------------------------------------------------------
   subroutine RelaxTimestep_RK4(ham,ham_fft,tstp,dt,field,T1,T2,Beta,Mu,Rhok,Peierls_only)
      type(wann90_tb_t),intent(in) :: ham
      type(wann_fft_t),intent(in)  :: ham_fft
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      procedure(vecpot_efield_func),pointer :: field
      real(dp),intent(in)          :: T1
      real(dp),intent(in)          :: T2
      real(dp),intent(in)          :: Beta
      real(dp),intent(in)          :: Mu
      complex(dp),intent(inout)    :: Rhok(:,:,:)
      logical,intent(in)           :: Peierls_only
      logical :: large_size
      integer :: nbnd,Nk,ik,idir
      real(dp) :: tn
      real(dp),dimension(3) :: Gmm
      real(dp),dimension(3,3) :: AF,EF,Ared
      complex(dp),allocatable,dimension(:,:) :: Rho_old,Dscatt
      complex(dp),allocatable,dimension(:,:,:) :: DRhok_step
      complex(dp),allocatable,dimension(:,:,:) :: Hk_1,Hk_2,Hk_3
      complex(dp),allocatable,dimension(:,:,:,:) :: Dk_1,Dk_2,Dk_3
      type(Batch_Diagonalize_t) :: batch_diag
      integer :: nthreads,tid

      nbnd = ham_fft%nwan
      Nk = ham_fft%nkpts

      large_size = get_large_size(nbnd)

      call assert_shape(Rhok, [nbnd,nbnd,Nk], "RelaxTimestep_RK4", "Rhok")

      Gmm(1) = 1.0_dp / T1
      Gmm(2) = 1.0_dp / T2
      Gmm(3) = Gmm(1) - Gmm(2)

      tn = tstp * dt
      call field(tn, AF(:,1), EF(:,1))
      tn = (tstp + 0.5_dp) * dt
      call field(tn, AF(:,2), EF(:,2))
      tn = (tstp + 1.0_dp) * dt
      call field(tn, AF(:,3), EF(:,3))

      Ared(:,1) = ham%get_kreduced(AF(:,1))
      Ared(:,2) = ham%get_kreduced(AF(:,2))
      Ared(:,3) = ham%get_kreduced(AF(:,3))

      allocate(Hk_1(nbnd,nbnd,Nk),Hk_2(nbnd,nbnd,Nk),Hk_3(nbnd,nbnd,Nk))

      call ham_fft%GetHam_Dressed(Ared(:,1), Hk_1)
      call ham_fft%GetHam_Dressed(Ared(:,2), Hk_2)
      call ham_fft%GetHam_Dressed(Ared(:,3), Hk_3)
      if(.not.Peierls_only) then
         allocate(Dk_1(nbnd,nbnd,3,Nk),Dk_2(nbnd,nbnd,3,Nk),Dk_3(nbnd,nbnd,3,Nk))
         call ham_fft%GetDipole_Dressed(Ared(:,1), Dk_1)
         call ham_fft%GetDipole_Dressed(Ared(:,2), Dk_2)
         call ham_fft%GetDipole_Dressed(Ared(:,3), Dk_3)
         !$OMP PARALLEL DO
         do ik=1,Nk
            do idir=1,3
               Hk_1(:,:,ik) = Hk_1(:,:,ik) - EF(idir,1) * Dk_1(:,:,idir,ik)
               Hk_2(:,:,ik) = Hk_2(:,:,ik) - EF(idir,2) * Dk_2(:,:,idir,ik)
               Hk_3(:,:,ik) = Hk_3(:,:,ik) - EF(idir,3) * Dk_3(:,:,idir,ik)
            end do
         end do
         deallocate(Dk_1,Dk_2,Dk_3)
      end if

      call batch_diag%Init(nbnd,nthreads=ham_fft%nthreads)     

      call omp_set_num_threads(ham_fft%nthreads)
      !$OMP PARALLEL PRIVATE(tid,Rho_old,Dscatt,DRhok_step)
      tid = omp_get_thread_num()
      allocate(Rho_old(nbnd,nbnd),Dscatt(nbnd,nbnd))
      allocate(DRhok_step(nbnd,nbnd,4)); DRhok_step=zero

      Dscatt = zero

      !$OMP DO
      do ik=1,Nk
         Rho_old = Rhok(:,:,ik)

         DRhok_step = zero

         call GetScattTerm(nbnd,large_size,Gmm,Beta,Mu,Hk_1(:,:,ik),Rho_old,batch_diag,tid,Dscatt)
         call util_matmul(Hk_1(:,:,ik), Rho_old, DRhok_step(:,:,1), alpha=-iu, large_size=large_size)
         call util_matmul(Rho_old, Hk_1(:,:,ik), DRhok_step(:,:,1), alpha=iu, beta=one, large_size=large_size)         
         DRhok_step(:,:,1) = DRhok_step(:,:,1) + Dscatt

         Rhok(:,:,ik) = Rho_old + 0.5_dp * dt * DRhok_step(:,:,1)
         call GetScattTerm(nbnd,large_size,Gmm,Beta,Mu,Hk_2(:,:,ik),Rhok(:,:,ik),batch_diag,tid,Dscatt)
         call util_matmul(Hk_2(:,:,ik), Rhok(:,:,ik), DRhok_step(:,:,2), alpha=-iu, large_size=large_size)
         call util_matmul(Rhok(:,:,ik), Hk_2(:,:,ik), DRhok_step(:,:,2), alpha=iu, beta=one, large_size=large_size)         
         DRhok_step(:,:,2) = DRhok_step(:,:,2) + Dscatt         

         Rhok(:,:,ik) = Rho_old + 0.5_dp * dt * DRhok_step(:,:,2)
         call GetScattTerm(nbnd,large_size,Gmm,Beta,Mu,Hk_2(:,:,ik),Rhok(:,:,ik),batch_diag,tid,Dscatt)
         call util_matmul(Hk_2(:,:,ik), Rhok(:,:,ik), DRhok_step(:,:,3), alpha=-iu, large_size=large_size)
         call util_matmul(Rhok(:,:,ik), Hk_2(:,:,ik), DRhok_step(:,:,3), alpha=iu, beta=one, large_size=large_size)         
         DRhok_step(:,:,3) = DRhok_step(:,:,3) + Dscatt        

         Rhok(:,:,ik) = Rho_old + dt * DRhok_step(:,:,3)
         call GetScattTerm(nbnd,large_size,Gmm,Beta,Mu,Hk_3(:,:,ik),Rhok(:,:,ik),batch_diag,tid,Dscatt)
         call util_matmul(Hk_3(:,:,ik), Rhok(:,:,ik), DRhok_step(:,:,4), alpha=-iu, large_size=large_size)
         call util_matmul(Rhok(:,:,ik), Hk_3(:,:,ik), DRhok_step(:,:,4), alpha=iu, beta=one, large_size=large_size)         
         DRhok_step(:,:,4) = DRhok_step(:,:,4) + Dscatt     

         Rhok(:,:,ik) = Rho_old + dt/6.0_dp * (DRhok_step(:,:,1) + 2.0_dp * DRhok_step(:,:,2) &
            + 2.0_dp * DRhok_step(:,:,3) + DRhok_step(:,:,4))         

      end do
      !$OMP END DO

      deallocate(Rho_old,Dscatt)
      deallocate(DRhok_step)
      !$OMP END PARALLEL

      deallocate(Hk_1,Hk_2,Hk_3)

      call batch_diag%Clean()

   end subroutine RelaxTimestep_RK4
!--------------------------------------------------------------------------------------
  subroutine RelaxTimestep_Hyb(ham,ham_fft,tstp,dt,field,T1,T2,Beta,Mu,Rhok,Peierls_only,comm_order)
      type(wann90_tb_t),intent(in) :: ham
      type(wann_fft_t),intent(in)  :: ham_fft
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      procedure(vecpot_efield_func),pointer :: field
      real(dp),intent(in)          :: T1
      real(dp),intent(in)          :: T2
      real(dp),intent(in)          :: Beta
      real(dp),intent(in)          :: Mu
      complex(dp),intent(inout)    :: Rhok(:,:,:)
      logical,intent(in)           :: Peierls_only
      integer,intent(in)           :: comm_order
      logical :: large_size
      integer :: nbnd,Nk,ik,icomm,idir
      real(dp) :: tn
      real(dp),dimension(3) :: Gmm
      real(dp),dimension(3) :: AF,EF,Ared
      complex(dp),dimension(comm_order) :: zh,mzh
      complex(dp),allocatable,dimension(:,:) :: Rho_old,Dscatt,Udt
      complex(dp),allocatable,dimension(:,:,:) :: Hk_12,Ckt
      complex(dp),allocatable,dimension(:,:,:,:) :: Dk_12
      type(Batch_Diagonalize_t) :: batch_diag
      integer :: nthreads,tid

      nbnd = ham_fft%nwan
      Nk = ham_fft%nkpts

      large_size = get_large_size(nbnd)

      call assert_shape(Rhok, [nbnd,nbnd,Nk], "RelaxTimestep_RK4", "Rhok")

      Gmm(1) = 1.0_dp / T1
      Gmm(2) = 1.0_dp / T2
      Gmm(3) = Gmm(1) - Gmm(2)


      tn = (tstp + 0.5_dp) * dt
      call field(tn, AF, EF)

      do icomm=1,comm_order
         zh(icomm) = iu * dt / (icomm + 1.0_dp)
         mzh(icomm) = -zh(icomm)
      end do

      Ared = ham%get_kreduced(AF)

      allocate(Hk_12(nbnd,nbnd,Nk))

      call ham_fft%GetHam_Dressed(Ared, Hk_12)
      if(.not.Peierls_only) then
         allocate(Dk_12(nbnd,nbnd,3,Nk))
         call ham_fft%GetDipole_Dressed(Ared, Dk_12)
         do ik=1,Nk
            do idir=1,3
               Hk_12(:,:,ik) = Hk_12(:,:,ik) - EF(idir) * Dk_12(:,:,idir,ik)
            end do
         end do
         deallocate(Dk_12)
      end if
      call batch_diag%Init(nbnd,nthreads=ham_fft%nthreads)     

      call omp_set_num_threads(ham_fft%nthreads)

      !$OMP PARALLEL PRIVATE(tid,Rho_old,Dscatt,Ckt,Udt)
      tid = omp_get_thread_num()
      allocate(Rho_old(nbnd,nbnd),Dscatt(nbnd,nbnd),Udt(nbnd,nbnd))
      allocate(Ckt(nbnd,nbnd,comm_order)); Ckt=zero

      !$OMP DO
      do ik=1,Nk

         call GetScattTerm(nbnd,large_size,Gmm,Beta,Mu,Hk_12(:,:,ik),Rhok(:,:,ik),batch_diag,tid,Dscatt)

         Ckt(:,:,1) = -iu * dt * Dscatt
         do icomm=1,comm_order-1
            call ZGEMM('N','N',nbnd,nbnd,nbnd,zh(icomm),Hk_12(1,1,ik),nbnd,Ckt(1,1,icomm),&
               nbnd,zero,Ckt(1,1,icomm+1),nbnd)
            call ZGEMM('N','N',nbnd,nbnd,nbnd,mzh,Ckt(1,1,icomm),nbnd,Hk_12(1,1,ik),&
               nbnd,one,Ckt(1,1,icomm+1),nbnd)                 
         end do

         Rho_old = Rhok(:,:,ik)
         do icomm=1,comm_order
            Rho_old = Rho_old + iu * Ckt(:,:,icomm)
         end do

         call GenU_CF2(dt,Hk_12(:,:,ik),Udt)
         call UnitaryStepFBW(nbnd,Udt,Rho_old,Rhok(:,:,ik),large_size=large_size)

      end do
      !$OMP END DO

      deallocate(Rho_old,Dscatt,Udt)
      deallocate(Ckt)
      !$OMP END PARALLEL

      deallocate(Hk_12)

      call batch_diag%Clean()

   end subroutine RelaxTimestep_Hyb
!--------------------------------------------------------------------------------------
   subroutine GetScattTerm(nbnd,large_size,Gmm,Beta,Mu,Hk,Rhok,batch_diag,tid,Dscatt)
      integer,intent(in)  :: nbnd
      logical,intent(in)  :: large_size
      real(dp),intent(in) :: Gmm(3)
      real(dp),intent(in) :: Beta
      real(dp),intent(in) :: Mu
      complex(dp),intent(in) :: Hk(:,:)
      complex(dp),intent(in) :: Rhok(:,:)
      type(Batch_Diagonalize_t) :: batch_diag
      integer,intent(in)  :: tid
      complex(dp),intent(inout) :: Dscatt(:,:)
      integer :: i
      real(dp),dimension(nbnd) :: Ek,Occk
      complex(dp),dimension(nbnd,nbnd) :: Rho_eq,Rho_tmp,Rho_off,rotk

      call batch_diag%Diagonalize(Hk, epsk=Ek, vectk=rotk, tid=tid)
      Occk = nfermi(Beta, Ek - Mu)

      Rho_tmp = zero
      do i=1,nbnd
         Rho_tmp(i,i) = Occk(i)
      end do
      Rho_eq = util_rotate_cc(nbnd,rotk,Rho_tmp,large_size=large_size)

      Rho_tmp = util_rotate(nbnd,rotk,Rhok,large_size=large_size)
      do i=1,nbnd
         Rho_tmp(i,i) = zero
      end do
      Rho_off = util_rotate_cc(nbnd,rotk,Rho_tmp,large_size=large_size)

      Dscatt = - Gmm(1) * (Rhok - Rho_eq) + Gmm(3) * Rho_off

   end subroutine GetScattTerm
!--------------------------------------------------------------------------------------

!======================================================================================
end module wan_fft_propagation
module wan_unitrk4
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp, zero, one, iu, nfermi
   use scitools_linalg,only: util_matmul,util_rotate,util_rotate_cc,get_large_size
   use wan_hamiltonian,only: wann90_tb_t
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: Init_RungeKutta, Clean_RungeKutta, RK5_TimeStep_Dip_k, RK5_TimeStep_velo_k

   integer,parameter :: qc = 1

   logical :: large_size
   integer :: nbnd_
   real(dp),allocatable,dimension(:) :: epsk, occk
   integer,allocatable,dimension(:) :: iwork,ifail
   real(dp),allocatable,dimension(:) :: rwork
   complex(dp),allocatable,dimension(:) :: mat_pack, cwork
   complex(dp),allocatable,dimension(:,:) :: Hkt, rotk, Rhok_old, Rho_eq, Rho_off, Rho_tmp
   complex(dp),allocatable,dimension(:,:,:) :: DRhok_step
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init_RungeKutta(nbnd)
      integer,intent(in) :: nbnd

      nbnd_ = nbnd
      large_size = get_large_size(nbnd)

      allocate(iwork(5*nbnd))
      allocate(ifail(nbnd))
      allocate(rwork(7*nbnd))
      allocate(cwork(5*nbnd))
      allocate(epsk(nbnd))
      allocate(mat_pack((nbnd*(nbnd + 1))/2))
      allocate(Hkt(nbnd,nbnd))
      allocate(rotk(nbnd,nbnd))
      allocate(Rhok_old(nbnd,nbnd))
      allocate(Rho_eq(nbnd,nbnd))
      allocate(Rho_tmp(nbnd,nbnd))
      allocate(Rho_off(nbnd,nbnd))
      allocate(DRhok_step(nbnd,nbnd,6))

   end subroutine Init_RungeKutta
!--------------------------------------------------------------------------------------
   subroutine Clean_RungeKutta()

      deallocate(iwork)
      deallocate(ifail)
      deallocate(rwork)
      deallocate(cwork)
      deallocate(epsk)
      deallocate(mat_pack)
      deallocate(Hkt)
      deallocate(rotk)
      deallocate(Rhok_old)
      deallocate(Rho_eq)
      deallocate(Rho_tmp)
      deallocate(Rho_off)
      deallocate(DRhok_step)

   end subroutine Clean_RungeKutta
!--------------------------------------------------------------------------------------
   subroutine RK5_TimeStep_Dip_k(w90,kpt,tstp,dt,field,T1,T2,beta,mu,Rhok,empirical)
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
      interface
         subroutine field(t,A,E)
            import :: dp
            implicit none
            real(dp),intent(in) :: t
            real(dp),intent(out) :: A(3),E(3)
         end subroutine field
      end interface
      logical :: empirical_
      integer :: k,j
      real(dp) :: tk
      real(dp) :: gm(3),AF(3),EF(3)

      empirical_ = .false.
      if(present(empirical)) empirical_ = empirical

      gm(1) = 1.0_dp / T1
      gm(2) = 1.0_dp / T2
      gm(3) = gm(1) - gm(2)

      Rhok_old = Rhok

      tk = tstp * dt
      call field(tk, AF, EF)
      DRhok_step(:,:,1) = deriv_dipole_gauge(nbnd_,w90,kpt,AF,EF,beta,mu,gm,Rhok_old,empirical_)
      do k=1,5
         Rhok = Rhok_old
         do j=1,k
            Rhok = Rhok + dt * CRK5(k,j)*DRhok_step(:,:,j)
         end do

         tk = (tstp + ARK5(k))*dt
         call field(tk, AF, EF)
         DRhok_step(:,:,k+1) = deriv_dipole_gauge(nbnd_,w90,kpt,AF,EF,beta,mu,gm,Rhok,empirical_)
      end do

      Rhok = Rhok_old
      do j=1,6
         Rhok = Rhok + dt*BRK5(j)*DRhok_step(:,:,j)
      end do

   end subroutine RK5_TimeStep_Dip_k
!--------------------------------------------------------------------------------------
   subroutine RK5_TimeStep_velo_k(ik,tstp,dt,field,T1,T2,beta,mu,Hk,vk,Rhok)
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
      interface
         subroutine field(t,A,E)
            import :: dp
            implicit none
            real(dp),intent(in) :: t
            real(dp),intent(out) :: A(3),E(3)
         end subroutine field
      end interface
      integer :: k,j
      real(dp) :: tk
      real(dp) :: gm(3),AF(3),EF(3)

      gm(1) = 1.0_dp / T1
      gm(2) = 1.0_dp / T2
      gm(3) = gm(1) - gm(2)

      Rhok_old = Rhok

      tk = tstp * dt
      call field(tk, AF, EF)
      DRhok_step(:,:,1) = deriv_velo_gauge(nbnd_,ik,Hk,vk,AF,EF,beta,mu,gm,Rhok_old)
      do k=1,5
         Rhok = Rhok_old
         do j=1,k
            Rhok = Rhok + dt * CRK5(k,j)*DRhok_step(:,:,j)
         end do

         tk = (tstp + ARK5(k))*dt
         call field(tk, AF, EF)
         DRhok_step(:,:,k+1) = deriv_velo_gauge(nbnd_,ik,Hk,vk,AF,EF,beta,mu,gm,Rhok_old)
      end do

      Rhok = Rhok_old
      do j=1,6
         Rhok = Rhok + dt*BRK5(j)*DRhok_step(:,:,j)
      end do

   end subroutine RK5_TimeStep_velo_k
!--------------------------------------------------------------------------------------
   function GetDeriv(nst,beta,mu,gm,yt) result(dydt)
      integer,intent(in)           :: nst
      real(dp),intent(in)          :: beta
      real(dp),intent(in)          :: mu
      real(dp),intent(in)          :: gm(3)
      complex(dp),intent(in)       :: yt(:,:)    
      complex(dp) :: dydt(nst,nst) 
      integer :: i

      call diagonalize(Hkt)
      occk = nfermi(Beta,epsk - Mu)

      Rho_tmp = zero
      do i=1,nst
         Rho_tmp(i,i) = occk(i)
      end do
      rho_eq = util_rotate_cc(nst,rotk,Rho_tmp,large_size=large_size)

      Rho_tmp = util_rotate(nst,rotk,yt,large_size=large_size)
      do i=1,nst
         Rho_tmp(i,i) = zero
      end do
      rho_off = util_rotate_cc(nst,rotk,Rho_tmp,large_size=large_size)

      dydt = zero
      call util_matmul(Hkt, yt, dydt, alpha=-iu, large_size=large_size)
      call util_matmul(yt, Hkt, dydt, alpha=iu, beta=one, large_size=large_size)
      dydt = dydt - gm(1) * (yt - Rho_eq)
      dydt = dydt + gm(3) * rho_off

   end function GetDeriv
!--------------------------------------------------------------------------------------
   function deriv_dipole_gauge(nst,w90,kpt,AF,EF,beta,mu,gm,yt,empirical) result(dydt)
      integer,intent(in)           :: nst
      type(wann90_tb_t),intent(in) :: w90
      real(dp),intent(in)          :: kpt(3)
      real(dp),intent(in)          :: AF(3),EF(3)
      real(dp),intent(in)          :: beta
      real(dp),intent(in)          :: mu
      real(dp),intent(in)          :: gm(3)
      complex(dp),intent(in)       :: yt(:,:)
      logical,intent(in)           :: empirical
      complex(dp) :: dydt(nst,nst)

      Hkt = Wann_GetHk_dip(w90,AF,EF,kpt,reducedA=.false.,&
         Peierls_only=empirical)

      dydt = GetDeriv(nst,beta,mu,gm,yt)

   end function deriv_dipole_gauge
!--------------------------------------------------------------------------------------
   function deriv_velo_gauge(nst,ik,Hk,vk,AF,EF,beta,mu,gm,yt) result(dydt)
      integer,intent(in)           :: nst
      integer,intent(in)           :: ik
      complex(dp),intent(in)       :: Hk(:,:,:)
      complex(dp),intent(in)       :: vk(:,:,:,:)
      real(dp),intent(in)          :: AF(3),EF(3)
      real(dp),intent(in)          :: beta
      real(dp),intent(in)          :: mu
      real(dp),intent(in)          :: gm(3)
      complex(dp),intent(in)       :: yt(:,:)
      complex(dp) :: dydt(nst,nst)
      integer :: idir

      Hkt = Hk(:,:,ik)
      do idir=1,3
         Hkt = Hkt - AF(idir) * vk(:,:,ik,idir)
      end do

      dydt = GetDeriv(nst,beta,mu,gm,yt)

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
   subroutine diagonalize(Hk)
      complex(dp),intent(in)  :: Hk(:,:)
      integer :: i,j
      integer :: info,nfound

      do j = 1, nbnd_
         do i = 1, j
            mat_pack(i + ((j - 1)*j)/2) = Hk(i, j)
         end do
      end do
      rotk = zero; epsk = 0.0_dp; cwork = zero; rwork = 0.0_dp; iwork = 0
      call ZHPEVX('V', 'A', 'U', nbnd_, mat_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
         nfound, epsk(1), rotk, nbnd_, cwork, rwork, iwork, ifail, info)
      if (info < 0) then
         write (output_unit, '(a,i3,a)') 'THE ', -info, &
         ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
         write(error_unit,*) 'Error in utility_diagonalize'
      end if
      if (info > 0) then
         write (output_unit, '(i3,a)') info, ' EIGENVECTORS FAILED TO CONVERGE'
         write(error_unit,*) 'Error in utility_diagonalize'
      end if

   end subroutine diagonalize
!--------------------------------------------------------------------------------------

!======================================================================================
end module wan_unitrk4
module wan_scissor
!====================================================================================== 
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
#ifdef OMP
   use omp_lib
#endif 
   use scitools_def,only: dp,zero
   use scitools_linalg,only: EigH, util_rotate_cc
   use wan_hamiltonian,only: wann90_tb_t
   use wan_latt_kpts,only: kpoints_t, GenKgrid
   use wan_fourier,only: fourier_k_to_R
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: ScissorHamiltonian

!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine ScissorHamiltonian(scissor_index,scissor_energy,wann,wann_sci)
      integer,intent(in)            :: scissor_index
      real(dp),intent(in)           :: scissor_energy
      type(wann90_tb_t),intent(in)  :: wann
      type(wann90_tb_t),intent(out) :: wann_sci
      integer :: nx,ny,nz,Nk,Nr,ik,j,ir
      integer,dimension(2) :: ix_bound,iy_bound,iz_bound
      integer,dimension(3) :: irvec
      real(dp),dimension(3) :: kpt
      real(dp),allocatable,dimension(:) :: epsk
      real(dp),allocatable,dimension(:,:) :: kpts
      complex(dp),allocatable,dimension(:,:) :: rotk,diagk
      complex(dp),allocatable,dimension(:,:,:) :: Hk

      call wann_sci%Set(wann)

      if(scissor_index <= 0) return

      ix_bound(1) = minval(wann%irvec(:,1))
      ix_bound(2) = maxval(wann%irvec(:,1))
      iy_bound(1) = minval(wann%irvec(:,2))
      iy_bound(2) = maxval(wann%irvec(:,2))
      iz_bound(1) = minval(wann%irvec(:,3))
      iz_bound(2) = maxval(wann%irvec(:,3))  

      nx = ix_bound(2) - ix_bound(1) + 1   
      ny = iy_bound(2) - iy_bound(1) + 1 
      nz = iz_bound(2) - iz_bound(1) + 1 

      call GenKgrid(kpts,nx,ny,nz,center_bz=.false.)
      Nk = size(kpts, dim=1)

      allocate(Hk(wann%num_wann,wann%num_wann,Nk))

      !$OMP PARALLEL PRIVATE(kpt,epsk,rotk,diagk,irvec)
      allocate(epsk(wann%num_wann))
      allocate(rotk(wann%num_wann,wann%num_wann))
      allocate(diagk(wann%num_wann,wann%num_wann))

      !$OMP DO
      do ik = 1, Nk
         kpt = kpts(ik,:)
         Hk(:,:,ik) = wann%get_ham(kpt)

         call EigH(Hk(:,:,ik), epsk, rotk)

         do j = scissor_index, wann%num_wann
            epsk(j) = epsk(j) + scissor_energy
         end do

         diagk = zero
         do j = 1, wann%num_wann
            diagk(j,j) = epsk(j)
         end do

         Hk(:,:,ik) = util_rotate_cc(wann%num_wann, rotk, diagk, large_size=.true.)

      end do
      !$OMP END DO

      deallocate(epsk)
      deallocate(rotk)
      deallocate(diagk)

      !$OMP DO
      do ir = 1, wann%nrpts
         irvec = wann%irvec(ir,:)
         call fourier_k_to_R(kpts, irvec, wann%ndegen(ir), Hk, wann_sci%ham_r(:,:,ir))
      end do
      !$OMP END DO

      !$OMP END PARALLEL

      deallocate(Hk)
      deallocate(kpts)

   end subroutine ScissorHamiltonian
!--------------------------------------------------------------------------------------

!======================================================================================   
end module wan_scissor

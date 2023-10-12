module wan_compress
!====================================================================================== 
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use scitools_def,only: dp
   use wan_hamiltonian,only: wann90_tb_t
   use wan_overlap,only: wann90_ovlp_t
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: PruneHoppings

   interface PruneHoppings
      module procedure PruneHoppings_ham, PruneHoppings_ovlp
   end interface PruneHoppings
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine PruneHoppings_ham(thresh,wann,wann_pruned,comp_rate)
      real(dp),intent(in)           :: thresh
      type(wann90_tb_t),intent(in)  :: wann
      type(wann90_tb_t),intent(out) :: wann_pruned
      real(dp),intent(out)          :: comp_rate
      integer :: ir,ir_old
      integer,allocatable :: ir_indc(:),ir_indc_pruned(:)
      real(dp),allocatable :: hop(:)

      wann_pruned%num_wann = wann%num_wann
      wann_pruned%real_lattice = wann%real_lattice
      wann_pruned%recip_lattice = wann%recip_lattice
      wann_pruned%recip_reduced = wann%recip_reduced
      wann_pruned%coords_present = wann%coords_present

      allocate(ir_indc(wann%nrpts))
      ir_indc = [ (ir, ir=1,wann%nrpts) ]

      allocate(hop(wann%nrpts))
      do ir=1,wann%nrpts
         hop(ir) = maxval(abs(wann%ham_r(:,:,ir)))
      end do
      wann_pruned%nrpts = count(hop > thresh)
      allocate(ir_indc_pruned(wann_pruned%nrpts))
      ir_indc_pruned = pack(ir_indc, hop > thresh)

      if(allocated(wann_pruned%ndegen)) deallocate(wann_pruned%ndegen)
      if(allocated(wann_pruned%irvec)) deallocate(wann_pruned%irvec)
      if(allocated(wann_pruned%crvec)) deallocate(wann_pruned%crvec)
      if(allocated(wann_pruned%ham_r)) deallocate(wann_pruned%ham_r)
      if(allocated(wann_pruned%pos_r)) deallocate(wann_pruned%pos_r)
      if(allocated(wann_pruned%coords)) deallocate(wann_pruned%coords)

      allocate(wann_pruned%ndegen(wann_pruned%nrpts))
      allocate(wann_pruned%irvec(wann_pruned%nrpts,3))
      allocate(wann_pruned%crvec(wann_pruned%nrpts,3))
      allocate(wann_pruned%ham_r(wann%num_wann,wann%num_wann,wann_pruned%nrpts))
      allocate(wann_pruned%pos_r(wann%num_wann,wann%num_wann,3,wann_pruned%nrpts))
      allocate(wann_pruned%coords(wann%num_wann,3))
      
      wann_pruned%coords = wann%coords

      do ir=1,wann_pruned%nrpts
         ir_old = ir_indc_pruned(ir)
         wann_pruned%ndegen(ir) = wann%ndegen(ir_old)
         wann_pruned%irvec(ir,1:3) = wann%irvec(ir_old,1:3)
         wann_pruned%crvec(ir,1:3) = wann%crvec(ir_old,1:3)
         wann_pruned%ham_r(:,:,ir) = wann%ham_r(:,:,ir_old)
         wann_pruned%pos_r(:,:,:,ir) = wann%pos_r(:,:,:,ir_old)
      end do

      comp_rate = 1.0_dp - wann_pruned%nrpts / dble(wann%nrpts)

      deallocate(ir_indc,hop,ir_indc_pruned)

   end subroutine PruneHoppings_ham
!--------------------------------------------------------------------------------------
   subroutine PruneHoppings_ovlp(thresh,ovlp,ovlp_pruned,comp_rate)
      real(dp),intent(in)             :: thresh
      type(wann90_ovlp_t),intent(in)  :: ovlp
      type(wann90_ovlp_t),intent(out) :: ovlp_pruned
      real(dp),intent(out)            :: comp_rate
      integer :: ir,ir_old
      integer,allocatable :: ir_indc(:),ir_indc_pruned(:)
      real(dp),allocatable :: hop(:)

      ovlp_pruned%num_wann = ovlp%num_wann
      ovlp_pruned%real_lattice = ovlp%real_lattice

      allocate(ir_indc(ovlp%nrpts))
      ir_indc = [ (ir, ir=1,ovlp%nrpts) ]

      allocate(hop(ovlp%nrpts))
      do ir=1,ovlp%nrpts
         hop(ir) = maxval(abs(ovlp%S_r(:,:,ir)))
      end do
      ovlp_pruned%nrpts = count(hop > thresh)
      allocate(ir_indc_pruned(ovlp_pruned%nrpts))
      ir_indc_pruned = pack(ir_indc, hop > thresh)

      if(allocated(ovlp_pruned%ndegen)) deallocate(ovlp_pruned%ndegen)
      if(allocated(ovlp_pruned%irvec)) deallocate(ovlp_pruned%irvec)
      if(allocated(ovlp_pruned%S_r)) deallocate(ovlp_pruned%S_r)
  
      allocate(ovlp_pruned%ndegen(ovlp_pruned%nrpts))
      allocate(ovlp_pruned%irvec(ovlp_pruned%nrpts,3))
      allocate(ovlp_pruned%S_r(ovlp%num_wann,ovlp%num_wann,ovlp_pruned%nrpts))

      do ir=1,ovlp_pruned%nrpts
         ir_old = ir_indc_pruned(ir)
         ovlp_pruned%ndegen(ir) = ovlp%ndegen(ir_old)
         ovlp_pruned%irvec(ir,1:3) = ovlp%irvec(ir_old,1:3)
         ovlp_pruned%S_r(:,:,ir) = ovlp%S_r(:,:,ir_old)
      end do

      comp_rate = 1.0_dp - ovlp_pruned%nrpts / dble(ovlp%nrpts)

      deallocate(ir_indc,hop,ir_indc_pruned)

   end subroutine PruneHoppings_ovlp
!--------------------------------------------------------------------------------------

!======================================================================================   
end module wan_compress

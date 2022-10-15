module Mwann_slab
!======================================================================================
   use Mdebug
   use scitools_def,only: dp,zero
   use Mlatt_utils,only: utility_recip_lattice, utility_recip_reduced
   use Mham_w90,only: wann90_tb_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: Wannier_BulkToSlab
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Wannier_BulkToSlab(bulk_w90,nlayer,slab_w90,ijmax)
      type(wann90_tb_t),intent(in)  :: bulk_w90
      integer,intent(in)            :: nlayer
      type(wann90_tb_t),intent(out) :: slab_w90
      integer,intent(in),optional   :: ijmax
      integer :: ijmax_
      integer :: nwan,nrpts_2d,irpt,irpt2d
      integer :: i1,i2,i3,i,j,idir
      real(dp) :: e3p(3),oop_vec(3),Ua(3,3)
      complex(dp),allocatable :: Hij(:,:,:,:),pos_r(:,:,:,:)
      complex(dp),allocatable :: Dij(:,:,:,:,:)

      ijmax_ = 10
      if(present(ijmax)) ijmax_ = ijmax

      nwan = bulk_w90%num_wann
      slab_w90%real_lattice = bulk_w90%real_lattice
      call GetDipoleRotation(bulk_w90%real_lattice,e3p,Ua)
      slab_w90%real_lattice(3,:) = e3p

      ! call GetOutOfPlaneVector(bulk_w90%real_lattice,oop_vec)
      !!! CHECK
      oop_vec = e3p * dot_product(bulk_w90%real_lattice(3,:),  e3p)

      call utility_recip_lattice(slab_w90%real_lattice, slab_w90%recip_lattice)
      call utility_recip_reduced(slab_w90%recip_lattice, slab_w90%recip_reduced)

      slab_w90%num_wann = nlayer * bulk_w90%num_wann

      nrpts_2d =&
         (maxval(bulk_w90%irvec(:,1))-minval(bulk_w90%irvec(:,1))+1) &
         *(maxval(bulk_w90%irvec(:,2))-minval(bulk_w90%irvec(:,2))+1) 

      slab_w90%nrpts = nrpts_2d
      allocate(slab_w90%irvec(nrpts_2d,3))
      allocate(slab_w90%ndegen(nrpts_2d))

      !! H_{ij} [num_wann, num_wann, (2*ijmax+1), nrpts_2d]
      allocate(Hij(bulk_w90%num_wann,bulk_w90%num_wann,2*ijmax_+1,nrpts_2d))
      allocate(Dij(bulk_w90%num_wann,bulk_w90%num_wann,2*ijmax_+1,nrpts_2d,3))
      allocate(pos_r(bulk_w90%num_wann,bulk_w90%num_wann,bulk_w90%nrpts,3))
      Hij = zero
      Dij = zero

      ! rotate dipoles
      do irpt=1,bulk_w90%nrpts
         do j=1,bulk_w90%num_wann
            do i=1,bulk_w90%num_wann
               pos_r(i,j,irpt,1:3) = matmul(bulk_w90%pos_r(i,j,irpt,1:3),Ua)
            end do
         end do
      end do

      slab_w90%irvec = 0
      slab_w90%ndegen = 1
      do irpt=1,bulk_w90%nrpts
         i1 = bulk_w90%irvec(irpt,1)
         i2 = bulk_w90%irvec(irpt,2)
         i3 = bulk_w90%irvec(irpt,3)    
         if(abs(i3) < ijmax_) then
            irpt2d = get_index_i3(bulk_w90,i1,i2)
            Hij(:,:,i3+ijmax_+1,irpt2d) = bulk_w90%ham_r(:,:,irpt)
            Dij(:,:,i3+ijmax_+1,irpt2d,1) = pos_r(:,:,irpt,1)
            slab_w90%irvec(irpt2d,1) = i1
            slab_w90%irvec(irpt2d,2) = i2
            slab_w90%irvec(irpt2d,3) = i3
         end if     
      end do

      allocate(slab_w90%ham_r(slab_w90%num_wann,slab_w90%num_wann,nrpts_2d))
      allocate(slab_w90%pos_r(slab_w90%num_wann,slab_w90%num_wann,nrpts_2d,3))
      do irpt2d=1,nrpts_2d
         do i=1,nlayer
            do j=1,nlayer
               slab_w90%ham_r((j-1)*nwan+1:j*nwan,(i-1)*nwan+1:i*nwan,irpt2d) = &
                  Hij(:,:,i-j+ijmax_+1,irpt2d)
               slab_w90%pos_r((j-1)*nwan+1:j*nwan,(i-1)*nwan+1:i*nwan,irpt2d,1) = &
                  Dij(:,:,i-j+ijmax_+1,irpt2d,1)
               slab_w90%pos_r((j-1)*nwan+1:j*nwan,(i-1)*nwan+1:i*nwan,irpt2d,2) = &
                  Dij(:,:,i-j+ijmax_+1,irpt2d,2)
               slab_w90%pos_r((j-1)*nwan+1:j*nwan,(i-1)*nwan+1:i*nwan,irpt2d,3) = &
                  Dij(:,:,i-j+ijmax_+1,irpt2d,3)
            end do
         end do
      end do

      if(bulk_w90%coords_present) then
         slab_w90%coords_present = .true.

         allocate(slab_w90%coords(slab_w90%num_wann,3))
         do i=1,nlayer
            slab_w90%coords((i-1)*nwan+1:i*nwan,1) = bulk_w90%coords(1:nwan,1) &
               - (i-1) * oop_vec(1)
            slab_w90%coords((i-1)*nwan+1:i*nwan,2) = bulk_w90%coords(1:nwan,2) &
               - (i-1) * oop_vec(2)
            slab_w90%coords((i-1)*nwan+1:i*nwan,3) = bulk_w90%coords(1:nwan,3) &
               - (i-1) * oop_vec(3)
         end do
      end if

      deallocate(pos_r,Hij,Dij)

   end subroutine Wannier_BulkToSlab
!--------------------------------------------------------------------------------------
   subroutine GetDipoleRotation(real_lattice,e3p,Ua)
      real(dp),intent(in)    :: real_lattice(3,3)
      real(dp),intent(out)   :: e3p(3)
      real(dp),intent(out)   :: Ua(3,3)
      real(dp),dimension(3) :: e1,e2,e3,e1p,e2p

      e1 = real_lattice(1,:) / norm2(real_lattice(1,:))
      e2 = real_lattice(2,:) / norm2(real_lattice(2,:))
      e3 = real_lattice(3,:) / norm2(real_lattice(3,:))

      e1p = e1
      e2p = e2 - dot_product(e1p,e2) * e1p
      e3p = cross(e1p,e2p)
      e3p = e3p / norm2(e3p)

      Ua(1:3,1) = e1p
      Ua(1:3,2) = e2p
      Ua(1:3,3) = e3p

   end subroutine GetDipoleRotation
!--------------------------------------------------------------------------------------
   function get_index_i3(wann,i1,i2) result(index_2d_i3)
      !! For some i3, index the 2D nrpts by following method
      !! 1 : (min_1,min_2); 2:(min_1, min_2+1) .....
      !! max_1-min_1+2 : (min_1+1,min_2)....
      type(wann90_tb_t),intent(in) :: wann
      integer,intent(in) :: i1,i2
      integer :: min_1,max_1,min_2,max_2,index_2d_i3

      min_1 = minval(wann%irvec(:,1))
      max_1 = maxval(wann%irvec(:,1))
      min_2 = minval(wann%irvec(:,2))
      max_2 = maxval(wann%irvec(:,2))

      index_2d_i3 = (i1 - min_1) * (max_2 - min_2 + 1) + (i2 - min_2 + 1)

    end function get_index_i3
!--------------------------------------------------------------------------------------
   function cross(a, b) result(cross_vec)
      real(dp),intent(in) :: a(3), b(3)
      real(dp) :: cross_vec(3)

      cross_vec(1) = a(2) * b(3) - a(3) * b(2)
      cross_vec(2) = a(3) * b(1) - a(1) * b(3)
      cross_vec(3) = a(1) * b(2) - a(2) * b(1)

   end function cross
!--------------------------------------------------------------------------------------

!======================================================================================
end module Mwann_slab
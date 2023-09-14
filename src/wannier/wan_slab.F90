module wan_slab
!! Provides functions to extend a bulk Wannier Hamiltonian to construct a slab.
!======================================================================================
   use Mdebug
   use scitools_def,only: dp,zero
   use wan_utils,only: utility_recip_lattice, utility_recip_reduced
   use wan_hamiltonian,only: wann90_tb_t
   use wan_overlap,only: wann90_ovlp_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: Wannier_BulkToSlab, Overlap_BulkToSlab
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Wannier_BulkToSlab(bulk_w90,nlayer,slab_w90)
   !! Given a bulk Wannier Hamiltonian, constructs a slab Wannier Hamiltonian.
   !! .We assume the slab is constructed along the 
   !! \(\mathbf{c} = \mathbf{a}\_1\times \mathbf{a}\_2\) direction, where 
   !! \(\mathbf{a}_{1,2}\) are the in-plane lattice vectors of the bulk Hamiltonian.
   !! Restriction: we assume the \(\mathbf{a}_{1,2}\) are in the \(x\)-\(y\) plane:
   !! \(\mathbf{a}_{1,2} \cdot \hat{z} = 0\), where \(\hat{z}\) is the unit vector
   !! in \(z\)-direction.
   !! 
   !! The slab Hamiltonian is equivalent to a two-dimensional system with enlarged 
   !! orbital space. By defining the new orbital indicies \(m = (j l)\)  with 
   !! \(j\) denoting the orbitals of the bulk Hamiltonian and \(l\) representing the 
   !! layer index, the slab Hamiltonian can be represented in the same standard format as
   !! any Wannier Hamiltonian
      type(wann90_tb_t),intent(in)  :: bulk_w90 !! The bulk Wannier Hamiltonian
      integer,intent(in)            :: nlayer !! number of layers
      type(wann90_tb_t),intent(out) :: slab_w90 !! The slab Wannier Hamiltonian with 
                                                !! enlarged orbital space.
      integer :: ijmax
      integer :: nwan,nrpts_2d,irpt,irpt2d
      integer :: i1,i2,i3,i,j,idir
      real(dp) :: e3p(3),oop_vec(3),Ua(3,3),pos_r_slab(3)
      complex(dp),allocatable :: Hij(:,:,:,:),pos_r(:,:,:,:)
      complex(dp),allocatable :: Dij(:,:,:,:,:)

      ijmax = nint(nlayer/2.0_dp) + 1

      nwan = bulk_w90%num_wann
      slab_w90%real_lattice = bulk_w90%real_lattice
      call GetDipoleRotation(bulk_w90%real_lattice,e3p,Ua)
      ! slab_w90%real_lattice = matmul(bulk_w90%real_lattice, Ua)
      ! slab_w90%real_lattice(3,:) = e3p

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
      allocate(Hij(bulk_w90%num_wann,bulk_w90%num_wann,2*ijmax+1,nrpts_2d))
      allocate(Dij(bulk_w90%num_wann,bulk_w90%num_wann,2*ijmax+1,nrpts_2d,3))
      allocate(pos_r(bulk_w90%num_wann,bulk_w90%num_wann,bulk_w90%nrpts,3))
      Hij = zero
      Dij = zero

      ! rotate dipoles
      do irpt=1,bulk_w90%nrpts
         do j=1,bulk_w90%num_wann
            do i=1,bulk_w90%num_wann
               ! pos_r(i,j,irpt,1:3) = matmul(bulk_w90%pos_r(i,j,irpt,1:3),Ua)
               pos_r(i,j,irpt,1:3) = bulk_w90%pos_r(i,j,irpt,1:3)
            end do
         end do
      end do

      slab_w90%irvec = 0
      slab_w90%ndegen = 1
      do irpt=1,bulk_w90%nrpts
         i1 = bulk_w90%irvec(irpt,1)
         i2 = bulk_w90%irvec(irpt,2)
         i3 = bulk_w90%irvec(irpt,3)    
         if(abs(i3) < ijmax) then
            irpt2d = get_index_i3(bulk_w90%irvec,i1,i2)
            Hij(:,:,i3+ijmax+1,irpt2d) = bulk_w90%ham_r(:,:,irpt)
            Dij(:,:,i3+ijmax+1,irpt2d,1) = pos_r(:,:,irpt,1)
            Dij(:,:,i3+ijmax+1,irpt2d,2) = pos_r(:,:,irpt,2)
            Dij(:,:,i3+ijmax+1,irpt2d,3) = pos_r(:,:,irpt,3)
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
                  Hij(:,:,i-j+ijmax+1,irpt2d)
               do idir=1,3
                  slab_w90%pos_r((j-1)*nwan+1:j*nwan,(i-1)*nwan+1:i*nwan,irpt2d,idir) = &
                     Dij(:,:,i-j+ijmax+1,irpt2d,idir)
               end do
            end do
         end do
      end do

      do irpt=1,slab_w90%nrpts
         do i=1,nlayer
            slab_w90%pos_r((i-1)*nwan+1:i*nwan,(i-1)*nwan+1:i*nwan,irpt,3) = &
               slab_w90%pos_r((i-1)*nwan+1:i*nwan,(i-1)*nwan+1:i*nwan,irpt,3) - (i-1) * norm2(oop_vec)
         end do
      end do

      if(bulk_w90%coords_present) then
         slab_w90%coords_present = .true.

         allocate(slab_w90%coords(slab_w90%num_wann,3))
         do i=1,nlayer
            ! slab_w90%coords((i-1)*nwan+1:i*nwan,1:3) = matmul(bulk_w90%coords(1:nwan,1:3),Ua)
            slab_w90%coords((i-1)*nwan+1:i*nwan,1:3) = bulk_w90%coords(1:nwan,1:3)
            do idir=1,3
               slab_w90%coords((i-1)*nwan+1:i*nwan,idir) = slab_w90%coords((i-1)*nwan+1:i*nwan,idir) &
                  - (i-1) * bulk_w90%real_lattice(3,idir)
            end do
         end do
      end if

      deallocate(pos_r,Hij,Dij)

   end subroutine Wannier_BulkToSlab
!--------------------------------------------------------------------------------------
   subroutine Overlap_BulkToSlab(bulk_w90,nlayer,slab_w90)
   !! Given a bulk Wannier overlap matrix, constructs a slab Wannier overlap matrix.
   !! .We assume the slab is constructed along the 
   !! \(\mathbf{c} = \mathbf{a}\_1\times \mathbf{a}\_2\) direction, where 
   !! \(\mathbf{a}_{1,2}\) are the in-plane lattice vectors of the bulk Hamiltonian.
   !! Restriction: we assume the \(\mathbf{a}_{1,2}\) are in the \(x\)-\(y\) plane:
   !! \(\mathbf{a}_{1,2} \cdot \hat{z} = 0\), where \(\hat{z}\) is the unit vector
   !! in \(z\)-direction.
      type(wann90_ovlp_t),intent(in)  :: bulk_w90 !! The bulk Wannier overlap
      integer,intent(in)              :: nlayer   !! number of layers
      type(wann90_ovlp_t),intent(out) :: slab_w90 !! The slab Wannier overlap with 
                                                  !! enlarged orbital space.
      integer :: ijmax
      integer :: nwan,nrpts_2d,irpt,irpt2d
      integer :: i1,i2,i3,i,j,idir
      real(dp) :: e3p(3),oop_vec(3),Ua(3,3),pos_r_slab(3)
      complex(dp),allocatable :: Sij(:,:,:,:)

      ijmax = nint(nlayer/2.0_dp) + 1

      nwan = bulk_w90%num_wann
      slab_w90%real_lattice = bulk_w90%real_lattice
      call GetDipoleRotation(bulk_w90%real_lattice,e3p,Ua)
      ! slab_w90%real_lattice = matmul(bulk_w90%real_lattice, Ua)
      ! slab_w90%real_lattice(3,:) = e3p

      ! call GetOutOfPlaneVector(bulk_w90%real_lattice,oop_vec)
      !!! CHECK
      oop_vec = e3p * dot_product(bulk_w90%real_lattice(3,:),  e3p)

      slab_w90%num_wann = nlayer * bulk_w90%num_wann

      nrpts_2d =&
         (maxval(bulk_w90%irvec(:,1))-minval(bulk_w90%irvec(:,1))+1) &
         *(maxval(bulk_w90%irvec(:,2))-minval(bulk_w90%irvec(:,2))+1) 

      slab_w90%nrpts = nrpts_2d
      allocate(slab_w90%irvec(nrpts_2d,3))
      allocate(slab_w90%ndegen(nrpts_2d))

      !! S_{ij} [num_wann, num_wann, (2*ijmax+1), nrpts_2d]
      allocate(Sij(bulk_w90%num_wann,bulk_w90%num_wann,2*ijmax+1,nrpts_2d))
      Sij = zero

      slab_w90%irvec = 0
      slab_w90%ndegen = 1
      do irpt=1,bulk_w90%nrpts
         i1 = bulk_w90%irvec(irpt,1)
         i2 = bulk_w90%irvec(irpt,2)
         i3 = bulk_w90%irvec(irpt,3)    
         if(abs(i3) < ijmax) then
            irpt2d = get_index_i3(bulk_w90%irvec,i1,i2)
            Sij(:,:,i3+ijmax+1,irpt2d) = bulk_w90%S_r(:,:,irpt)
            slab_w90%irvec(irpt2d,1) = i1
            slab_w90%irvec(irpt2d,2) = i2
            slab_w90%irvec(irpt2d,3) = i3
         end if     
      end do

      allocate(slab_w90%S_r(slab_w90%num_wann,slab_w90%num_wann,nrpts_2d))
      do irpt2d=1,nrpts_2d
         do i=1,nlayer
            do j=1,nlayer
               slab_w90%S_r((j-1)*nwan+1:j*nwan,(i-1)*nwan+1:i*nwan,irpt2d) = &
                  Sij(:,:,i-j+ijmax+1,irpt2d)
            end do
         end do
      end do

      deallocate(Sij)

   end subroutine Overlap_BulkToSlab
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
   function get_index_i3(irvec,i1,i2) result(index_2d_i3)
      !! For some i3, index the 2D nrpts by following method
      !! 1 : (min_1,min_2); 2:(min_1, min_2+1) .....
      !! max_1-min_1+2 : (min_1+1,min_2)....
      integer,intent(in) :: irvec(:,:)
      integer,intent(in) :: i1,i2
      integer :: min_1,max_1,min_2,max_2,index_2d_i3

      min_1 = minval(irvec(:,1))
      max_1 = maxval(irvec(:,1))
      min_2 = minval(irvec(:,2))
      max_2 = maxval(irvec(:,2))

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
end module wan_slab
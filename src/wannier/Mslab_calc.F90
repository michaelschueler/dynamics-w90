module Mslab_calc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp,zero,one
   use Mlinalg,only: EigHE,inv
   use Mlatt_kpts,only: Read_Kpoints
   use Mwannier_calc,only: wannier_calc_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: BulkToSlab
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine BulkToSlab(wann,file_inp,U,rot_option)
      type(wannier_calc_t),intent(inout) :: wann
      character(len=*),intent(in) :: file_inp
      real(dp),intent(in),optional :: U(3,3)
      integer,intent(in),optional :: rot_option
      integer :: nslab, gauge_slab
      integer :: unit_inp
      logical :: calc_slab,calc_rotation,calc_slab_berry,calc_slab_evecs,calc_slab_oam,&
          calc_slab_orbweight,calc_slab_spin
      namelist/SLABCALCOPT/calc_slab,calc_rotation,calc_slab_berry,calc_slab_evecs,&
          calc_slab_oam,calc_slab_orbweight,calc_slab_spin,gauge_slab,nslab
      integer :: ik
      complex(dp),allocatable :: Hk_s(:,:)
      real(dp) :: Uinv(3,3)

      !! -------------------------------------------------------
      !! Construct the slab hamiltonian from the bulk hamiltonian
      !! Both belong to the class wann90_tb_t
      !! -------------------------------------------------------

      open(newunit=unit_inp, file=trim(file_inp), status='OLD', action='READ')
      read(unit_inp, nml=SLABCALCOPT)
      close(unit_inp)

      !! Initiallize all the wann%ham%parameters
      wann%ham_s%real_lattice      = wann%ham%real_lattice
      wann%ham_s%num_wann          = wann%ham%num_wann * nslab
      wann%ham_s%nslab             = nslab
      wann%ham_s%ijmax             = 10

      !! rotation operation
      if (present(U)) then
          call Rotation_Slab(wann, U, rot_option)
      end if

      !! Initiallize all the slab_class%parameters
      wann%nwan     = wann%nwan * wann%ham_s%nslab
      wann%nbnd     = wann%ham_s%num_wann

      if(allocated(wann%epsk)) deallocate(wann%epsk)
      if(allocated(wann%vectk)) deallocate(wann%vectk)

      !! Initialize slab Hamiltonian
      call Init_Slab_Ham_r(wann)
      
      !! Initialize slab Dipole Matrix element
      call Init_Slab_Pos_r(wann)

      !! Slab k path : use the first two terms in kpts 
      call Read_Kpoints(file_inp,wann%kpts,print_info=.true.)
      wann%Nk = size(wann%kpts,1) 

      allocate(wann%epsk(wann%nbnd,wann%Nk),&
          wann%vectk(wann%nbnd,wann%nbnd,wann%Nk))

      allocate(Hk_s(wann%nbnd,wann%nbnd))
      do ik=1,wann%Nk
        !Hk_s = bulk_class%Ham%get_ham_slab(slab_class%kpts(ik,:),slab_class%ham%nslab)
        Hk_s = wann%Ham_s%get_ham(wann%kpts(ik,:))
        call EigHE(Hk_s,wann%epsk(:,ik),wann%vectk(:,:,ik))
      end do
      deallocate(Hk_s)
      
   end subroutine BulkToSlab
!--------------------------------------------------------------------------------------
   subroutine Rotation_Slab(wann,U,rot_option)
      type(wannier_calc_t),intent(inout) :: wann
      real(dp),intent(in) :: U(3,3)
      integer,intent(in) :: rot_option
      real(dp) :: Uinv(3,3),a2p(3),a3p(3),Ua(3,3)
      integer :: irpt,i1,i2,i
      real(dp) :: n1, n2, n3, n2p, n3p
      
      !! ------------------------------------------------------
      !! rotation operation for slab construction
      !! Initialize all the rotation in this subroutine
      !! rot_option = 1
      !! Rotate the frame / Construct the slab along a1' x a2'
      !! 
      !! ------------------------------------------------------

      !! irvec
      if (rot_option == 1) then
          Uinv = inv(U)
          wann%ham%irvec(:,1:3) = matmul(wann%ham%irvec(:,1:3), Uinv)
          wann%ham_s%real_lattice(:,1:3) = matmul(wann%ham_s%real_lattice(:,1:3),U)
          wann%ham%real_lattice(:,1:3) = matmul(wann%ham%real_lattice(:,1:3),U)
      else if (rot_option == 2) then
          wann%ham_s%real_lattice(:,1:3) = matmul(wann%ham_s%real_lattice(:,1:3),U)
          wann%ham%real_lattice(:,1:3) = matmul(wann%ham%real_lattice(:,1:3),U)
      end if
      
      n1 = sqrt(dot_product(wann%ham%real_lattice(1,:),wann%ham%real_lattice(1,:)))
      n2 = sqrt(dot_product(wann%ham%real_lattice(2,:),wann%ham%real_lattice(2,:)))
      n3 = sqrt(dot_product(wann%ham%real_lattice(3,:),wann%ham%real_lattice(3,:)))

      wann%ham%real_lattice(1,:) = wann%ham%real_lattice(1,:) / n1
      wann%ham%real_lattice(2,:) = wann%ham%real_lattice(2,:) / n2
      wann%ham%real_lattice(3,:) = wann%ham%real_lattice(3,:) / n3

      a2p(:) = wann%ham%real_lattice(2,:) -&
          dot_product(wann%ham%real_lattice(1,:),wann%ham%real_lattice(2,:))*&
          wann%ham%real_lattice(1,:)
      n2p = sqrt(dot_product(a2p, a2p))
      a2p(:) = a2p(:) / n2p

      a3p(:) = cross(wann%ham%real_lattice(1,:),a2p(:))
      n3p = sqrt(dot_product(a3p, a3p))
      a3p(:) = a3p(:) / n3p

      !! Construct Ua matrix
      if (rot_option == 1) then
          do i=1,3
            Ua(i,1) = wann%ham%real_lattice(1,i)
            Ua(i,2) = a2p(i)
            Ua(i,3) = a3p(i)
          end do
      else if (rot_option == 2) then
          Ua = U
      end if

      !! dipole matrix element
      do irpt=1,wann%ham%nrpts
        do i1=1,wann%ham%num_wann
            do i2=1,wann%ham%num_wann
                wann%ham%pos_r(i1,i2,irpt,1:3) = &
                    matmul(wann%ham%pos_r(i1,i2,irpt,1:3),Ua)
            end do
        end do
      end do

      end subroutine Rotation_Slab
!--------------------------------------------------------------------------------------
   subroutine Init_Slab_Ham_r(wann)
      type(wannier_calc_t),intent(inout) :: wann
      complex(kind=dp),allocatable :: Hij(:,:,:,:)
      integer :: irpt,irpt2d,i1,i2,i3,i,j,nrpts_2d,bnum_wann

      !! -------------------------------------------------------
      !! Initialize the slab hamiltonian and store in ham%ham_r
      !! To do :
      !! (1) Combine two steps, no need for Hij
      !! (2) Smarter wat to get irpt2d?
      !! (3) get rid of all the _slab procedure in ham.
      !! (4) use different ham instead of wann_calc?
      !! -------------------------------------------------------

      !! Number of nrpts at some certain i3
      nrpts_2d =&
          (maxval(wann%ham%irvec(:,1))-minval(wann%ham%irvec(:,1))+1)*&
          (maxval(wann%ham%irvec(:,2))-minval(wann%ham%irvec(:,2))+1)
      wann%ham_s%nrpts = nrpts_2d

      !! H_{ij} [num_wann, num_wann, (2*ijmax+1), nrpts_2d]
      allocate(Hij(wann%ham%num_wann,wann%ham%num_wann,2*wann%ham_s%ijmax+1,nrpts_2d))
      !! Initiallize the irvec
      allocate(wann%ham_s%irvec(nrpts_2d,3))
      allocate(wann%ham_s%ndegen(nrpts_2d))
      
      Hij(:,:,:,:) = zero
      wann%ham_s%irvec(:,:) = zero
      wann%ham_s%ndegen(:) = one    !! This is an assumption

      do irpt=1,wann%ham%nrpts
         i1 = wann%ham%irvec(irpt,1)
         i2 = wann%ham%irvec(irpt,2)
         i3 = wann%ham%irvec(irpt,3)
         if (abs(i3) < wann%ham_s%ijmax) then
             irpt2d = get_index_i3(wann,i1,i2)
             Hij(:,:,i3+wann%ham_s%ijmax+1,irpt2d) = wann%ham%ham_r(:,:,irpt)
             wann%ham_s%irvec(irpt2d,1) = i1
             wann%ham_s%irvec(irpt2d,2) = i2
             wann%ham_s%irvec(irpt2d,3) = i3
             !wann%ham_s%ndegen(irpt2d) = 1    
         end if
      end do
            
      !! Now you have to take care the order of the irpt2d_i3
      !! is the same for every i3, for example ;
      !! irpt2d_i3(1) : (-9,-4,-ijmax), (-9,-3,-ijmax),....
      !! irpt2d_i3(2) : (-9,-4,-ijmax+1), (-9,-3,-ijmax+1),....
      !! irpt2d_i3(3) : (-9,-4,-ijmax+2), (-9,-3,-ijmax+2),....
      !! ...etc.

      bnum_wann = wann%ham%num_wann
      allocate(wann%ham_s%ham_r(wann%nbnd,wann%nbnd,nrpts_2d))
      do irpt2d=1,nrpts_2d
        do i=1,wann%ham_s%nslab
            do j=1,wann%ham_s%nslab
                wann%ham_s%ham_r(&
                    (j-1)*bnum_wann+1:j*bnum_wann,(i-1)*bnum_wann+1:i*bnum_wann,irpt2d)=&
                    Hij(:,:,i-j+wann%ham_s%ijmax+1,irpt2d)
            end do
        end do
       end do
      
   end subroutine Init_Slab_Ham_r
!--------------------------------------------------------------------------------------
   subroutine Init_Slab_Pos_r(wann)
      type(wannier_calc_t),intent(inout) :: wann
      complex(kind=dp),allocatable :: Dij(:,:,:,:,:)
      integer :: irpt,irpt2d,i1,i2,i3,i,j,nrpts_2d,bnum_wann
      
      !! -------------------------------------------------------
      !! Initialize the slab dipole ME and store in ham%pos_r
      !! To do :
      !! (1) Combine two steps, no need for Dij
      !! (2) Smarter wat to get irpt2d?
      !! (3) get rid of all the _slab procedure in ham.
      !! (4) use different ham instead of wann_calc?
      !! -------------------------------------------------------

      !! Number of nrpts at some certain i3
      nrpts_2d =&
          (maxval(wann%ham%irvec(:,1))-minval(wann%ham%irvec(:,1))+1)*&
          (maxval(wann%ham%irvec(:,2))-minval(wann%ham%irvec(:,2))+1)
      
      !! D_{ij} [num_wann, num_wann, (2*ijmax+1), nrpts_2d, 3]
      allocate(Dij(wann%ham%num_wann,wann%ham%num_wann,2*wann%ham_s%ijmax+1,nrpts_2d,3))

      Dij(:,:,:,:,:) = zero

      do irpt=1,wann%ham%nrpts
         i1 = wann%ham%irvec(irpt,1)
         i2 = wann%ham%irvec(irpt,2)
         i3 = wann%ham%irvec(irpt,3)
         if (abs(i3) < wann%ham_s%ijmax) then
             irpt2d = get_index_i3(wann,i1,i2)
             Dij(:,:,i3+wann%ham_s%ijmax+1,irpt2d,1) = wann%ham%pos_r(:,:,irpt,1)
             Dij(:,:,i3+wann%ham_s%ijmax+1,irpt2d,2) = wann%ham%pos_r(:,:,irpt,2)
             Dij(:,:,i3+wann%ham_s%ijmax+1,irpt2d,3) = wann%ham%pos_r(:,:,irpt,3)
         end if
      end do

      bnum_wann = wann%ham%num_wann
      allocate(wann%ham_s%pos_r(wann%nbnd,wann%nbnd,nrpts_2d,3))
      do irpt2d=1,nrpts_2d
        do i=1,wann%ham_s%nslab
            do j=1,wann%ham_s%nslab
                wann%ham_s%pos_r(&
                    (j-1)*bnum_wann+1:j*bnum_wann,(i-1)*bnum_wann+1:i*bnum_wann,irpt2d,1)=&
                    Dij(:,:,i-j+wann%ham_s%ijmax+1,irpt2d,1)
                wann%ham_s%pos_r(&
                    (j-1)*bnum_wann+1:j*bnum_wann,(i-1)*bnum_wann+1:i*bnum_wann,irpt2d,2)=&
                    Dij(:,:,i-j+wann%ham_s%ijmax+1,irpt2d,2)
                wann%ham_s%pos_r(&
                    (j-1)*bnum_wann+1:j*bnum_wann,(i-1)*bnum_wann+1:i*bnum_wann,irpt2d,3)=&
                    Dij(:,:,i-j+wann%ham_s%ijmax+1,irpt2d,3)
            end do
        end do
       end do
      
   end subroutine Init_Slab_Pos_r
!--------------------------------------------------------------------------------------
   function get_index_i3(wann,i1,i2) result(index_2d_i3)
      type(wannier_calc_t),intent(in) :: wann
      integer,intent(in) :: i1,i2
      integer :: min_1,max_1,min_2,max_2,index_2d_i3

      !! -------------------------------------------------------
      !! For some i3, index the 2D nrpts by following method
      !! 1 : (min_1,min_2); 2:(min_1, min_2+1) .....
      !! max_1-min_1+2 : (min_1+1,min_2)....
      !! -------------------------------------------------------

      min_1 = minval(wann%ham%irvec(:,1))
      max_1 = maxval(wann%ham%irvec(:,1))
      min_2 = minval(wann%ham%irvec(:,2))
      max_2 = maxval(wann%ham%irvec(:,2))

      index_2d_i3 = (i1 - min_1) * (max_2 - min_2 + 1) + (i2 - min_2 + 1)

    end function get_index_i3
!--------------------------------------------------------------------------------------
    FUNCTION cross(a, b) result(cross_vec)
      real(dp),INTENT(IN) :: a(3), b(3)
      real(dp) :: cross_vec(3)

        cross_vec(1) = a(2) * b(3) - a(3) * b(2)
        cross_vec(2) = a(3) * b(1) - a(1) * b(3)
        cross_vec(3) = a(1) * b(2) - a(2) * b(1)

    END FUNCTION cross
!--------------------------------------------------------------------------------------

!======================================================================================   
end module Mslab_calc

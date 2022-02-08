module Mham_w90
!======================================================================================
   use Mdef,only: dp,iu,zero,one
   use Munits,only: DPI,BohrAngstrom,HreV
   use Mlatt_utils,only: utility_recip_lattice, utility_recip_reduced
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: wann90_tb_t,ReadTB_from_w90,utility_recip_lattice
!--------------------------------------------------------------------------------------
   type :: wann90_tb_t
      integer                                    :: num_wann,nrpts,Nk      
      integer,allocatable,dimension(:)           :: ndegen
      integer,allocatable,dimension(:,:)         :: irvec
      real(dp),dimension(3,3)                    :: real_lattice,recip_lattice
      real(dp),dimension(3,3)                    :: recip_reduced
      complex(dp),allocatable,dimension(:,:,:)   :: ham_r
      complex(dp),allocatable,dimension(:,:,:,:) :: pos_r
   contains
      procedure,public  :: ReadFromW90
      procedure,public  :: Set
      procedure,public  :: get_eig
      procedure,public  :: get_ham_diag
      procedure,public  :: get_ham
      procedure,public  :: get_gradk_ham
      procedure,public  :: get_hess_ham
      procedure,public  :: get_dipole
      procedure,public  :: get_velocity
      procedure,public  :: get_emp_velocity
      procedure,public  :: get_berrycurv
      procedure,public  :: get_berrycurv_dip
      procedure,public  :: get_oam
      procedure,public  :: get_oam_dip
#if WITHHDF5
      procedure,public  :: ReadFromHDF5
      procedure,public  :: SaveToHDF5
#endif
   end type wann90_tb_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Set(me,w90)
      class(wann90_tb_t) :: me
      type(wann90_tb_t),intent(in) :: w90

      me%num_wann = w90%num_wann
      me%nrpts = w90%nrpts
      me%real_lattice = w90%real_lattice
      me%recip_lattice = w90%recip_lattice
      me%recip_reduced = w90%recip_reduced

      if(allocated(me%ndegen)) deallocate(me%ndegen)
      if(allocated(me%irvec)) deallocate(me%irvec)
      if(allocated(me%ham_r)) deallocate(me%ham_r)
      if(allocated(me%pos_r)) deallocate(me%pos_r)

      allocate(me%ndegen(me%nrpts))
      allocate(me%irvec(me%nrpts,3))
      allocate(me%ham_r(me%num_wann,me%num_wann,me%nrpts))
      allocate(me%pos_r(me%num_wann,me%num_wann,me%nrpts,3))

      me%ndegen = w90%ndegen
      me%irvec = w90%irvec
      me%ham_r = w90%ham_r
      me%pos_r = w90%pos_r

   end subroutine Set
!--------------------------------------------------------------------------------------
   function get_eig(me,kpt) result(Ek)
      use Mlinalg,only: EigValsHE
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      real(dp)            :: Ek(me%num_wann)
      complex(dp) :: Hk(me%num_wann,me%num_wann)

      Hk = me%get_ham(kpt)
      call EigValsHE(Hk,Ek)

   end function get_eig
!--------------------------------------------------------------------------------------
   function get_ham_diag(me,kpt) result(Hk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      complex(dp) :: Hk(me%num_wann,me%num_wann)
      integer :: i
      real(dp)            :: Ek(me%num_wann)

      Ek = me%get_eig(kpt)
      Hk = zero
      do i=1,me%num_wann
         Hk(i,i) = Ek(i)
      end do

   end function get_ham_diag
!--------------------------------------------------------------------------------------
   function get_ham(me,kpt) result(Hk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      complex(dp) :: Hk(me%num_wann,me%num_wann)
   
      call fourier_R_to_k(kpt, me, me%ham_r, Hk, 0)

   end function get_ham
!--------------------------------------------------------------------------------------
   function get_gradk_ham(me,kpt) result(grad_Hk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      complex(dp) :: grad_Hk(me%num_wann,me%num_wann,3)
      integer :: idir
   
      do idir=1,3
         call fourier_R_to_k(kpt, me, me%ham_r, grad_Hk(:,:,idir), idir)
      end do

   end function get_gradk_ham
!--------------------------------------------------------------------------------------
   function get_hess_ham(me,kpt) result(hess_Hk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      complex(dp) :: hess_Hk(me%num_wann,me%num_wann,3,3)
      integer :: idir,jdir
   
      do idir=1,3
         do jdir=1,3
            call fourier_D2_R_to_k(kpt, me, me%ham_r, hess_Hk(:,:,idir,jdir), idir, jdir)
         end do
      end do

   end function get_hess_ham
!--------------------------------------------------------------------------------------
   function get_dipole(me,kpt,herm_part,band_basis) result(Dk)
      use Mlinalg,only: EigHE
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      logical,intent(in),optional :: herm_part
      logical,intent(in),optional :: band_basis
      logical :: herm_part_, band_basis_
      integer :: idir
      real(dp)    :: Ek(me%num_wann)
      complex(dp) :: Hk(me%num_wann,me%num_wann),Qk(me%num_wann,me%num_wann)
      complex(dp) :: Dk(me%num_wann,me%num_wann,3)
     
      herm_part_ = .true.
      if(present(herm_part)) herm_part_ = herm_part
      band_basis_ = .false.
      if(present(band_basis)) band_basis_ = band_basis

      call fourier_R_to_k_vec(kpt, me, me%pos_r, OO_true=Dk)

      if(herm_part_) then
         do idir=1,3
            Dk(:,:,idir) = 0.5_dp * (Dk(:,:,idir) + conjg(transpose(Dk(:,:,idir))))
         end do
      end if

      if(band_basis_) then
         Hk = me%get_ham(kpt)
         call EigHE(Hk,Ek,Qk)
         do idir=1,3
            Dk(:,:,idir) = matmul(conjg(transpose(Qk)), matmul(Dk(:,:,idir), Qk))
         end do
      end if

   end function get_dipole
!--------------------------------------------------------------------------------------
   function get_emp_velocity(me,kpt) result(vk)
      use Mlinalg,only: EigHE
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      complex(dp)         :: vk(me%num_wann,me%num_wann,3)
      integer :: idir
      real(dp)    :: Ek(me%num_wann)
      complex(dp) :: Hk(me%num_wann,me%num_wann),Qk(me%num_wann,me%num_wann)

      vk = me%get_gradk_ham(kpt)
      Hk = me%get_ham(kpt)
      call EigHE(Hk,Ek,Qk)
      do idir=1,3
         vk(:,:,idir) = matmul(conjg(transpose(Qk)), matmul(vk(:,:,idir), Qk))
      end do

   end function get_emp_velocity
!--------------------------------------------------------------------------------------
   function get_velocity(me,kpt,herm_part) result(Vk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      logical,intent(in),optional :: herm_part
      logical :: herm_part_
      integer :: idir
      complex(dp)         :: vk(me%num_wann,me%num_wann,3)

      integer :: i,j
      real(dp) :: eig(me%num_wann)
      real(dp) :: del_eig(me%num_wann,3)
      complex(dp),allocatable :: HH(:,:)
      complex(dp),allocatable :: delHH(:,:,:)
      complex(dp),allocatable :: UU(:,:)
      complex(dp),allocatable :: D_h(:,:,:)
      complex(dp),allocatable :: AA(:,:,:)

      herm_part_ = .true.
      if(present(herm_part)) herm_part_ = herm_part

      allocate(HH(me%num_wann, me%num_wann))
      allocate(delHH(me%num_wann, me%num_wann, 3))
      allocate(UU(me%num_wann, me%num_wann))
      allocate(D_h(me%num_wann, me%num_wann, 3))
      allocate(AA(me%num_wann, me%num_wann, 3))

      call wham_get_eig_deleig(kpt, me, eig, del_eig, HH, delHH, UU)

      call wham_get_D_h(me%num_wann, delHH, UU, eig, D_h)

      call fourier_R_to_k_vec(kpt, me, me%pos_r, OO_true=AA)
      do i = 1, 3
         AA(:, :, i) = utility_rotate(AA(:, :, i), UU, me%num_wann)
      end do
      AA = AA + iu*D_h ! Eq.(25) WYSV06

      vk = zero
      do i=1,me%num_wann
         vk(i,i,:) = del_eig(i,:)
      end do

      do i=1,me%num_wann
         do j=1,me%num_wann
            vk(j,i,:) = vk(j,i,:) - iu * (eig(i) - eig(j)) * AA(j,i,:)
         end do
      end do


      if(herm_part_) then
         do i=1,3
            vk(:,:,i) = 0.5_dp * (vk(:,:,i) + conjg(transpose(vk(:,:,i))))
         end do
      end if

      deallocate(HH,delHH,UU,D_h,AA)

   end function get_velocity
!--------------------------------------------------------------------------------------
   function get_berrycurv(me, kpt) result(Wk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      real(dp)            :: Wk(me%num_wann,3)

      integer :: i,j,idir
      real(dp) :: eig(me%num_wann)
      real(dp) :: del_eig(me%num_wann,3)
      complex(dp),allocatable :: HH(:,:)
      complex(dp),allocatable :: delHH(:,:,:)
      complex(dp),allocatable :: UU(:,:)
      complex(dp),allocatable :: D_h(:,:,:)
      complex(dp),allocatable :: AA(:,:,:)

      allocate(HH(me%num_wann, me%num_wann))
      allocate(delHH(me%num_wann, me%num_wann, 3))
      allocate(UU(me%num_wann, me%num_wann))
      allocate(D_h(me%num_wann, me%num_wann, 3))
      allocate(AA(me%num_wann, me%num_wann, 3))

      call wham_get_eig_deleig(kpt, me, eig, del_eig, HH, delHH, UU)

      call wham_get_D_h(me%num_wann, delHH, UU, eig, D_h)

      call fourier_R_to_k_vec(kpt, me, me%pos_r, OO_true=AA)
      do idir = 1, 3
         AA(:, :, idir) = utility_rotate(AA(:, :, idir), UU, me%num_wann)
      end do
      AA = AA + iu*D_h ! Eq.(25) WYSV06

      Wk = 0.0_dp
      do i=1,me%num_wann
         do j=1,me%num_wann
            if(i == j) cycle
            Wk(i,3) = Wk(i,3) + 2.0_dp * aimag(AA(i,j,1)*AA(j,i,2))
            Wk(i,1) = Wk(i,1) + 2.0_dp * aimag(AA(i,j,2)*AA(j,i,3))
            Wk(i,2) = Wk(i,2) + 2.0_dp * aimag(AA(i,j,3)*AA(j,i,1))
         end do
      end do

   end function get_berrycurv
!--------------------------------------------------------------------------------------
   subroutine get_berrycurv_dip(me, kpt, Bhk, Bdip) 
      class(wann90_tb_t)   :: me
      real(dp),intent(in)  :: kpt(3)
      real(dp),intent(out) :: Bhk(me%num_wann,3)
      real(dp),intent(out) :: Bdip(me%num_wann,3)

      integer :: i,j,idir
      real(dp) :: ediff
      real(dp) :: eig(me%num_wann)
      complex(dp),dimension(me%num_wann,me%num_wann) :: Hk,UU
      complex(dp),dimension(me%num_wann,me%num_wann,3) :: grad_Hk,Dk

      Hk = me%get_ham(kpt)
      grad_Hk = me%get_gradk_ham(kpt)
      Dk = me%get_dipole(kpt)

      call utility_diagonalize(Hk, me%num_wann, eig, UU)

      do idir=1,3
         grad_Hk(:,:,idir) = utility_rotate(grad_Hk(:,:,idir), UU, me%num_wann)
         Dk(:,:,idir) = utility_rotate(Dk(:,:,idir), UU, me%num_wann)
      end do

      Bhk = 0.0_dp; Bdip = 0.0_dp
      do i=1,me%num_wann
         do j=1,me%num_wann
            if(i == j) cycle
            ediff = eig(i) - eig(j)
            if(abs(ediff) < 1.0e-7_dp) cycle

            Bhk(i,3) = Bhk(i,3) - aimag(grad_Hk(i,j,1) * grad_Hk(j,i,2) - grad_Hk(i,j,2) * grad_Hk(j,i,1)) / &
               ediff**2
            Bhk(i,1) = Bhk(i,1) - aimag(grad_Hk(i,j,2) * grad_Hk(j,i,3) - grad_Hk(i,j,3) * grad_Hk(j,i,2)) / &
               ediff**2
            Bhk(i,2) = Bhk(i,2) - aimag(grad_Hk(i,j,3) * grad_Hk(j,i,1) - grad_Hk(i,j,1) * grad_Hk(j,i,3)) / &
               ediff**2

            Bdip(i,3) = Bdip(i,3) + dble(grad_Hk(i,j,1) * Dk(j,i,2) - grad_Hk(i,j,2) * Dk(j,i,1)) / &
               ediff
            Bdip(i,1) = Bdip(i,1) + dble(grad_Hk(i,j,2) * Dk(j,i,3) - grad_Hk(i,j,3) * Dk(j,i,2)) / &
               ediff
            Bdip(i,2) = Bdip(i,2) + dble(grad_Hk(i,j,3) * Dk(j,i,1) - grad_Hk(i,j,1) * Dk(j,i,3)) / &
               ediff
         end do
      end do 

      Bdip = 2.0_dp * Bdip

   end subroutine get_berrycurv_dip
!--------------------------------------------------------------------------------------
   function get_oam(me, kpt) result(Lk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      real(dp)            :: Lk(me%num_wann,3)

      integer :: i,j,idir
      real(dp) :: eig(me%num_wann)
      real(dp) :: del_eig(me%num_wann,3)
      complex(dp),allocatable :: HH(:,:)
      complex(dp),allocatable :: delHH(:,:,:)
      complex(dp),allocatable :: UU(:,:)
      complex(dp),allocatable :: D_h(:,:,:)
      complex(dp),allocatable :: AA(:,:,:)

      allocate(HH(me%num_wann, me%num_wann))
      allocate(delHH(me%num_wann, me%num_wann, 3))
      allocate(UU(me%num_wann, me%num_wann))
      allocate(D_h(me%num_wann, me%num_wann, 3))
      allocate(AA(me%num_wann, me%num_wann, 3))

      call wham_get_eig_deleig(kpt, me, eig, del_eig, HH, delHH, UU)

      call wham_get_D_h(me%num_wann, delHH, UU, eig, D_h)

      call fourier_R_to_k_vec(kpt, me, me%pos_r, OO_true=AA)
      do idir = 1, 3
         AA(:, :, idir) = utility_rotate(AA(:, :, idir), UU, me%num_wann)
      end do
      AA = AA + iu*D_h ! Eq.(25) WYSV06

      Lk = 0.0_dp
      do i=1,me%num_wann
         do j=1,me%num_wann
            if(i == j) cycle
            Lk(i,3) = Lk(i,3) + 2.0_dp * (eig(i) - eig(j))*aimag(AA(i,j,1)*AA(j,i,2))
            Lk(i,1) = Lk(i,1) + 2.0_dp * (eig(i) - eig(j))*aimag(AA(i,j,2)*AA(j,i,3))
            Lk(i,2) = Lk(i,2) + 2.0_dp * (eig(i) - eig(j))*aimag(AA(i,j,3)*AA(j,i,1))
         end do
      end do

   end function get_oam
!--------------------------------------------------------------------------------------
   function get_oam_dip(me, kpt) result(Lk)
      class(wann90_tb_t)   :: me
      real(dp),intent(in)  :: kpt(3)
      real(dp)             :: Lk(me%num_wann,3)

      integer :: i,j,idir
      real(dp) :: ediff
      real(dp) :: eig(me%num_wann)
      real(dp),dimension(me%num_wann,3) :: Lk_disp,Lk_dip
      complex(dp),dimension(me%num_wann,me%num_wann) :: Hk,UU
      complex(dp),dimension(me%num_wann,me%num_wann,3) :: grad_Hk,Dk

      Hk = me%get_ham(kpt)
      grad_Hk = me%get_gradk_ham(kpt)
      Dk = me%get_dipole(kpt)

      call utility_diagonalize(Hk, me%num_wann, eig, UU)

      do idir=1,3
         grad_Hk(:,:,idir) = utility_rotate(grad_Hk(:,:,idir), UU, me%num_wann)
         Dk(:,:,idir) = utility_rotate(Dk(:,:,idir), UU, me%num_wann)
      end do

      Lk_disp = 0.0_dp; Lk_dip = 0.0_dp
      do i=1,me%num_wann
         do j=1,me%num_wann
            if(i == j) cycle
            ediff = eig(i) - eig(j)
            if(abs(ediff) < 1.0e-7_dp) cycle

            Lk_disp(i,3) = Lk_disp(i,3) - aimag(grad_Hk(i,j,1) * grad_Hk(j,i,2) - grad_Hk(i,j,2) * grad_Hk(j,i,1)) / &
               ediff
            Lk_disp(i,1) = Lk_disp(i,1) - aimag(grad_Hk(i,j,2) * grad_Hk(j,i,3) - grad_Hk(i,j,3) * grad_Hk(j,i,2)) / &
               ediff
            Lk_disp(i,2) = Lk_disp(i,2) - aimag(grad_Hk(i,j,3) * grad_Hk(j,i,1) - grad_Hk(i,j,1) * grad_Hk(j,i,3)) / &
               ediff

            Lk_dip(i,3) = Lk_dip(i,3) + dble(grad_Hk(i,j,1) * Dk(j,i,2) - grad_Hk(i,j,2) * Dk(j,i,1))
            Lk_dip(i,1) = Lk_dip(i,1) + dble(grad_Hk(i,j,2) * Dk(j,i,3) - grad_Hk(i,j,3) * Dk(j,i,2))
            Lk_dip(i,2) = Lk_dip(i,2) + dble(grad_Hk(i,j,3) * Dk(j,i,1) - grad_Hk(i,j,1) * Dk(j,i,3))
         end do
      end do 

      Lk_dip = 2.0_dp * Lk_dip

      Lk = Lk_disp + Lk_dip

   end function get_oam_dip
!--------------------------------------------------------------------------------------
   subroutine ReadFromW90(me,fname)
      class(wann90_tb_t)  :: me
      character(len=*),intent(in) :: fname

      call ReadTB_from_w90(fname,me%real_lattice,me%num_wann,me%nrpts,me%ndegen,&
         me%irvec,me%ham_r,me%pos_r)

      call utility_recip_lattice(me%real_lattice, me%recip_lattice)

      call utility_recip_reduced(me%recip_lattice, me%recip_reduced)

   end subroutine ReadFromW90
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine SaveToHDF5(me,fname,atomic_units)
      use Mhdf5_utils
      class(wann90_tb_t)  :: me
      character(len=*),intent(in) :: fname
      logical,intent(in),optional :: atomic_units
      integer :: au_flag = 1
      integer(HID_t) :: file_id
      real(dp),allocatable :: d_ham_r(:,:,:)
      real(dp),allocatable :: d_pos_r(:,:,:,:)

      if(present(atomic_units)) then
         if(atomic_units) then
            au_flag = 1
         else
            au_flag = 0
         end if
      end if

      call hdf_open_file(file_id, trim(fname), STATUS='NEW')

      call hdf_write_attribute(file_id,'','atomic_units', au_flag)
      call hdf_write_attribute(file_id,'','num_wann', me%num_wann)
      call hdf_write_attribute(file_id,'','nrpts', me%nrpts)

      call hdf_write_dataset(file_id,'ndegen',me%ndegen)
      call hdf_write_dataset(file_id,'real_lattice',me%real_lattice)
      call hdf_write_dataset(file_id,'irvec',me%irvec)

      allocate(d_ham_r(me%num_wann,me%num_wann,me%nrpts))

      d_ham_r = dble(me%ham_r)
      call hdf_write_dataset(file_id,'ham_r_real',d_ham_r)
      d_ham_r = aimag(me%ham_r)
      call hdf_write_dataset(file_id,'ham_r_imag',d_ham_r)

      deallocate(d_ham_r)

      if(allocated(me%pos_r)) then
         call hdf_write_attribute(file_id,'','pos_stored', 1)

         allocate(d_pos_r(me%num_wann,me%num_wann,me%nrpts,3))

         d_pos_r = dble(me%pos_r)
         call hdf_write_dataset(file_id,'pos_r_real',d_pos_r)
         d_pos_r = aimag(me%pos_r)
         call hdf_write_dataset(file_id,'pos_r_imag',d_pos_r)

         deallocate(d_pos_r)
      else
         call hdf_write_attribute(file_id,'','pos_stored', 0)
      end if

      call hdf_close_file(file_id)

   end subroutine SaveToHDF5
#endif
!--------------------------------------------------------------------------------------
#if WITHHDF5
   subroutine ReadFromHDF5(me,fname)
      use Mhdf5_utils
      class(wann90_tb_t)  :: me
      character(len=*),intent(in) :: fname
      integer(HID_t) :: file_id
      integer :: atomic_units,pos_stored
      real(dp),allocatable :: d_ham_r(:,:,:)
      real(dp),allocatable :: d_pos_r(:,:,:,:)

      call hdf_open_file(file_id, trim(fname), STATUS='OLD', ACTION='READ')

      call hdf_read_attribute(file_id,'','atomic_units', atomic_units)
      call hdf_read_attribute(file_id,'','num_wann', me%num_wann)
      call hdf_read_attribute(file_id,'','nrpts', me%nrpts)

      allocate(me%ndegen(me%nrpts))
      call hdf_read_dataset(file_id,'ndegen',me%ndegen)
      call hdf_read_dataset(file_id,'real_lattice',me%real_lattice)

      if(atomic_units == 0) me%recip_lattice = me%recip_lattice / BohrAngstrom

      call utility_recip_lattice(me%real_lattice, me%recip_lattice)
      call utility_recip_reduced(me%recip_lattice, me%recip_reduced)

      allocate(me%irvec(me%nrpts,3))
      call hdf_read_dataset(file_id,'irvec',me%irvec)

      allocate(d_ham_r(me%num_wann,me%num_wann,me%nrpts))
      allocate(me%ham_r(me%num_wann,me%num_wann,me%nrpts))

      call hdf_read_dataset(file_id,'ham_r_real',d_ham_r)
      me%ham_r = d_ham_r
      call hdf_read_dataset(file_id,'ham_r_imag',d_ham_r)
      me%ham_r = me%ham_r + iu * d_ham_r

      if(atomic_units == 0) me%ham_r = me%ham_r / HreV

      deallocate(d_ham_r)

      call hdf_read_attribute(file_id,'','pos_stored', pos_stored)

      if(pos_stored == 1) then
         allocate(d_pos_r(me%num_wann,me%num_wann,me%nrpts,3))
         allocate(me%pos_r(me%num_wann,me%num_wann,me%nrpts,3))

         call hdf_read_dataset(file_id,'pos_r_real',d_pos_r)
         me%pos_r = d_pos_r
         call hdf_read_dataset(file_id,'pos_r_imag',d_pos_r)
         me%pos_r = me%pos_r + iu * d_pos_r

         deallocate(d_pos_r)

         if(atomic_units == 0) me%pos_r = me%pos_r / BohrAngstrom
      end if

      call hdf_close_file(file_id)

   end subroutine ReadFromHDF5
#endif
!--------------------------------------------------------------------------------------
   subroutine ReadTB_from_w90(fname,real_lattice,num_wann,nrpts,ndegen,irvec,ham_r,pos_r)
      character(len=*),intent(in) :: fname
      real(dp),dimension(3,3),intent(out) :: real_lattice
      integer,intent(out) :: num_wann,nrpts
      integer,allocatable,intent(out) :: ndegen(:),irvec(:,:)
      complex(dp),allocatable,intent(out) :: ham_r(:,:,:)
      complex(dp),allocatable,intent(out) :: pos_r(:,:,:,:)
      integer :: rst,qst,i,j,irpt,ndx1,ndx2
      integer :: file_unit
      real(dp) :: a,b,pos_real(3),pos_imag(3)
      complex(dp) :: pos(3)
      
      open(newunit=file_unit, file=trim(fname), status='old', action='read')

      read(file_unit, *)  ! Date and time
      !
      ! lattice vectors
      !
      read(file_unit, *) real_lattice(1, :) !a_1
      read(file_unit, *) real_lattice(2, :) !a_2
      read(file_unit, *) real_lattice(3, :) !a_3
      
      ! convert to atomic units
      real_lattice = real_lattice / BohrAngstrom

      read(file_unit, *) num_wann
      read(file_unit, *) nrpts
      rst=mod(nrpts,15)
      qst=int(nrpts/15)
      allocate(ndegen(nrpts),irvec(nrpts,3))

      ! read WS degeneracies
      do i=1,qst
         read(file_unit,*) (ndegen(j+(i-1)*15),j=1,15)
      end do
      if(rst.ne.0) read(file_unit,*) (ndegen(j+qst*15),j=1,rst)

      ! read real-space Hamiltonian (no spinup-spindw hybridizations assumed)
      allocate(ham_r(num_wann,num_wann,nrpts))
      do irpt = 1,nrpts
         read(file_unit,*) irvec(irpt,1),irvec(irpt,2),irvec(irpt,3)
         do i=1,num_wann
            do j=1,num_wann
               read(file_unit,*) ndx1, ndx2, a, b
               ham_r(j,i,irpt) = (a + iu*b) / HreV
            end do
         end do
      end do

      ! read dipole matrix elements 
      allocate(pos_r(num_wann,num_wann,nrpts,3))
      do irpt = 1,nrpts
         read(file_unit,*) irvec(irpt,1),irvec(irpt,2),irvec(irpt,3)
         do i=1,num_wann
            do j=1,num_wann
               ! read(file_unit,*) ndx1, ndx2, pos_real, pos_imag
               ! pos_r(j,i,irpt,:) = (pos_real + iu*pos_imag) / BohrAngstrom
               read(file_unit,'(2I5,3x,6(E15.8,1x))') ndx1, ndx2, pos
               pos_r(j,i,irpt,:) = pos / BohrAngstrom
            end do
         end do
      end do

      close(file_unit)

   end subroutine ReadTB_from_w90
!--------------------------------------------------------------------------------------
   function utility_rotate(mat, rot, dim)
      !==========================================================!
      !                                                           !
      !! Rotates the dim x dim matrix 'mat' according to
      !! (rot)^dagger.mat.rot, where 'rot' is a unitary matrix
      !                                                           !
      !===========================================================!
      integer          :: dim
      complex(kind=dp) :: utility_rotate(dim, dim)
      complex(kind=dp) :: mat(dim, dim)
      complex(kind=dp) :: rot(dim, dim)

      utility_rotate = matmul(matmul(transpose(conjg(rot)), mat), rot)

   end function utility_rotate
!--------------------------------------------------------------------------------------
   subroutine utility_diagonalize(mat, dim, eig, rot)
      !============================================================!
      !                                                            !
      !! Diagonalize the dim x dim  hermitian matrix 'mat' and
      !! return the eigenvalues 'eig' and the unitary rotation 'rot'
      !                                                            !
      !============================================================!
      use,intrinsic::iso_fortran_env,only: output_unit, error_unit

      integer, intent(in)           :: dim
      complex(kind=dp), intent(in)  :: mat(dim, dim)
      real(kind=dp), intent(out)    :: eig(dim)
      complex(kind=dp), intent(out) :: rot(dim, dim)

      complex(kind=dp)   :: mat_pack((dim*(dim + 1))/2), cwork(2*dim)
      real(kind=dp)      :: rwork(7*dim)
      integer            :: i, j, info, nfound, iwork(5*dim), ifail(dim)

      do j = 1, dim
         do i = 1, j
            mat_pack(i + ((j - 1)*j)/2) = mat(i, j)
         end do
      end do
      rot = zero; eig = 0.0_dp; cwork = zero; rwork = 0.0_dp; iwork = 0
      call ZHPEVX('V', 'A', 'U', dim, mat_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
         nfound, eig(1), rot, dim, cwork, rwork, iwork, ifail, info)
      if (info < 0) then
         write (output_unit, '(a,i3,a)') 'THE ', -info, &
         ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
         write(error_unit,*) 'Error in utility_diagonalize'
      end if
      if (info > 0) then
         write (output_unit, '(i3,a)') info, ' EIGENVECTORS FAILED TO CONVERGE'
         write(error_unit,*) 'Error in utility_diagonalize'
      end if

  end subroutine utility_diagonalize
!--------------------------------------------------------------------------------------
   function utility_rotate_diag(mat, rot, dim)
      !===========================================================!
      !                                                           !
      !! Rotates the dim x dim matrix 'mat' according to
      !! (rot)^dagger.mat.rot, where 'rot' is a unitary matrix.
      !! Computes only the diagonal elements of rotated matrix.
      !                                                           !
      !===========================================================!
      integer          :: dim
      complex(kind=dp) :: utility_rotate_diag(dim)
      complex(kind=dp) :: mat(dim, dim)
      complex(kind=dp) :: rot(dim, dim)
      complex(kind=dp) :: tmp(dim, dim)

      call utility_zgemm_new(rot, mat, tmp, 'C', 'N')
      utility_rotate_diag = utility_matmul_diag(tmp, rot, dim)

   end function utility_rotate_diag
!--------------------------------------------------------------------------------------
 subroutine utility_zgemm_new(a, b, c, transa_opt, transb_opt)
   !=============================================================!
   !                                                             !
   ! Return matrix product of complex matrices a and b:          !
   !                                                             !
   !                       C = Op(A) Op(B)                       !
   !                                                             !
   ! transa = 'N'  ==> Op(A) = A                                 !
   ! transa = 'T'  ==> Op(A) = transpose(A)                      !
   ! transa = 'C'  ==> Op(A) = congj(transpose(A))               !
   !                                                             !
   ! similarly for B                                             !
   !                                                             !
   ! Due to the use of assumed shape arrays, this routine is a   !
   ! safer and more general replacement for the above routine    !
   ! utility_zgemm. Consider removing utility_zgemm and using    !
   ! utility_zgemm_new throughout.                               !
   !                                                             !
   !=============================================================!
      complex(kind=dp), intent(in)            :: a(:, :)
      complex(kind=dp), intent(in)            :: b(:, :)
      complex(kind=dp), intent(out)           :: c(:, :)
      character(len=1), intent(in), optional  :: transa_opt
      character(len=1), intent(in), optional  :: transb_opt

      integer          :: m, n, k
      character(len=1) :: transa, transb

      transa = 'N'
      transb = 'N'
      if (present(transa_opt)) transa = transa_opt
      if (present(transb_opt)) transb = transb_opt

      ! m ... number of rows in Op(A) and C
      ! n ... number of columns in Op(B) and C
      ! k ... number of columns in Op(A) resp. rows in Op(B)
      m = size(c, 1)
      n = size(c, 2)

      if (transa /= 'N') then
         k = size(a, 1)
      else
         k = size(a, 2)
      end if

      call zgemm(transa, transb, m, n, k, one, a, size(a, 1), b, size(b, 1), zero, c, m)

   end subroutine utility_zgemm_new
!--------------------------------------------------------------------------------------
   function utility_matmul_diag(mat1, mat2, dim)
      !===========================================================!
      !                                                           !
      !! Computes the diagonal elements of the matrix mat1.mat2
      !                                                           !
      !===========================================================!
      integer          :: dim
      complex(kind=dp) :: utility_matmul_diag(dim)
      complex(kind=dp) :: mat1(dim, dim)
      complex(kind=dp) :: mat2(dim, dim)

      integer i, j

      utility_matmul_diag = zero
      do i = 1, dim
         do j = 1, dim
            utility_matmul_diag(i) = utility_matmul_diag(i) + mat1(i, j)*mat2(j, i)
         end do
      end do

   end function utility_matmul_diag
!--------------------------------------------------------------------------------------
   subroutine wham_get_D_h(num_wann, delHH, UU, eig, D_h)
      !=========================================!
      !                                         !
      !! Compute D^H_a=UU^dag.del_a UU (a=x,y,z)
      !! using Eq.(24) of WYSV06
      !                                         !
      !=========================================!

      ! TO DO: Implement version where energy denominators only connect
      !        occupied and empty states. In this case probably do not need
      !        to worry about avoiding small energy denominators

      ! Arguments
      !
      integer,intent(in) :: num_wann
      complex(kind=dp), dimension(:, :, :), intent(in)  :: delHH
      complex(kind=dp), dimension(:, :), intent(in)    :: UU
      real(kind=dp), dimension(:), intent(in)    :: eig
      complex(kind=dp), dimension(:, :, :), intent(out) :: D_h

      complex(kind=dp), allocatable :: delHH_bar_i(:, :)
      integer                       :: n, m, i

      allocate (delHH_bar_i(num_wann, num_wann))
      D_h = zero
      do i = 1, 3
         delHH_bar_i(:, :) = utility_rotate(delHH(:, :, i), UU, num_wann)
         do m = 1, num_wann
            do n = 1, num_wann
               if (n == m .or. abs(eig(m) - eig(n)) < 1.0e-7_dp) cycle
               D_h(n, m, i) = delHH_bar_i(n, m)/(eig(m) - eig(n))
            end do
         end do
      enddo

      deallocate(delHH_bar_i)

   end subroutine wham_get_D_h
!--------------------------------------------------------------------------------------
  subroutine wham_get_eig_deleig(kpt, w90, eig, del_eig, HH, delHH, UU)
    !! Given a k point, this function returns eigenvalues E and
    !! derivatives of the eigenvalues dE/dk_a, using wham_get_deleig_a

    real(kind=dp), dimension(3), intent(in)         :: kpt
    !! the three coordinates of the k point vector (in relative coordinates)
    type(wann90_tb_t), intent(in)                   :: w90
    real(kind=dp), intent(out)                      :: eig(w90%num_wann)
    !! the calculated eigenvalues at kpt
    real(kind=dp), intent(out)                      :: del_eig(w90%num_wann, 3)
    !! the calculated derivatives of the eigenvalues at kpt [first component: band; second component: 1,2,3
    !! for the derivatives along the three k directions]
    complex(kind=dp), dimension(:, :), intent(out)   :: HH
    !! the Hamiltonian matrix at kpt
    complex(kind=dp), dimension(:, :, :), intent(out) :: delHH
    !! the delHH matrix (derivative of H) at kpt
    complex(kind=dp), dimension(:, :), intent(out)   :: UU
    !! the rotation matrix that gives the eigenvectors of HH

    call fourier_R_to_k(kpt, w90, w90%ham_r, HH, 0)
    call utility_diagonalize(HH, w90%num_wann, eig, UU)
    call fourier_R_to_k(kpt,  w90, w90%ham_r, delHH(:, :, 1), 1)
    call fourier_R_to_k(kpt,  w90, w90%ham_r, delHH(:, :, 2), 2)
    call fourier_R_to_k(kpt,  w90, w90%ham_r, delHH(:, :, 3), 3)
    call wham_get_deleig_a(w90%num_wann, del_eig(:, 1), w90, eig, delHH(:, :, 1), UU)
    call wham_get_deleig_a(w90%num_wann, del_eig(:, 2), w90, eig, delHH(:, :, 2), UU)
    call wham_get_deleig_a(w90%num_wann, del_eig(:, 3), w90, eig, delHH(:, :, 3), UU)

  end subroutine wham_get_eig_deleig
!--------------------------------------------------------------------------------------
   subroutine wham_get_deleig_a(num_wann, deleig_a, w90, eig, delHH_a, UU, use_degen_pert, degen_thr)
      !==========================!
      !                          !
      !! Band derivatives dE/dk_a
      !                          !
      !==========================!
      ! Arguments
      !
      integer,intent(in) :: num_wann
      real(kind=dp), intent(out) :: deleig_a(num_wann)
      type(wann90_tb_t),intent(in) :: w90
      real(kind=dp), intent(in)  :: eig(num_wann)
      complex(kind=dp), dimension(:, :), intent(in)  :: delHH_a
      complex(kind=dp), dimension(:, :), intent(in)  :: UU
      logical,intent(in),optional  :: use_degen_pert
      real(dp),intent(in),optional :: degen_thr

      ! Misc/Dummy
      !
      logical :: use_degen_pert_=.false.
      real(dp) :: degen_thr_ = 1.0e-4_dp
      integer                       :: i, degen_min, degen_max, dim
      real(kind=dp)                 :: diff
      complex(kind=dp), allocatable :: delHH_bar_a(:, :), U_deg(:, :)

      if(present(use_degen_pert)) use_degen_pert_ = use_degen_pert
      if(present(degen_thr)) degen_thr_ = degen_thr

      allocate (delHH_bar_a(num_wann, num_wann))
      allocate (U_deg(num_wann, num_wann))

      if (use_degen_pert_) then

         delHH_bar_a = utility_rotate(delHH_a, UU, num_wann)

         ! Assuming that the energy eigenvalues are stored in eig(:) in
         ! increasing order (diff >= 0)

         i = 0
         do
            i = i + 1
            if (i > num_wann) exit
            if (i + 1 <= num_wann) then
               diff = eig(i + 1) - eig(i)
            else
               !
               ! i-th is the highest band, and it is non-degenerate
               !
               diff = degen_thr_ + 1.0_dp
            end if
            if (diff < degen_thr_) then
               !
               ! Bands i and i+1 are degenerate
               !
               degen_min = i
               degen_max = degen_min + 1
               !
               ! See if any higher bands are in the same degenerate group
               !
               do
                  if (degen_max + 1 > num_wann) exit
                  diff = eig(degen_max + 1) - eig(degen_max)
                  if (diff < degen_thr_) then
                     degen_max = degen_max + 1
                  else
                     exit
                  end if
               end do
               !
               ! Bands from degen_min to degen_max are degenerate. Diagonalize
               ! the submatrix in Eq.(31) YWVS07 over this degenerate subspace.
               ! The eigenvalues are the band gradients
               !
               !
               dim = degen_max - degen_min + 1
               call utility_diagonalize(delHH_bar_a(degen_min:degen_max, &
                 degen_min:degen_max), dim, &
               deleig_a(degen_min:degen_max), U_deg(1:dim, 1:dim))
               !
               ! Scanned bands up to degen_max
               !
               i = degen_max
            else
               !
               ! Use non-degenerate form [Eq.(27) YWVS07] for current (i-th) band
               !
               deleig_a(i) = real(delHH_bar_a(i, i), dp)
            end if
         end do

      else

         ! Use non-degenerate form for all bands
         !
         deleig_a(:) = real(utility_rotate_diag(delHH_a(:, :), UU, num_wann), dp)

      end if

   end subroutine wham_get_deleig_a
!--------------------------------------------------------------------------------------
   subroutine get_crvec(w90,crvec)
      type(wann90_tb_t),intent(in)     :: w90
      real(dp),allocatable,intent(out) :: crvec(:,:)
      integer :: ir

      if(allocated(crvec)) deallocate(crvec)
      allocate(crvec(3,w90%nrpts))
      do ir = 1, w90%nrpts
        ! Note that 'real_lattice' stores the lattice vectors as *rows*
        crvec(:, ir) = matmul(transpose(w90%real_lattice), w90%irvec(ir, :))
      end do

   end subroutine
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k(kpt, w90, OO_R, OO, alpha)
      !=========================================================!
      !                                                         !
      !! For alpha=0:
      !! O_ij(R) --> O_ij(k) = sum_R e^{+ik.R}*O_ij(R)
      !!
      !! For alpha=1,2,3:
      !! sum_R [cmplx_i*R_alpha*e^{+ik.R}*O_ij(R)]
      !! where R_alpha is a Cartesian component of R
      !! ***REMOVE EVENTUALLY*** (replace with pw90common_fourier_R_to_k_new)
      !                                                         !
      !=========================================================!

      ! Arguments
      !
      real(kind=dp)                                     :: kpt(3)
      type(wann90_tb_t),intent(in)                      :: w90
      complex(kind=dp), dimension(:, :, :), intent(in)  :: OO_R
      complex(kind=dp), dimension(:, :), intent(inout)  :: OO
      integer                                           :: alpha

      integer          :: ir, i, j, ideg
      real(kind=dp)    :: rdotk
      complex(kind=dp) :: phase_fac
      real(dp),allocatable :: crvec(:,:)

      if(alpha > 0) call get_crvec(w90, crvec)

      OO(:, :) = zero
      do ir = 1, w90%nrpts
         rdotk = DPI*dot_product(kpt(:), w90%irvec(ir, :))
         phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(w90%ndegen(ir), dp)
         if (alpha == 0) then
            OO(:, :) = OO(:, :) + phase_fac*OO_R(:, :, ir)
         elseif (alpha == 1 .or. alpha == 2 .or. alpha == 3) then
            OO(:, :) = OO(:, :) + &
               iu*crvec(alpha, ir)*phase_fac*OO_R(:, :, ir)
         else
            stop 'wrong value of alpha in fourier_R_to_k'
         end if

      end do

      if(allocated(crvec)) deallocate(crvec)

   end subroutine fourier_R_to_k
!--------------------------------------------------------------------------------------
   subroutine fourier_D2_R_to_k(kpt, w90, OO_R, OO, a, b)
      !=========================================================!
      !                                                         !
      !! sum_R [- R_a * R_b * e^{+ik.R}*O_ij(R)]
      !! where R_a is a Cartesian component of R
      !                                                         !
      !=========================================================!

      ! Arguments
      !
      real(kind=dp)                                     :: kpt(3)
      type(wann90_tb_t),intent(in)                      :: w90
      complex(kind=dp), dimension(:, :, :), intent(in)  :: OO_R
      complex(kind=dp), dimension(:, :), intent(inout)  :: OO
      integer,intent(in)                                :: a,b

      integer          :: ir, i, j, ideg
      real(kind=dp)    :: rdotk
      complex(kind=dp) :: phase_fac
      real(dp),allocatable :: crvec(:,:)

      call get_crvec(w90, crvec)

      OO(:, :) = zero
      do ir = 1, w90%nrpts
         rdotk = DPI*dot_product(kpt(:), w90%irvec(ir, :))
         phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(w90%ndegen(ir), dp)
         OO(:, :) = OO(:, :) - &
            crvec(a, ir)*crvec(b, ir)*phase_fac*OO_R(:, :, ir)

      end do

      if(allocated(crvec)) deallocate(crvec)

   end subroutine fourier_D2_R_to_k
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k_vec(kpt, w90, OO_R, OO_true, OO_pseudo)
      !====================================================================!
      !                                                                    !
      !! For OO_true (true vector):
      !! $${\vec O}_{ij}(k) = \sum_R e^{+ik.R} {\vec O}_{ij}(R)$$
      !                                                                    !
      !====================================================================!

      ! Arguments
      !
      real(kind=dp)                                     :: kpt(3)
      type(wann90_tb_t),intent(in)                      :: w90
      complex(kind=dp), dimension(:, :, :, :), intent(in)  :: OO_R
      complex(kind=dp), optional, dimension(:, :, :), intent(inout)   :: OO_true
      complex(kind=dp), optional, dimension(:, :, :), intent(inout)   :: OO_pseudo

      integer          :: ir, i, j, ideg
      real(kind=dp)    :: rdotk
      complex(kind=dp) :: phase_fac
      real(dp),allocatable :: crvec(:,:)

      if (present(OO_pseudo)) call get_crvec(w90, crvec)

      if (present(OO_true)) OO_true = zero
      if (present(OO_pseudo)) OO_pseudo = zero

      do ir = 1, w90%nrpts
         rdotk = DPI*dot_product(kpt(:), w90%irvec(ir, :))
         phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(w90%ndegen(ir), dp)

         if (present(OO_true)) then
            OO_true(:, :, 1) = OO_true(:, :, 1) + phase_fac*OO_R(:, :, ir, 1)
            OO_true(:, :, 2) = OO_true(:, :, 2) + phase_fac*OO_R(:, :, ir, 2)
            OO_true(:, :, 3) = OO_true(:, :, 3) + phase_fac*OO_R(:, :, ir, 3)
         end if
         if (present(OO_pseudo)) then
            OO_pseudo(:, :, 1) = OO_pseudo(:, :, 1) &
               + iu*crvec(2, ir)*phase_fac*OO_R(:, :, ir, 3) &
               - iu*crvec(3, ir)*phase_fac*OO_R(:, :, ir, 2)
            OO_pseudo(:, :, 2) = OO_pseudo(:, :, 2) &
               + iu*crvec(3, ir)*phase_fac*OO_R(:, :, ir, 1) &
               - iu*crvec(1, ir)*phase_fac*OO_R(:, :, ir, 3)
            OO_pseudo(:, :, 3) = OO_pseudo(:, :, 3) &
               + iu*crvec(1, ir)*phase_fac*OO_R(:, :, ir, 2) &
               - iu*crvec(2, ir)*phase_fac*OO_R(:, :, ir, 1)
         end if

      end do

      if(allocated(crvec)) deallocate(crvec)

   end subroutine fourier_R_to_k_vec
!--------------------------------------------------------------------------------------

!======================================================================================
end module Mham_w90
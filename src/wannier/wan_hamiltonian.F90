module wan_hamiltonian
!! Provides tools for reading/writing the Wannier Hamiltonian and for computing 
!! band structures, observables, and Berry-phase properties.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
   use scitools_linalg,only: get_large_size,util_zgemm,util_matmul,util_rotate,util_rotate_cc
   use wan_utils,only: utility_recip_lattice, utility_recip_reduced, utility_diagonalize, &
      utility_rotate_diag, utility_matmul_diag
   use wan_read_xyz,only: ReadXYZ
   implicit none
   include '../units_inc.f90'
   include '../formats.h'
!--------------------------------------------------------------------------------------
   private
   public :: wann90_tb_t,ReadTB_from_w90,utility_recip_lattice
!--------------------------------------------------------------------------------------
   integer,parameter :: field_mode_positions=0,field_mode_dipole=1,field_mode_berry=2

   integer,parameter :: blocksize=32
!--------------------------------------------------------------------------------------
   type :: wann90_tb_t
      !! Wannier Hamiltonian class. Reads/writes the Wannier Hamiltonian from/to file,
      !! calculates the Hamiltonian and Berry-phase properties
      logical                                    :: coords_present=.false.
      integer                                    :: num_wann,nrpts
      integer,allocatable,dimension(:)           :: ndegen
      integer,allocatable,dimension(:,:)         :: irvec
      real(dp),dimension(3,3)                    :: real_lattice,recip_lattice
      real(dp),dimension(3,3)                    :: recip_reduced
      real(dp),allocatable,dimension(:,:)        :: coords
      real(dp),allocatable,dimension(:,:)        :: crvec
      complex(dp),allocatable,dimension(:,:,:)   :: ham_r
      complex(dp),allocatable,dimension(:,:,:)   :: OO_r
      complex(dp),allocatable,dimension(:,:,:,:) :: pos_r
      ! .. internal parameters for tuning the calculation ..
      logical  :: use_degen_pert=.false.,force_herm=.true.,force_antiherm=.true.
      real(dp) :: degen_thresh=1.0e-5_dp
   contains
      procedure,public  :: ReadFromW90
      procedure,public  :: SaveToW90
      procedure,public  :: SetExpertParams
      procedure,public  :: Set
      procedure,public  :: get_kcart
      procedure,public  :: get_kreduced
      procedure,public  :: get_eig
      procedure,public  :: get_ham_diag
      procedure,public  :: get_ham
      procedure,public  :: get_ham_field 
      procedure,public  :: get_ham_Peierls_Dipole
      procedure,public  :: get_ham_elpot
      procedure,public  :: get_gradk_ham
      procedure,public  :: get_hess_ham
      procedure,public  :: get_dipole
      procedure,public  :: get_velocity
      procedure,public  :: get_berry_connection
      procedure,public  :: get_emp_velocity
      procedure,public  :: get_berrycurv
      procedure,public  :: get_berrycurv_dip
      procedure,public  :: get_oam
      procedure,public  :: get_oam_dip
      procedure,public  :: get_metric
      procedure,public  :: get_spin_velocity
      procedure,public  :: get_spin_berrycurv
      procedure,public  :: get_spin_berrycurv_dip      
      procedure,public  :: Clean   
#if WITHHDF5
      procedure,public  :: ReadFromHDF5
      procedure,public  :: SaveToHDF5
#endif
   end type wann90_tb_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Clean(me)
      class(wann90_tb_t) :: me

      if(allocated(me%ndegen)) deallocate(me%ndegen)
      if(allocated(me%irvec)) deallocate(me%irvec)
      if(allocated(me%ham_r)) deallocate(me%ham_r)
      if(allocated(me%OO_r)) deallocate(me%OO_r)   
      if(allocated(me%pos_r)) deallocate(me%pos_r)      
      if(allocated(me%coords)) deallocate(me%coords)      

   end subroutine Clean
!--------------------------------------------------------------------------------------
   subroutine SetExpertParams(me,use_degen_pert,degen_thresh,force_herm,force_antiherm)
      class(wann90_tb_t)  :: me
      logical,intent(in),optional  :: use_degen_pert
      real(dp),intent(in),optional :: degen_thresh
      logical,intent(in),optional  :: force_herm
      logical,intent(in),optional  :: force_antiherm

      if(present(use_degen_pert)) me%use_degen_pert = use_degen_pert
      if(present(degen_thresh)) me%degen_thresh = degen_thresh
      if(present(force_herm)) me%force_herm = force_herm
      if(present(force_antiherm)) me%force_antiherm = force_antiherm

   end subroutine SetExpertParams
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
      if(allocated(me%OO_r)) deallocate(me%OO_r)
      if(allocated(me%pos_r)) deallocate(me%pos_r)
      if(allocated(me%coords)) deallocate(me%coords)
      if(allocated(me%crvec)) deallocate(me%crvec)

      allocate(me%ndegen(me%nrpts))
      allocate(me%irvec(me%nrpts,3))
      allocate(me%ham_r(me%num_wann,me%num_wann,me%nrpts))
      allocate(me%OO_r(me%num_wann,me%num_wann,me%nrpts))
      allocate(me%pos_r(me%num_wann,me%num_wann,3,me%nrpts))
      allocate(me%coords(me%num_wann,3))
      allocate(me%crvec(me%nrpts,3))

      me%ndegen = w90%ndegen
      me%irvec = w90%irvec
      me%ham_r = w90%ham_r
      me%pos_r = w90%pos_r
      me%coords = w90%coords
      me%coords_present = w90%coords_present
      me%crvec = w90%crvec

   end subroutine Set
!--------------------------------------------------------------------------------------
   function get_kcart(me,kpt) result(kvec)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(:)
      real(dp) :: kvec(size(kpt, dim=1))
      integer :: ndim

      ndim = size(kpt, dim=1)
      select case(ndim)
      case(1)
         kvec(1) = me%recip_lattice(1,1) * kpt(1)
      case(2)
         kvec(1:2) = me%recip_lattice(1,1:2) * kpt(1) + me%recip_lattice(2,1:2) * kpt(2)
      case(3)
         kvec(1:3) = me%recip_lattice(1,1:3) * kpt(1) + me%recip_lattice(2,1:3) * kpt(2) &
            + me%recip_lattice(3,1:3) * kpt(3) 
      case default
         kvec = 0.0_dp
      end select

   end function get_kcart
!--------------------------------------------------------------------------------------
   function get_kreduced(me,kvec) result(kpt)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kvec(:)
      real(dp) :: kpt(size(kvec, dim=1))
      integer :: ndim

      ndim = size(kpt, dim=1)
      select case(ndim)
      case(1)
         kpt(1) = kpt(1) / me%recip_lattice(1,1)
      case(2)
         kpt(1:2) = matmul(me%recip_reduced(1:2,1:2), kvec(1:2))
      case(3)
         kpt(1:3) = matmul(me%recip_reduced(1:3,1:3), kvec(1:3))
      case default
         kpt = 0.0_dp
      end select
     
   end function get_kreduced
!--------------------------------------------------------------------------------------
   function get_eig(me,kpt) result(Ek)
      use scitools_linalg,only: EigValsHE
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
   
      call fourier_R_to_k(kpt, me, me%ham_r, Hk)

   end function get_ham
!--------------------------------------------------------------------------------------
   function get_ham_field(me,kpt,Ef,field_mode) result(Hk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      real(dp),intent(in) :: Ef(3)
      integer,intent(in)  :: field_mode
      complex(dp) :: Hk(me%num_wann,me%num_wann)
      integer :: i,idir
      complex(dp),dimension(me%num_wann,me%num_wann,3) :: Dk,AA
   
      Hk = me%get_ham(kpt)

      select case(field_mode)
      case(field_mode_positions)
         do i=1,me%num_wann
            Hk(i,i) = Hk(i,i) - Ef(1) * me%coords(i,1)
            Hk(i,i) = Hk(i,i) - Ef(2) * me%coords(i,2)
            Hk(i,i) = Hk(i,i) - Ef(3) * me%coords(i,3)
         end do
      case(field_mode_dipole)
         Dk = me%get_dipole(kpt,band_basis=.false.)
         do idir=1,3
            Hk(:,:) = Hk(:,:) - Ef(idir) * Dk(:,:,idir)
         end do
      case(field_mode_berry)
         AA = me%get_berry_connection(kpt)
         do idir=1,3
            Hk(:,:) = Hk(:,:) - Ef(idir) * AA(:,:,idir)
         end do         
      end select 


   end function get_ham_field
!--------------------------------------------------------------------------------------
   function get_ham_Peierls_Dipole(me,kpt,Av,Ef,Peierls_only) result(Hk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      real(dp),intent(in) :: Av(3)
      real(dp),intent(in) :: Ef(3)
      logical,intent(in),optional :: Peierls_only
      logical :: peierls_
      complex(dp) :: Hk(me%num_wann,me%num_wann)
      integer :: idir,ir
      real(dp) :: Ared(3),kA(3)
      complex(dp) :: za

      peierls_ = .false.
      if(present(Peierls_only)) peierls_ = Peierls_only

      me%OO_r = me%ham_r

      if(.not.peierls_) then
         do ir=1,me%nrpts
            do idir=1,3
               za = -Ef(idir)
               ! me%OO_r(:,:,ir) = me%OO_r(:,:,ir) - Ef(idir) * me%pos_r(:,:,idir,ir)
               call ZAXPY(me%num_wann**2,za,me%pos_r(1,1,idir,ir),1,me%OO_r(1,1,ir),1)
            end do
         end do
      end if

      Ared = me%get_kreduced(Av)
      kA = kpt - Ared

      call fourier_R_to_k(kA, me, me%OO_r, Hk)

   end function get_ham_Peierls_Dipole
!--------------------------------------------------------------------------------------
   function get_ham_elpot(me,kpt,elpot) result(Hk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      interface 
         function elpot(z) result(v)
            use scitools_def,only: dp
            real(dp),intent(in) :: z
            real(dp) :: v
         end function elpot
      end interface
      complex(dp) :: Hk(me%num_wann,me%num_wann)
      integer :: i,idir
      complex(dp),dimension(me%num_wann,me%num_wann,3) :: Dk,AA
   
      Hk = me%get_ham(kpt)
      do i=1,me%num_wann
         Hk(i,i) = Hk(i,i) - elpot(me%coords(i,3))
      end do

   end function get_ham_elpot
!--------------------------------------------------------------------------------------
   function get_gradk_ham(me,kpt) result(grad_Hk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      complex(dp) :: grad_Hk(me%num_wann,me%num_wann,3)
   
      call fourier_R_to_k_deriv(kpt, me, me%ham_r, grad_Hk)

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
   function get_dipole(me,kpt,band_basis) result(Dk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      logical,intent(in),optional :: band_basis
      logical :: band_basis_
      logical :: large_size
      integer :: idir
      real(dp)    :: Ek(me%num_wann)
      complex(dp) :: Hk(me%num_wann,me%num_wann),Qk(me%num_wann,me%num_wann)
      complex(dp) :: Dk(me%num_wann,me%num_wann,3)
     
      band_basis_ = .false.
      if(present(band_basis)) band_basis_ = band_basis

      large_size = get_large_size(me%num_wann)

      call fourier_R_to_k_truevec(kpt, me, me%pos_r, Dk)

      if(me%force_herm) then
         do idir=1,3
            Dk(:,:,idir) = 0.5_dp * (Dk(:,:,idir) + conjg(transpose(Dk(:,:,idir))))
         end do
      end if

      if(band_basis_) then
         Hk = me%get_ham(kpt)
         call utility_diagonalize(Hk, me%num_wann, Ek, Qk)
         do idir=1,3
            Dk(:,:,idir) = util_rotate(me%num_wann,Qk,Dk(:,:,idir),large_size=large_size)
         end do
      end if

   end function get_dipole
!--------------------------------------------------------------------------------------
   function get_emp_velocity(me,kpt,orbs_excl) result(vk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      complex(dp)         :: vk(me%num_wann,me%num_wann,3)
      integer,intent(in),optional :: orbs_excl(:)
      logical :: large_size
      integer :: idir
      real(dp)    :: Ek(me%num_wann)
      complex(dp) :: Hk(me%num_wann,me%num_wann),Qk(me%num_wann,me%num_wann)

      large_size = get_large_size(me%num_wann)

      vk = me%get_gradk_ham(kpt)
      Hk = me%get_ham(kpt)
      call utility_diagonalize(Hk, me%num_wann, Ek, Qk)

      if(present(orbs_excl)) then
         call RemoveOrbitals(orbs_excl, Qk)
      end if

      do idir=1,3
         vk(:,:,idir) = util_rotate(me%num_wann,Qk,vk(:,:,idir),large_size=large_size)
      end do

   end function get_emp_velocity
!--------------------------------------------------------------------------------------
   function get_velocity(me,kpt,orbs_excl) result(Vk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      integer,intent(in),optional :: orbs_excl(:)
      integer :: idir
      complex(dp)         :: vk(me%num_wann,me%num_wann,3)

      logical :: large_size
      integer :: i,j
      real(dp) :: eig(me%num_wann)
      real(dp) :: del_eig(me%num_wann,3)
      complex(dp),allocatable :: HH(:,:)
      complex(dp),allocatable :: delHH(:,:,:)
      complex(dp),allocatable :: UU(:,:)
      complex(dp),allocatable :: D_h(:,:,:)
      complex(dp),allocatable :: AA(:,:,:)

      large_size = get_large_size(me%num_wann)

      allocate(HH(me%num_wann, me%num_wann))
      allocate(delHH(me%num_wann, me%num_wann, 3))
      allocate(UU(me%num_wann, me%num_wann))
      allocate(D_h(me%num_wann, me%num_wann, 3))
      allocate(AA(me%num_wann, me%num_wann, 3))

      call wham_get_eig_deleig(kpt, me, eig, del_eig, HH, delHH, UU, &
         use_degen_pert=me%use_degen_pert, degen_thr=me%degen_thresh)

      if(present(orbs_excl)) then
         call RemoveOrbitals(orbs_excl, UU)
      end if

      call wham_get_D_h(me%num_wann, delHH, UU, eig, D_h, &
         degen_thr=me%degen_thresh, anti_herm=me%force_antiherm)

      call fourier_R_to_k_truevec(kpt, me, me%pos_r, AA)
      do i = 1, 3
         AA(:, :, i) = util_rotate(me%num_wann, UU, AA(:, :, i), large_size=large_size)
         ! AA(:, :, i) = utility_rotate(AA(:, :, i), UU, me%num_wann)
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


      if(me%force_herm) then
         do i=1,3
            vk(:,:,i) = 0.5_dp * (vk(:,:,i) + conjg(transpose(vk(:,:,i))))
         end do
      end if

      deallocate(HH,delHH,UU,D_h,AA)

   end function get_velocity
!--------------------------------------------------------------------------------------
   function get_spin_velocity(me,kpt,orbs_excl) result(Vk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      integer,intent(in),optional :: orbs_excl(:)
      complex(dp)         :: vk(me%num_wann,me%num_wann,3,3)

      logical :: large_size
      integer :: i,j,idir,ispin,norb
      real(dp) :: eig(me%num_wann)
      real(dp) :: del_eig(me%num_wann,3)
      complex(dp),allocatable :: HH(:,:)
      complex(dp),allocatable :: delHH(:,:,:)
      complex(dp),allocatable :: UU(:,:)
      complex(dp),allocatable :: Dh_spin(:,:,:,:)
      complex(dp),allocatable :: AA(:,:,:),AA_spin(:,:,:,:)

      if(mod(me%num_wann,2) /= 0) then
         write(error_unit,fmt900) "num_wann is odd - not compatible with spin orbit coupling"
         stop
      end if

      norb = nint(me%num_wann / 2.0_dp)

      large_size = get_large_size(me%num_wann)

      allocate(HH(me%num_wann, me%num_wann))
      allocate(delHH(me%num_wann, me%num_wann, 3))
      allocate(UU(me%num_wann, me%num_wann))
      allocate(Dh_spin(me%num_wann, me%num_wann, 3, 3))
      allocate(AA(me%num_wann, me%num_wann, 3))
      allocate(AA_spin(me%num_wann, me%num_wann, 3, 3))

      call wham_get_eig_deleig(kpt, me, eig, del_eig, HH, delHH, UU, &
         use_degen_pert=me%use_degen_pert, degen_thr=me%degen_thresh)

      if(present(orbs_excl)) then
         call RemoveOrbitals(orbs_excl, UU)
      end if

      call wham_get_D_h_spin(me%num_wann, delHH, UU, eig, Dh_spin, &
         degen_thr=me%degen_thresh, anti_herm=me%force_antiherm)

      call fourier_R_to_k_truevec(kpt, me, me%pos_r, AA)

      do idir = 1, 3
         call GetSpinElements(me%num_wann, AA(:,:,idir), AA_spin(:,:,idir,:))
         do ispin = 1, 3
            AA_spin(:, :, idir, ispin) = util_rotate(me%num_wann, UU, AA_spin(:, :, idir, ispin), &
               large_size=large_size)
         end do
      end do
      AA_spin = AA_spin + iu*Dh_spin ! Eq.(25) WYSV06

      vk = zero
      do idir=1,3
         do i=1,norb
            vk(2*i-1,2*i-1,idir,3) = del_eig(2*i-1,idir)
            vk(2*i,2*i,idir,3) = -del_eig(2*i,idir)
         end do
      end do

      do i=1,me%num_wann
         do j=1,me%num_wann
            vk(j,i,:,:) = vk(j,i,:,:) - iu * (eig(i) - eig(j)) * AA_spin(j,i,:,:)
         end do
      end do


      if(me%force_herm) then
         do ispin=1,3
            do idir=1,3
               vk(:,:,idir,ispin) = 0.5_dp * (vk(:,:,idir,ispin) + conjg(transpose(vk(:,:,idir,ispin))))
            end do
         end do
      end if

      deallocate(HH,delHH,UU,Dh_spin,AA,AA_spin)

   end function get_spin_velocity
!--------------------------------------------------------------------------------------
   function get_berry_connection(me,kpt,orbs_excl) result(AA)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      integer,intent(in),optional :: orbs_excl(:)
      complex(dp)         :: AA(me%num_wann,me%num_wann,3)
      logical :: large_size
      integer :: i,j,idir
      real(dp) :: eig(me%num_wann)
      real(dp) :: del_eig(me%num_wann,3)
      complex(dp),allocatable :: HH(:,:)
      complex(dp),allocatable :: delHH(:,:,:)
      complex(dp),allocatable :: UU(:,:)
      complex(dp),allocatable :: D_h(:,:,:)

      large_size = get_large_size(me%num_wann)

      allocate(HH(me%num_wann, me%num_wann))
      allocate(delHH(me%num_wann, me%num_wann, 3))
      allocate(UU(me%num_wann, me%num_wann))
      allocate(D_h(me%num_wann, me%num_wann, 3))

      call wham_get_eig_deleig(kpt, me, eig, del_eig, HH, delHH, UU, &
         use_degen_pert=me%use_degen_pert, degen_thr=me%degen_thresh)

      if(present(orbs_excl)) then
         call RemoveOrbitals(orbs_excl, UU)
      end if

      call wham_get_D_h(me%num_wann, delHH, UU, eig, D_h, &
         degen_thr=me%degen_thresh, anti_herm=me%force_antiherm)

      call fourier_R_to_k_truevec(kpt, me, me%pos_r, AA)
      do idir = 1, 3
         AA(:, :, idir) = util_rotate(me%num_wann, UU, AA(:, :, idir), large_size=large_size)
      end do
      AA = AA + iu*D_h ! Eq.(25) WYSV06

      deallocate(HH,delHH,UU,D_h)

   end function get_berry_connection
!--------------------------------------------------------------------------------------
   function get_berrycurv(me, kpt, muchem, orbs_excl) result(Wk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      real(dp),intent(in),optional :: muchem
      integer,intent(in),optional :: orbs_excl(:)
      real(dp)            :: Wk(me%num_wann,3)
      logical :: kubo_
      logical :: large_size
      integer :: i,j
      real(dp),dimension(me%num_wann) :: epsk
      complex(dp),allocatable :: AA(:,:,:)

      kubo_ = present(muchem)

      large_size = get_large_size(me%num_wann)

      allocate(AA(me%num_wann, me%num_wann, 3))

      if(present(orbs_excl)) then
         AA = me%get_berry_connection(kpt, orbs_excl=orbs_excl)
      else
         AA = me%get_berry_connection(kpt)
      end if

      Wk = 0.0_dp

      if(kubo_) then
         epsk = me%get_eig(kpt)
         do i=1,me%num_wann
            if(epsk(i) > muchem) cycle
            do j=1,me%num_wann
               if(epsk(j) < muchem) cycle
               Wk(i,3) = Wk(i,3) + 2.0_dp * aimag(AA(i,j,1)*AA(j,i,2))
               Wk(i,1) = Wk(i,1) + 2.0_dp * aimag(AA(i,j,2)*AA(j,i,3))
               Wk(i,2) = Wk(i,2) + 2.0_dp * aimag(AA(i,j,3)*AA(j,i,1))
            end do
         end do
      else
         do i=1,me%num_wann
            do j=1,me%num_wann
               if(i == j) cycle
               Wk(i,3) = Wk(i,3) + 2.0_dp * aimag(AA(i,j,1)*AA(j,i,2))
               Wk(i,1) = Wk(i,1) + 2.0_dp * aimag(AA(i,j,2)*AA(j,i,3))
               Wk(i,2) = Wk(i,2) + 2.0_dp * aimag(AA(i,j,3)*AA(j,i,1))
            end do
         end do
      end if

      deallocate(AA)

   end function get_berrycurv
!--------------------------------------------------------------------------------------
   function get_spin_berrycurv(me, kpt, muchem, orbs_excl) result(Wk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      real(dp),intent(in),optional :: muchem
      integer,intent(in),optional :: orbs_excl(:)
      real(dp)            :: Wk(me%num_wann,3,3)
      logical :: kubo_
      logical :: large_size
      integer :: i,j,idir,ispin
      real(dp) :: eig(me%num_wann)
      real(dp) :: del_eig(me%num_wann,3)
      complex(dp),allocatable :: HH(:,:)
      complex(dp),allocatable :: delHH(:,:,:)
      complex(dp),allocatable :: UU(:,:)
      complex(dp),allocatable :: Dh_spin(:,:,:,:)
      complex(dp),allocatable :: AA(:,:,:),Dh(:,:,:)
      complex(dp),allocatable :: AA_spin(:,:,:,:)

      kubo_ = present(muchem)

      if(mod(me%num_wann,2) /= 0) then
         write(error_unit,fmt900) "num_wann is odd - not compatible with spin orbit coupling"
         stop
      end if

      large_size = get_large_size(me%num_wann)

      allocate(HH(me%num_wann, me%num_wann))
      allocate(delHH(me%num_wann, me%num_wann, 3))
      allocate(UU(me%num_wann, me%num_wann))
      allocate(AA(me%num_wann, me%num_wann, 3))
      allocate(Dh(me%num_wann, me%num_wann, 3))
      allocate(AA_spin(me%num_wann, me%num_wann, 3, 3))
      allocate(Dh_spin(me%num_wann, me%num_wann, 3, 3))

      call wham_get_eig_deleig(kpt, me, eig, del_eig, HH, delHH, UU, &
         use_degen_pert=me%use_degen_pert, degen_thr=me%degen_thresh)

      if(present(orbs_excl)) then
         call RemoveOrbitals(orbs_excl, UU)
      end if

      call wham_get_D_h(me%num_wann, delHH, UU, eig, Dh, &
         degen_thr=me%degen_thresh, anti_herm=me%force_antiherm)

      call wham_get_D_h_spin(me%num_wann, delHH, UU, eig, Dh_spin, &
         degen_thr=me%degen_thresh, anti_herm=me%force_antiherm)

      call fourier_R_to_k_truevec(kpt, me, me%pos_r, AA)

      do idir = 1, 3
         call GetSpinElements(me%num_wann, AA(:,:,idir), AA_spin(:,:,idir,:))
         do ispin = 1, 3
            AA_spin(:, :, idir,ispin) = util_rotate(me%num_wann, UU, AA_spin(:, :, idir, ispin), &
               large_size=large_size)
         end do
      end do
      AA_spin = AA_spin + iu*Dh_spin ! Eq.(25) WYSV06

      do idir = 1, 3
         AA(:, :, idir) = util_rotate(me%num_wann, UU, AA(:, :, idir), large_size=large_size)
      end do
      AA = AA + iu*Dh ! Eq.(25) WYSV06


      Wk = 0.0_dp

      if(kubo_) then
         do i=1,me%num_wann
            if(eig(i) > muchem) cycle
            do j=1,me%num_wann
               if(eig(j) < muchem) cycle
               Wk(i,3,ispin) = Wk(i,3,ispin) + 2.0_dp * aimag(AA_spin(i,j,1,ispin)*AA(j,i,2))
               Wk(i,1,ispin) = Wk(i,1,ispin) + 2.0_dp * aimag(AA_spin(i,j,2,ispin)*AA(j,i,3))
               Wk(i,2,ispin) = Wk(i,2,ispin) + 2.0_dp * aimag(AA_spin(i,j,3,ispin)*AA(j,i,1))
            end do
         end do     
      else
         do i=1,me%num_wann
            do j=1,me%num_wann
               if(i == j) cycle
               do ispin=1,3
                  Wk(i,3,ispin) = Wk(i,3,ispin) + 2.0_dp * aimag(AA_spin(i,j,1,ispin)*AA(j,i,2))
                  Wk(i,1,ispin) = Wk(i,1,ispin) + 2.0_dp * aimag(AA_spin(i,j,2,ispin)*AA(j,i,3))
                  Wk(i,2,ispin) = Wk(i,2,ispin) + 2.0_dp * aimag(AA_spin(i,j,3,ispin)*AA(j,i,1))
               end do
            end do
         end do
      end if

      deallocate(HH,delHH,UU,Dh,Dh_spin,AA,AA_spin)

   end function get_spin_berrycurv
!--------------------------------------------------------------------------------------
   subroutine get_berrycurv_dip(me, kpt, Bhk, Bdip, muchem, orbs_excl) 
      class(wann90_tb_t)   :: me
      real(dp),intent(in)  :: kpt(3)
      real(dp),intent(out) :: Bhk(me%num_wann,3)
      real(dp),intent(out) :: Bdip(me%num_wann,3)
      real(dp),intent(in),optional :: muchem
      integer,intent(in),optional :: orbs_excl(:)
      logical :: kubo_
      logical :: large_size
      integer :: i,j,idir
      real(dp) :: ediff
      real(dp) :: eig(me%num_wann)
      complex(dp),dimension(me%num_wann,me%num_wann) :: Hk,UU
      complex(dp),dimension(me%num_wann,me%num_wann,3) :: grad_Hk,Dk

      kubo_ = present(muchem)
      large_size = get_large_size(me%num_wann)

      Hk = me%get_ham(kpt)
      grad_Hk = me%get_gradk_ham(kpt)
      Dk = me%get_dipole(kpt)

      call utility_diagonalize(Hk, me%num_wann, eig, UU)

      if(present(orbs_excl)) then
         call RemoveOrbitals(orbs_excl, UU)
      end if

      do idir=1,3
         grad_Hk(:, :, idir) = util_rotate(me%num_wann, UU, grad_Hk(:, :, idir), large_size=large_size)
         Dk(:, :, idir) = util_rotate(me%num_wann, UU, Dk(:, :, idir), large_size=large_size)
         ! grad_Hk(:,:,idir) = utility_rotate(grad_Hk(:,:,idir), UU, me%num_wann)
         ! Dk(:,:,idir) = utility_rotate(Dk(:,:,idir), UU, me%num_wann)
      end do

      Bhk = 0.0_dp; Bdip = 0.0_dp

      if(kubo_) then
         do i=1,me%num_wann
            if(eig(i) > muchem) cycle
            do j=1,me%num_wann
               if(eig(j) < muchem) cycle
               ediff = eig(i) - eig(j)
               if(abs(ediff) < me%degen_thresh) cycle

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
      else
         do i=1,me%num_wann
            do j=1,me%num_wann
               if(i == j) cycle
               ediff = eig(i) - eig(j)
               if(abs(ediff) < me%degen_thresh) cycle

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
      end if

      Bdip = 2.0_dp * Bdip

   end subroutine get_berrycurv_dip
!--------------------------------------------------------------------------------------
   subroutine get_spin_berrycurv_dip(me, kpt, Bhk, Bdip, muchem, orbs_excl) 
      class(wann90_tb_t)   :: me
      real(dp),intent(in)  :: kpt(3)
      real(dp),intent(out) :: Bhk(me%num_wann,3,3)
      real(dp),intent(out) :: Bdip(me%num_wann,3,3)
      real(dp),intent(in),optional :: muchem
      integer,intent(in),optional :: orbs_excl(:)
      logical :: kubo_
      logical :: large_size
      integer :: i,j,idir,isig
      real(dp) :: ediff
      real(dp) :: eig(me%num_wann)
      complex(dp),dimension(me%num_wann,me%num_wann) :: Hk,UU
      complex(dp),dimension(me%num_wann,me%num_wann,3) :: grad_Hk,Dk,Mspin
      complex(dp),allocatable,dimension(:,:,:,:)       :: gHk_spin,Dk_spin

      if(mod(me%num_wann,2) /= 0) then
         write(error_unit,fmt900) "num_wann is odd - not compatible with spin orbit coupling"
         stop
      end if

      kubo_ = present(muchem)
      large_size = get_large_size(me%num_wann)

      Hk = me%get_ham(kpt)
      grad_Hk = me%get_gradk_ham(kpt)
      Dk = me%get_dipole(kpt)

      call utility_diagonalize(Hk, me%num_wann, eig, UU)

      if(present(orbs_excl)) then
         call RemoveOrbitals(orbs_excl, UU)
      end if

      allocate(gHk_spin(me%num_wann,me%num_wann,3,3),Dk_spin(me%num_wann,me%num_wann,3,3))

      do idir=1,3
         call GetSpinElements(me%num_wann, grad_Hk(:,:,idir), Mspin)
         do isig=1,3
            gHk_spin(:, :, idir, isig) = util_rotate(me%num_wann, UU, Mspin(:,:,isig),&
               large_size=large_size)
         end do
         call GetSpinElements(me%num_wann, Dk(:,:,idir), Mspin)
         do isig=1,3
            Dk_spin(:, :, idir, isig) = util_rotate(me%num_wann, UU, Mspin(:,:,isig),&
               large_size=large_size)
         end do
      end do

      do idir=1,3
         grad_Hk(:, :, idir) = util_rotate(me%num_wann, UU, grad_Hk(:, :, idir), large_size=large_size)
      end do

      Bhk = 0.0_dp; Bdip = 0.0_dp

      if(kubo_) then
         do i=1,me%num_wann
            if(eig(i) > muchem) cycle
            do j=1,me%num_wann
               if(eig(j) < muchem) cycle
               ediff = eig(i) - eig(j)
               if(abs(ediff) < me%degen_thresh) cycle

               do isig=1,3

                  Bhk(i,3,isig) = Bhk(i,3,isig) + aimag(grad_Hk(i,j,1) * gHk_spin(j,i,2,isig) &
                     - grad_Hk(i,j,2) * gHk_spin(j,i,1,isig)) / ediff**2
                  Bhk(i,1,isig) = Bhk(i,1,isig) + aimag(grad_Hk(i,j,2) * gHk_spin(j,i,3,isig) &
                     - grad_Hk(i,j,3) * gHk_spin(j,i,2,isig)) / ediff**2
                  Bhk(i,2,isig) = Bhk(i,2,isig) + aimag(grad_Hk(i,j,3) * gHk_spin(j,i,1,isig) - &
                     grad_Hk(i,j,1) * gHk_spin(j,i,3,isig)) / ediff**2

                  Bdip(i,3,isig) = Bdip(i,3,isig) - dble(grad_Hk(i,j,1) * Dk_spin(j,i,2,isig) &
                     - grad_Hk(i,j,2) * Dk_spin(j,i,1,isig)) / ediff
                  Bdip(i,1,isig) = Bdip(i,1,isig) - dble(grad_Hk(i,j,2) * Dk_spin(j,i,3,isig) &
                     - grad_Hk(i,j,3) * Dk_spin(j,i,2,isig)) / ediff
                  Bdip(i,2,isig) = Bdip(i,2,isig) - dble(grad_Hk(i,j,3) * Dk_spin(j,i,1,isig) &
                     - grad_Hk(i,j,1) * Dk_spin(j,i,3,isig)) / ediff

               end do
            end do
         end do          
      else
         do i=1,me%num_wann
            do j=1,me%num_wann
               if(i == j) cycle
               ediff = eig(i) - eig(j)
               if(abs(ediff) < me%degen_thresh) cycle

               do isig=1,3

                  Bhk(i,3,isig) = Bhk(i,3,isig) + aimag(grad_Hk(i,j,1) * gHk_spin(j,i,2,isig) &
                     - grad_Hk(i,j,2) * gHk_spin(j,i,1,isig)) / ediff**2
                  Bhk(i,1,isig) = Bhk(i,1,isig) + aimag(grad_Hk(i,j,2) * gHk_spin(j,i,3,isig) &
                     - grad_Hk(i,j,3) * gHk_spin(j,i,2,isig)) / ediff**2
                  Bhk(i,2,isig) = Bhk(i,2,isig) + aimag(grad_Hk(i,j,3) * gHk_spin(j,i,1,isig) - &
                     grad_Hk(i,j,1) * gHk_spin(j,i,3,isig)) / ediff**2

                  Bdip(i,3,isig) = Bdip(i,3,isig) - dble(grad_Hk(i,j,1) * Dk_spin(j,i,2,isig) &
                     - grad_Hk(i,j,2) * Dk_spin(j,i,1,isig)) / ediff
                  Bdip(i,1,isig) = Bdip(i,1,isig) - dble(grad_Hk(i,j,2) * Dk_spin(j,i,3,isig) &
                     - grad_Hk(i,j,3) * Dk_spin(j,i,2,isig)) / ediff
                  Bdip(i,2,isig) = Bdip(i,2,isig) - dble(grad_Hk(i,j,3) * Dk_spin(j,i,1,isig) &
                     - grad_Hk(i,j,1) * Dk_spin(j,i,3,isig)) / ediff

               end do
            end do
         end do 
      end if

      Bdip = 2.0_dp * Bdip

      deallocate(gHk_spin,Dk_spin)

   end subroutine get_spin_berrycurv_dip
!--------------------------------------------------------------------------------------
   function get_oam(me, kpt, orbs_excl) result(Lk)
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3)
      integer,intent(in),optional :: orbs_excl(:)
      real(dp)            :: Lk(me%num_wann,3)

      logical :: large_size
      integer :: i,j,idir
      real(dp) :: eig(me%num_wann)
      real(dp) :: del_eig(me%num_wann,3)
      complex(dp),allocatable :: HH(:,:)
      complex(dp),allocatable :: delHH(:,:,:)
      complex(dp),allocatable :: UU(:,:)
      complex(dp),allocatable :: D_h(:,:,:)
      complex(dp),allocatable :: AA(:,:,:)

      large_size = get_large_size(me%num_wann)

      allocate(HH(me%num_wann, me%num_wann))
      allocate(delHH(me%num_wann, me%num_wann, 3))
      allocate(UU(me%num_wann, me%num_wann))
      allocate(D_h(me%num_wann, me%num_wann, 3))
      allocate(AA(me%num_wann, me%num_wann, 3))

      call wham_get_eig_deleig(kpt, me, eig, del_eig, HH, delHH, UU, &
         use_degen_pert=me%use_degen_pert, degen_thr=me%degen_thresh)

      if(present(orbs_excl)) then
         call RemoveOrbitals(orbs_excl, UU)
      end if

      call wham_get_D_h(me%num_wann, delHH, UU, eig, D_h, &
         degen_thr=me%degen_thresh, anti_herm=me%force_antiherm)

      call fourier_R_to_k_truevec(kpt, me, me%pos_r, AA)
      do idir = 1, 3
         ! AA(:, :, idir) = utility_rotate(AA(:, :, idir), UU, me%num_wann)
         AA(:, :, idir) = util_rotate(me%num_wann, UU, AA(:, :, idir), large_size=large_size)
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
   function get_oam_dip(me, kpt, orbs_excl) result(Lk)
      class(wann90_tb_t)   :: me
      real(dp),intent(in)  :: kpt(3)
      integer,intent(in),optional :: orbs_excl(:)
      real(dp)             :: Lk(me%num_wann,3)

      logical :: large_size
      integer :: i,j,idir
      real(dp) :: ediff
      real(dp) :: eig(me%num_wann)
      real(dp),dimension(me%num_wann,3) :: Lk_disp,Lk_dip
      complex(dp),dimension(me%num_wann,me%num_wann) :: Hk,UU
      complex(dp),dimension(me%num_wann,me%num_wann,3) :: grad_Hk,Dk

      large_size = get_large_size(me%num_wann)

      Hk = me%get_ham(kpt)
      grad_Hk = me%get_gradk_ham(kpt)
      Dk = me%get_dipole(kpt)

      call utility_diagonalize(Hk, me%num_wann, eig, UU)

      if(present(orbs_excl)) then
         call RemoveOrbitals(orbs_excl, UU)
      end if

      do idir=1,3
         grad_Hk(:, :, idir) = util_rotate(me%num_wann, UU, grad_Hk(:, :, idir), large_size=large_size)
         Dk(:, :, idir) = util_rotate(me%num_wann, UU, Dk(:, :, idir), large_size=large_size)
         ! grad_Hk(:,:,idir) = utility_rotate(grad_Hk(:,:,idir), UU, me%num_wann)
         ! Dk(:,:,idir) = utility_rotate(Dk(:,:,idir), UU, me%num_wann)
      end do

      Lk_disp = 0.0_dp; Lk_dip = 0.0_dp
      do i=1,me%num_wann
         do j=1,me%num_wann
            if(i == j) cycle
            ediff = eig(i) - eig(j)
            if(abs(ediff) < me%degen_thresh) cycle

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
   function get_metric(me,kpt,muchem,orbs_excl) result(gmet)
      !! computes the quantum metric \(g^\alpha_{\mu\nu}(k)\) from the Berry connections
      class(wann90_tb_t)  :: me
      real(dp),intent(in) :: kpt(3) !! k-point (reduced coordinates)
      real(dp),intent(in),optional :: muchem !! chemical potential
      integer,intent(in),optional :: orbs_excl(:)
      real(dp) :: gmet(me%num_wann,3,3)
      logical :: kubo_
      logical :: large_size
      integer :: i,j,idir,jdir
      real(dp) :: eig(me%num_wann)
      real(dp) :: del_eig(me%num_wann,3)
      complex(dp),allocatable :: HH(:,:)
      complex(dp),allocatable :: delHH(:,:,:)
      complex(dp),allocatable :: UU(:,:)
      complex(dp),allocatable :: D_h(:,:,:)
      complex(dp),allocatable :: AA(:,:,:)

      kubo_ = present(muchem)
      large_size = get_large_size(me%num_wann)

      allocate(HH(me%num_wann, me%num_wann))
      allocate(delHH(me%num_wann, me%num_wann, 3))
      allocate(UU(me%num_wann, me%num_wann))
      allocate(D_h(me%num_wann, me%num_wann, 3))
      allocate(AA(me%num_wann, me%num_wann, 3))

      call wham_get_eig_deleig(kpt, me, eig, del_eig, HH, delHH, UU, &
         use_degen_pert=me%use_degen_pert, degen_thr=me%degen_thresh)

      if(present(orbs_excl)) then
         call RemoveOrbitals(orbs_excl, UU)
      end if

      call wham_get_D_h(me%num_wann, delHH, UU, eig, D_h, &
         degen_thr=me%degen_thresh, anti_herm=me%force_antiherm)

      call fourier_R_to_k_truevec(kpt, me, me%pos_r, AA)
      do idir = 1, 3
         AA(:, :, idir) = util_rotate(me%num_wann, UU, AA(:, :, idir), large_size=large_size)
         ! AA(:, :, i) = utility_rotate(AA(:, :, i), UU, me%num_wann)
      end do
      AA = AA + iu*D_h ! Eq.(25) WYSV06

      gmet = 0.0_dp

      if(kubo_) then
         do jdir=1,3
            do idir=1,3
               do i=1,me%num_wann
                  if(eig(i) > muchem) cycle
                  do j=1,me%num_wann
                     if(eig(j) < muchem) cycle
                     gmet(i,idir,jdir) = gmet(i,idir,jdir) &
                        + 0.5_dp * dble(conjg(AA(i,j,idir)) * AA(i,j,jdir) + conjg(AA(i,j,jdir)) * AA(i,j,idir))
                  end do
               end do
            end do
         end do
      else
         do jdir=1,3
            do idir=1,3
               do i=1,me%num_wann
                  do j=1,me%num_wann
                     if(i == j) cycle
                     gmet(i,idir,jdir) = gmet(i,idir,jdir) &
                        + 0.5_dp * dble(conjg(AA(i,j,idir)) * AA(i,j,jdir) + conjg(AA(i,j,jdir)) * AA(i,j,idir))
                  end do
               end do
            end do
         end do
      end if

      deallocate(HH,delHH,UU,D_h,AA)

   end function get_metric
!--------------------------------------------------------------------------------------
   subroutine ReadFromW90(me,file_ham,file_xyz)
      class(wann90_tb_t)  :: me
      character(len=*),intent(in) :: file_ham
      character(len=*),intent(in),optional :: file_xyz
      integer :: ir

      call ReadTB_from_w90(file_ham,me%real_lattice,me%num_wann,me%nrpts,me%ndegen,&
         me%irvec,me%ham_r,me%pos_r)

      allocate(me%coords(me%num_wann,3)); me%coords = 0.0_dp
      if(present(file_xyz)) then
         if(len_trim(file_xyz) > 0) then
            call ReadXYZ(file_xyz,me%coords)
            me%coords = me%coords / BohrAngstrom
            me%coords_present = .true.
         end if
      end if

      call utility_recip_lattice(me%real_lattice, me%recip_lattice)

      call utility_recip_reduced(me%recip_lattice, me%recip_reduced)

      call get_crvec(me,me%crvec)

   end subroutine ReadFromW90
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine SaveToHDF5(me,fname,atomic_units)
      use scitools_hdf5_utils
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

         allocate(d_pos_r(me%num_wann,me%num_wann,3,me%nrpts))

         d_pos_r = dble(me%pos_r)
         call hdf_write_dataset(file_id,'pos_r_real',d_pos_r)
         d_pos_r = aimag(me%pos_r)
         call hdf_write_dataset(file_id,'pos_r_imag',d_pos_r)

         deallocate(d_pos_r)
      else
         call hdf_write_attribute(file_id,'','pos_stored', 0)
      end if

      if(me%coords_present) then
         call hdf_write_attribute(file_id,'','coords_stored', 1)
         call hdf_write_dataset(file_id,'coords',me%coords)
      else
         call hdf_write_attribute(file_id,'','coords_stored', 0)        
      end if

   end subroutine SaveToHDF5
#endif
!--------------------------------------------------------------------------------------
#if WITHHDF5
   subroutine ReadFromHDF5(me,fname)
      !! Reads the Wannier Hamiltonian from HDF5 binary format
      use scitools_hdf5_utils
      class(wann90_tb_t)  :: me
      character(len=*),intent(in) :: fname
      integer(HID_t) :: file_id
      integer :: atomic_units,pos_stored,coords_stored
      integer :: ir
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
         allocate(d_pos_r(me%num_wann,me%num_wann,3,me%nrpts))
         allocate(me%pos_r(me%num_wann,me%num_wann,3,me%nrpts))

         call hdf_read_dataset(file_id,'pos_r_real',d_pos_r)
         me%pos_r = d_pos_r
         call hdf_read_dataset(file_id,'pos_r_imag',d_pos_r)
         me%pos_r = me%pos_r + iu * d_pos_r

         deallocate(d_pos_r)

         if(atomic_units == 0) me%pos_r = me%pos_r / BohrAngstrom
      end if

      allocate(me%coords(me%num_wann,3)); me%coords = 0.0_dp
      call hdf_read_attribute(file_id,'','coords_stored', coords_stored)  
      if(coords_stored == 1) then
         me%coords_present = .true.
         call hdf_read_dataset(file_id,'coords',me%coords)
         if(atomic_units == 0) me%coords = me%coords / BohrAngstrom
      end if

      call hdf_close_file(file_id)

      call get_crvec(me,me%crvec)

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
      allocate(pos_r(num_wann,num_wann,3,nrpts))
      do irpt = 1,nrpts
         read(file_unit,*) irvec(irpt,1),irvec(irpt,2),irvec(irpt,3)
         do i=1,num_wann
            do j=1,num_wann
               ! read(file_unit,*) ndx1, ndx2, pos_real, pos_imag
               ! pos_r(j,i,irpt,:) = (pos_real + iu*pos_imag) / BohrAngstrom
               read(file_unit,'(2I5,3x,6(E15.8,1x))') ndx1, ndx2, pos
               pos_r(j,i,:,irpt) = pos / BohrAngstrom
            end do
         end do
      end do

      close(file_unit)

   end subroutine ReadTB_from_w90
!--------------------------------------------------------------------------------------
   subroutine SaveToW90(me,fname)
      !! Save the Wannier Hamiltonian to file, matching the format of Wannier90
      character(len=*),parameter :: fmt_time='(a,"/",a,"/",a," at ",a,":",a,":",a)'
      class(wann90_tb_t)  :: me
      character(len=*),intent(in) :: fname
      character(len=8) :: date
      character(len=10) :: time
      integer :: i,j,irpt
      integer :: file_unit
      
      open(newunit=file_unit, file=trim(fname), status='replace')

      call date_and_time(date=date,time=time)
      write(file_unit, fmt_time) date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
      !
      ! lattice vectors
      !
      write(file_unit, *) BohrAngstrom * me%real_lattice(1, :) !a_1
      write(file_unit, *) BohrAngstrom * me%real_lattice(2, :) !a_2
      write(file_unit, *) BohrAngstrom * me%real_lattice(3, :) !a_3
      
      write(file_unit, *) me%num_wann
      write(file_unit, *) me%nrpts
      write(file_unit, '(15I5)') (me%ndegen(i), i=1, me%nrpts)
      !
      ! <0n|H|Rm>
      !
      do irpt = 1, me%nrpts
         write(file_unit, '(/,3I5)') me%irvec(irpt,:)
         do i = 1, me%num_wann
            do j = 1, me%num_wann
               write(file_unit, '(2I5,3x,2(E15.8,1x))') j, i, me%ham_r(j, i, irpt) * HreV
            end do
         end do
      end do
      !
      ! <0n|r|Rm>
      !
      do irpt = 1, me%nrpts
         write(file_unit, '(/,3I5)') me%irvec(irpt,:)
         do i = 1, me%num_wann
            do j = 1, me%num_wann
               write(file_unit, '(2I5,3x,6(E15.8,1x))') j, i, me%pos_r(j, i, 1:3, irpt) * BohrAngstrom
            end do
         end do
      end do

      close(file_unit)

   end subroutine SaveToW90
!--------------------------------------------------------------------------------------
   subroutine wham_get_D_h(num_wann, delHH, UU, eig, D_h, degen_thr, anti_herm)
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
      complex(kind=dp), dimension(:, :, :), intent(inout) :: D_h
      real(dp), intent(in), optional :: degen_thr
      logical, intent(in), optional  :: anti_herm

      real(dp) :: degen_thr_
      logical :: anti_herm_
      logical :: large_size
      complex(kind=dp), allocatable :: delHH_bar_i(:, :)
      integer                       :: n, m, i

      degen_thr_ = 1.0e-7_dp
      if(present(degen_thr)) degen_thr_ = degen_thr

      anti_herm_ = .true.
      if(present(anti_herm)) anti_herm_ = anti_herm

      large_size = get_large_size(num_wann)

      allocate (delHH_bar_i(num_wann, num_wann))
      D_h = zero
      do i = 1, 3
         delHH_bar_i(:, :) = util_rotate(num_wann, UU, delHH(:, :, i), large_size=large_size)
         ! delHH_bar_i(:, :) = utility_rotate(delHH(:, :, i), UU, num_wann)
         do m = 1, num_wann
            do n = 1, num_wann
               if (n == m .or. abs(eig(m) - eig(n)) < degen_thr_) cycle
               D_h(n, m, i) = delHH_bar_i(n, m)/(eig(m) - eig(n))
            end do
         end do
      enddo

      if(anti_herm_) then
         do i=1,3
            D_h(:, :, i) = 0.5_dp*(D_h(:, :, i) - conjg(transpose(D_h(:, :, i))))
         end do
      end if

      deallocate(delHH_bar_i)

   end subroutine wham_get_D_h
!--------------------------------------------------------------------------------------
   subroutine wham_get_D_h_spin(num_wann, delHH, UU, eig, D_h, degen_thr, anti_herm)
      integer,intent(in) :: num_wann 
      complex(kind=dp), dimension(:, :, :), intent(in)  :: delHH
      complex(kind=dp), dimension(:, :), intent(in)    :: UU
      real(kind=dp), dimension(:), intent(in)    :: eig
      complex(kind=dp), dimension(:, :, :, :), intent(inout) :: D_h
      real(dp), intent(in), optional :: degen_thr
      logical, intent(in), optional  :: anti_herm

      real(dp) :: degen_thr_
      logical :: anti_herm_
      logical :: large_size
      complex(kind=dp), allocatable :: delHH_bar_i(:, :), delHH_spin(:,:,:)
      integer                       :: n, m, i, mu


      degen_thr_ = 1.0e-7_dp
      if(present(degen_thr)) degen_thr_ = degen_thr

      anti_herm_ = .true.
      if(present(anti_herm)) anti_herm_ = anti_herm

      large_size = get_large_size(num_wann)

      allocate(delHH_bar_i(num_wann, num_wann), delHH_spin(num_wann, num_wann, 3))
      D_h = zero
      do i = 1, 3
         call GetSpinElements(num_wann, delHH(:, :, i), delHH_spin)
         do mu = 1, 3
            delHH_bar_i(:, :) = util_rotate(num_wann, UU, delHH_spin(:, :, mu), large_size=large_size)
            ! delHH_bar_i(:, :) = utility_rotate(delHH(:, :, i), UU, num_wann)
            do m = 1, num_wann
               do n = 1, num_wann
                  if (n == m .or. abs(eig(m) - eig(n)) < degen_thr_) cycle
                  D_h(n, m, i, mu) = delHH_bar_i(n, m)/(eig(m) - eig(n))
               end do
            end do
         end do    
      end do    

   end subroutine wham_get_D_h_spin
!--------------------------------------------------------------------------------------
   subroutine wham_get_eig_deleig(kpt, w90, eig, del_eig, HH, delHH, UU, use_degen_pert, degen_thr)
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
      complex(kind=dp), dimension(:, :), intent(inout)   :: HH
      !! the Hamiltonian matrix at kpt
      complex(kind=dp), dimension(:, :, :), intent(inout) :: delHH
      !! the delHH matrix (derivative of H) at kpt
      complex(kind=dp), dimension(:, :), intent(inout)   :: UU
      !! the rotation matrix that gives the eigenvectors of HH
      logical, intent(in), optional :: use_degen_pert
      real(dp),intent(in), optional :: degen_thr
      logical :: use_degen_pert_
      real(dp) :: degen_thr_

      use_degen_pert_ = .false.
      if(present(use_degen_pert)) use_degen_pert_ = use_degen_pert

      degen_thr_ = 1.0e-5_dp
      if(present(degen_thr)) degen_thr_ = degen_thr

      call fourier_R_to_k(kpt, w90, w90%ham_r, HH)
      call utility_diagonalize(HH, w90%num_wann, eig, UU)
      call fourier_R_to_k_deriv(kpt, w90, w90%ham_r, delHH)
      ! call fourier_R_to_k(kpt, w90, w90%ham_r, delHH(:, :, 1), 1)
      ! call fourier_R_to_k(kpt, w90, w90%ham_r, delHH(:, :, 2), 2)
      ! call fourier_R_to_k(kpt, w90, w90%ham_r, delHH(:, :, 3), 3)
      call wham_get_deleig_a(w90%num_wann, del_eig(:, 1), w90, eig, delHH(:, :, 1), UU, use_degen_pert_, degen_thr_)
      call wham_get_deleig_a(w90%num_wann, del_eig(:, 2), w90, eig, delHH(:, :, 2), UU, use_degen_pert_, degen_thr_)
      call wham_get_deleig_a(w90%num_wann, del_eig(:, 3), w90, eig, delHH(:, :, 3), UU, use_degen_pert_, degen_thr_)

   end subroutine wham_get_eig_deleig
!--------------------------------------------------------------------------------------
   subroutine wham_get_deleig_a(num_wann, deleig_a, w90, eig, delHH_a, UU, use_degen_pert, degen_thr)
      !! Band derivatives \( dE/dk_a \)
      integer,intent(in) :: num_wann !! number of Wannier functions 
      real(kind=dp), intent(out) :: deleig_a(num_wann)
      type(wann90_tb_t),intent(in) :: w90
      real(kind=dp), intent(in)  :: eig(num_wann)
      complex(kind=dp), dimension(:, :), intent(in)  :: delHH_a
      complex(kind=dp), dimension(:, :), intent(in)  :: UU
      logical,intent(in),optional  :: use_degen_pert
      real(dp),intent(in),optional :: degen_thr

      ! Misc/Dummy
      !
      logical :: use_degen_pert_
      real(dp) :: degen_thr_
      logical :: large_size
      integer                       :: i, degen_min, degen_max, dim
      real(kind=dp)                 :: diff
      complex(kind=dp), allocatable :: delHH_bar_a(:, :), U_deg(:, :)

      use_degen_pert_=.false.
      if(present(use_degen_pert)) use_degen_pert_ = use_degen_pert

      degen_thr_ = 1.0e-5_dp
      if(present(degen_thr)) degen_thr_ = degen_thr

      large_size = get_large_size(num_wann)

      allocate (delHH_bar_a(num_wann, num_wann))
      allocate (U_deg(num_wann, num_wann))

      if (use_degen_pert_) then

         ! delHH_bar_a = utility_rotate(delHH_a, UU, num_wann)
         delHH_bar_a = util_rotate(num_wann, UU, delHH_a, large_size)

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
      allocate(crvec(w90%nrpts,3))
      do ir = 1, w90%nrpts
        ! Note that 'real_lattice' stores the lattice vectors as *rows*
        crvec(ir, :) = matmul(transpose(w90%real_lattice), w90%irvec(ir, :))
      end do

   end subroutine
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k(kpt, w90, OO_R, OO)
      !! Performs the Fourier transformation R -> k
      !! For \(\alpha=0\): 
      !! \(O_{ij}(R) \rightarrow O_{ij}(k) = \sum_R e^{i k \cdot R} O_{ij}(R)\)
      !! For \(\alpha=1,2,3\):
      !! \(i \sum_R R_\alpha e^{i k \cdot R} O_{ij}(R) \)
      real(kind=dp)                                     :: kpt(3) !! k-point (reduced coordinates)
      type(wann90_tb_t),intent(in)                      :: w90 !! Wannier90 object
      complex(kind=dp), dimension(:, :, :), intent(in)  :: OO_R !! operator in real space O(R)
      complex(kind=dp), dimension(:, :), intent(inout)  :: OO !! operator in k-space O(k)

      integer          :: ir,nn,nr,m
      complex(kind=dp) :: phase_fac(blocksize)
      integer :: numblock,imin,imax

      ! compute the number of chuncks
      numblock  = (w90%nrpts+blocksize-1)/blocksize
      nn = (w90%num_wann)**2

      OO(:, :) = zero
      do m = 1, numblock
         imin = (m-1)*blocksize+1
         imax = min(m * blocksize, w90%nrpts)
         nr = imax-imin+1
         call GetPhase(nr,kpt,w90%irvec(imin:imax,:),w90%ndegen(imin:imax),phase_fac)

         call ZGEMV("N",nn,nr,one,OO_R(1,1,imin),nn,phase_fac(1),1,one,OO(1,1),1)

         ! do ir = imin, imax
         !    OO(:, :) = OO(:, :) + phase_fac(ir - imin + 1)*OO_R(:, :, ir)
         ! end do      
      end do

   end subroutine fourier_R_to_k
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k_deriv(kpt, w90, OO_R, OO)
      !! Performs the Fourier transformation R -> k
      !! For \(\alpha=0\): 
      !! \(O_{ij}(R) \rightarrow O_{ij}(k) = \sum_R e^{i k \cdot R} O_{ij}(R)\)
      !! For \(\alpha=1,2,3\):
      !! \(i \sum_R R_\alpha e^{i k \cdot R} O_{ij}(R) \)
      real(kind=dp)                                     :: kpt(3) !! k-point (reduced coordinates)
      type(wann90_tb_t),intent(in)                      :: w90 !! Wannier90 object
      complex(kind=dp), dimension(:, :, :), intent(in)  :: OO_R !! operator in real space O(R)
      complex(kind=dp), dimension(:, :, :), intent(inout)  :: OO !! operator in k-space O(k)

      integer          :: m,nn,nr,ir,i,j,idir
      integer :: numblock,imin,imax
      complex(kind=dp) :: phase_fac(blocksize),r_phase(blocksize)

      ! compute the number of chuncks
      numblock  = (w90%nrpts+blocksize-1)/blocksize
      nn = (w90%num_wann)**2

      OO(:, :, :) = zero
      do m = 1, numblock
         imin = (m-1)*blocksize+1
         imax = min(m * blocksize, w90%nrpts)
         nr = imax-imin+1
         call GetPhase(nr,kpt,w90%irvec(imin:imax,:),w90%ndegen(imin:imax),phase_fac)

         do idir=1,3
            r_phase(1:nr) = iu * w90%crvec(imin:imax,idir) * phase_fac(1:nr)
            call ZGEMV("N",nn,nr,one,OO_R(1,1,imin),nn,r_phase(1),1,one,OO(1,1,idir),1)
         end do

         ! do concurrent(ir=imin:imax, idir=1:3, j=1:w90%num_wann, i=1:w90%num_wann)
         !    OO(i, j, idir) = OO(i, j, idir) + iu * w90%crvec(ir,idir) &
         !       * phase_fac(ir-imin+1) * OO_R(i,j,ir)
         ! end do

      end do

   end subroutine fourier_R_to_k_deriv
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k_slab(ijmax, kpt, w90, OO_R, OO, alpha)
      !=========================================================!
      !
      !! 2D Fourier Transformation for the slab calculation
      !!
      !! For alpha=0:
      !! O_ij(R) --> O_ij(k) = sum_R e^{+ik.R}*O_ij(R)
      !!
      !! For alpha=1,2,3:
      !! sum_R [cmplx_i*R_alpha*e^{+ik.R}*O_ij(R)]
      !! where R_alpha is a Cartesian component of R
      !! ***REMOVE EVENTUALLY*** (replace with pw90common_fourier_R_to_k_new)
      !!
      !! TO DO:
      !! (1) transformation for random direction of irvec (U)
      !! (2) slab(logic) --> input parameters into the type w90
      !
      !=========================================================!

      ! Arguments
      !
      integer,intent(in)                                            :: ijmax
      real(kind=dp)                                                 :: kpt(2)
      type(wann90_tb_t),intent(in)                                  :: w90
      complex(kind=dp), dimension(:, :, :), intent(in)              :: OO_R
      complex(kind=dp), dimension(:, :, :), intent(inout)           :: OO
      integer                                                       :: alpha

      integer          :: ir, i, j, ideg, i3
      real(kind=dp)    :: rdotk
      complex(kind=dp) :: phase_fac
      real(dp),allocatable :: crvec(:,:)

      if(alpha > 0) call get_crvec(w90, crvec)

      OO(:, :, :) = zero
      do ir = 1, w90%nrpts
         rdotk = DPI*dot_product(kpt(:), w90%irvec(ir, :))
         phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(w90%ndegen(ir), dp)
         i3 = w90%irvec(ir, 3)
         if (abs(i3) < ijmax) then
             if (alpha == 0) then
                OO(:, :,i3+ijmax+1) = OO(:, :,i3+ijmax+1) + &
                    phase_fac*OO_R(:, :, ir)
             elseif (alpha == 1 .or. alpha == 2 .or. alpha == 3) then
                OO(:, :,i3+ijmax+1) = OO(:, :,i3+ijmax+1) + &
                   iu*crvec(ir, alpha)*phase_fac*OO_R(:, :, ir)
             else
                stop 'wrong value of alpha in fourier_R_to_k_2D'
             end if
         end if

      end do

      if(allocated(crvec)) deallocate(crvec)

   end subroutine fourier_R_to_k_slab
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

      OO(:, :) = zero
      do ir = 1, w90%nrpts
         rdotk = DPI*dot_product(kpt(:), w90%irvec(ir, :))
         phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(w90%ndegen(ir), dp)
         OO(:, :) = OO(:, :) - &
            w90%crvec(ir, a)*w90%crvec(ir, b)*phase_fac*OO_R(:, :, ir)

      end do

   end subroutine fourier_D2_R_to_k
!--------------------------------------------------------------------------------------
   subroutine GetPhase(nr,kpt,irvecs,dgens,exp_iphase) 
      integer,intent(in)  :: nr
      real(dp),intent(in) :: kpt(3)
      integer,intent(in)  :: irvecs(:,:)
      integer,intent(in)  :: dgens(:)
      complex(dp),intent(inout) :: exp_iphase(:)
      real(dp),dimension(blocksize) :: s,c,rdotk

      rdotk(1:nr) = DPI*(kpt(1) * irvecs(1:nr,1) + kpt(2) * irvecs(1:nr,2) &
         + kpt(3) * irvecs(1:nr,3))
      c(1:nr) = cos(rdotk(1:nr))
      s(1:nr) = sin(rdotk(1:nr))
      exp_iphase(1:nr) = cmplx(c(1:nr), s(1:nr), kind=dp) / dgens(1:nr)

   end subroutine GetPhase
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k_truevec(kpt, w90, OO_R, OO_true)
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
      complex(kind=dp), dimension(:, :, :), intent(inout)   :: OO_true

      integer          :: nr,ir,m,nn
      complex(kind=dp) :: phase_fac(blocksize)
      integer :: numblock,imin,imax

      ! compute the number of chuncks
      numblock  = (w90%nrpts+blocksize-1)/blocksize
      nn = 3 * (w90%num_wann)**2

      OO_true = zero

      do m = 1, numblock
         imin = (m-1)*blocksize+1
         imax = min(m * blocksize, w90%nrpts)
         nr = imax-imin+1
         call GetPhase(nr,kpt,w90%irvec(imin:imax,:),w90%ndegen(imin:imax),phase_fac)
         ! do ir=imin,imax
         !    OO_true(:,:,:) = OO_true(:,:,:) + phase_fac(ir-imin+1) * OO_R(:,:,:,ir)
         ! end do
         call ZGEMV("N",nn,nr,one,OO_R(1,1,1,imin),nn,phase_fac(1),1,one,OO_true(1,1,1),1)
      end do

   end subroutine fourier_R_to_k_truevec
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
            OO_true(:, :, 1) = OO_true(:, :, 1) + phase_fac*OO_R(:, :, 1, ir)
            OO_true(:, :, 2) = OO_true(:, :, 2) + phase_fac*OO_R(:, :, 2, ir)
            OO_true(:, :, 3) = OO_true(:, :, 3) + phase_fac*OO_R(:, :, 3, ir)
         end if
         if (present(OO_pseudo)) then
            OO_pseudo(:, :, 1) = OO_pseudo(:, :, 1) &
               + iu*crvec(ir, 2)*phase_fac*OO_R(:, :, 3, ir) &
               - iu*crvec(ir, 3)*phase_fac*OO_R(:, :, 2, ir)
            OO_pseudo(:, :, 2) = OO_pseudo(:, :, 2) &
               + iu*crvec(ir, 3)*phase_fac*OO_R(:, :, 1, ir) &
               - iu*crvec(ir, 1)*phase_fac*OO_R(:, :, 3, ir)
            OO_pseudo(:, :, 3) = OO_pseudo(:, :, 3) &
               + iu*crvec(ir, 1)*phase_fac*OO_R(:, :, 2, ir) &
               - iu*crvec(ir, 2)*phase_fac*OO_R(:, :, 1, ir)
         end if

      end do

      if(allocated(crvec)) deallocate(crvec)

   end subroutine fourier_R_to_k_vec
!--------------------------------------------------------------------------------------
   subroutine GetSpinElements(num_wann,A,Aspin)
      integer,intent(in)        :: num_wann
      complex(dp),intent(in)    :: A(:,:)
      complex(dp),intent(inout) :: Aspin(:,:,:)
      integer :: norb
      integer :: i,j

      call assert_shape(A, [num_wann, num_wann], "Mham_w90: GetSpinElements", "A")
      call assert_shape(Aspin, [num_wann, num_wann, 3], "Mham_w90: GetSpinElements", "Aspin")

      norb = nint(num_wann / 2.0_dp)

      Aspin = zero

      do j=1,norb
         do i=1,norb
            Aspin(2*i-1, 2*j, 1) = A(2*i-1, 2*j) ! <up|\sigma_x|dn>
            Aspin(2*i, 2*j-1, 1) = A(2*i, 2*j-1) ! <dn|\sigma_x|up>
            Aspin(2*i-1, 2*j, 2) = -iu*A(2*i-1, 2*j) ! <up|\sigma_y|dn>
            Aspin(2*i, 2*j-1, 2) = iu*A(2*i, 2*j-1) ! <dn|\sigma_y|up>
            Aspin(2*i-1, 2*j-1, 3) = A(2*i-1, 2*j-1) ! <up|\sigma_z|up>
            Aspin(2*i, 2*j, 3) = -A(2*i, 2*j) ! <dn|\sigma_z|dn>
         end do
      end do

   end subroutine GetSpinElements
!--------------------------------------------------------------------------------------
   subroutine RemoveOrbitals(orbs,UU)
      integer,intent(in) :: orbs(:)
      complex(dp),intent(inout) :: UU(:,:)
      integer :: nexc,i

      nexc = size(orbs, dim=1)
      if(nexc <= 0) return

      if(nexc == 1 .and. orbs(1) == 0) return

      do i=1,nexc
         UU(orbs(i),:) = zero
      end do

   end subroutine RemoveOrbitals
!--------------------------------------------------------------------------------------

!======================================================================================
end module wan_hamiltonian

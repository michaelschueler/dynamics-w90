module Mwannier_calc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero
   use Mlinalg,only: EigHE
   use Mlatt_kpts,only: Read_Kpoints
   use Mham_w90,only: wann90_tb_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: wannier_calc_t
!--------------------------------------------------------------------------------------
   type wannier_calc_t
      logical :: soc_mode=.false.
      integer :: nbnd,nwan,Nk
      real(dp),allocatable,dimension(:,:) :: kpts,epsk
      complex(dp),allocatable,dimension(:,:,:) :: vectk
      type(wann90_tb_t) :: ham
   contains
      procedure,public  :: Init
      procedure,public  :: GetOrbitalWeight
      procedure,public  :: GetSpin
      procedure,public  :: GetBerryCurvature
      procedure,public  :: GetOAM
   end type wannier_calc_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,file_inp)
      class(wannier_calc_t) :: me
      character(len=*),intent(in) :: file_inp
      character(len=255) :: file_ham
      logical :: w90_with_soc=.false.
      integer :: unit_inp,ppos
      namelist/HAMILTONIAN/file_ham,w90_with_soc
      complex(dp),allocatable :: Hk(:,:)

      open(newunit=unit_inp,file=trim(file_inp),status='OLD',action='READ')
      read(unit_inp,nml=HAMILTONIAN) 
      close(unit_inp)

      ppos = scan(trim(file_ham),".", BACK= .true.)
      if(trim(file_ham(ppos+1:)) == "h5") then
         call me%Ham%ReadFromHDF5(file_ham)
      else
         call me%Ham%ReadFromW90(file_ham)
      end if

      me%nwan = me%ham%num_wann
      me%nbnd = me%nwan

      call Read_Kpoints(file_inp,me%kpts,print_info=.true.)

      me%Nk = size(kpts,1)
      allocate(me%epsk(me%nbnd,me%Nk),me%vectk(me%nbnd,me%nbnd,me%Nk))

      allocate(Hk(me%nwan,me%nwan))
      do ik=1,me%Nk
         Hk = me%Ham%get_ham(kpts(ik,:))
         call EigHE(Hk,me%epsk(:,ik),me%vectk(:,:,ik))
      end do
      deallocate(Hk)

      me%soc_mode = w90_with_soc

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine GetOrbitalWeight(me,orb_weight)
      class(wannier_calc_t) :: me
      real(dp),allocatable,intent(out) :: orb_weight(:,:,:)
      integer :: iorb

      allocate(orb_weight(me%nwan,me%nbnd,me%Nk)); orb_weight = 0.0_dp
      orb_weight = abs(me%vectk)**2

   end subroutine GetOrbitalWeight
!--------------------------------------------------------------------------------------
   subroutine GetSpin(me,spin)
      class(wannier_calc_t) :: me
      real(dp),allocatable,intent(out) :: spin(:,:,:)
      integer :: ibnd

      allocate(spin(3,me%nbnd,me%Nk)); spin = 0.0_dp

      ! ... TODO ....
      if(me%soc_mode) then
         write(output_unit,fmt700) "Spin projection with w90_with_soc not implemented yet!"
      end if

   end subroutine GetSpin
!--------------------------------------------------------------------------------------
   subroutine GetBerryCurvature(me,berry,gauge)
      class(wannier_calc_t) :: me
      real(dp),allocatable,intent(out) :: berry(:,:,:)
      integer,intent(in),optional :: gauge
      integer :: gauge_
      integer :: ibnd
      complex(dp),allocatable :: berry_disp(:,:),berry_dip(:,:)

      gauge_ = 0
      if(present(gauge)) gauge_ = gauge

      allocate(berry(me%nbnd,3,me%Nk)); berry = 0.0_dp

      select case(gauge_)
      case(0) 
         do ik=1,me%Nk
            berry(:,:,ik) = me%Ham%get_berrycurv(me%kpts(ik,:))
         end do
      case(1)
         allocate(berry_disp(me%nbnd,3),berry_dip(me%nbnd,3))
         do ik=1,me%Nk
            call me%Ham%get_berrycurv_dip(me%kpts(ik,:),berry_disp,berry_dip)
            berry(:,:,ik) = berry_disp + berry_dip
         deallocate(berry_disp,berry_dip)
      case default
         write(output_unit,fmt900) "GetBerryCurvature: unrecognized gauge!"
         return
      end select

   end subroutine GetBerryCurvature
!--------------------------------------------------------------------------------------

   subroutine ham_soc_ReadFromHDF5(me,fname)
      use Mhdf5_utils
      class(ham_soc_t) :: me
      character(len=*),intent(in) :: fname
      integer(HID_t) :: file_id
      real(dp),allocatable :: rdata(:,:,:)
      logical :: file_ok

      inquire(file=trim(fname),exist=file_ok)
      if(.not.file_ok) then
         write(error_unit,fmt900) 'Input file does not exist: '//trim(fname)
         stop
      end if

      call hdf_open_file(file_id, trim(fname), STATUS='OLD', ACTION='READ')

      call hdf_read_attribute(file_id,'','ngroups', me%ngroups)
      allocate(me%ndim(me%ngroups),me%Lorb(me%ngroups))
      call hdf_read_dataset(file_id,'ndim',me%ndim)
      call hdf_read_dataset(file_id,'lorb',me%Lorb)

      me%norb = sum(me%ndim)

      allocate(me%Lmat(me%norb,me%norb,3))
      allocate(rdata(me%norb,me%norb,3))

      call hdf_read_dataset(file_id,'hsoc_real',rdata)
      me%Lmat = rdata
      call hdf_read_dataset(file_id,'hsoc_imag',rdata)
      me%Lmat = me%Lmat + iu * rdata

      deallocate(rdata)

      call hdf_close_file(file_id) 

   end subroutine ham_soc_ReadFromHDF5
!--------------------------------------------------------------------------------------
   subroutine Read_SOC_lambda(fname,lam)
      character(len=*),intent(in) :: fname
      real(dp),allocatable,intent(out) :: lam(:)
      integer :: ngroups
      logical :: file_ok
      integer :: iunit

      inquire(file=trim(fname),exist=file_ok)
      if(.not.file_ok) then
         write(error_unit,fmt900) 'Input file does not exist: '//trim(fname)
         stop
      end if

      open(newunit=iunit,file=trim(fname),status='OLD',action='READ')
      read(iunit,*) ngroups
      allocate(lam(ngroups))
      read(iunit,*) lam
      close(iunit)

   end subroutine Read_SOC_lambda
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
   subroutine GetHam_SOC(kpt,wann,soc,lam,Hk) 
      real(dp),intent(in)           :: kpt(3)
      type(wann90_tb_t),intent(in)  :: wann
      type(ham_soc_t),intent(in)    :: soc
      real(dp),intent(in)           :: lam(:)
      complex(dp),intent(inout)     :: Hk(:,:)
      integer :: norb,nbnd,ig,imin,imax
      complex(dp),allocatable :: Hk_nosoc(:,:),H0_soc(:,:)

      norb = soc%norb
      nbnd = 2 * norb

      call assert(soc%norb == wann%num_wann, "GetHam_SOC: soc%norb == wann%num_wann")
      call assert_shape(lam, [soc%ngroups], "GetHam_SOC", "lam")
      call assert_shape(Hk, [nbnd,nbnd], "GetHam_SOC", "Hk")

      allocate(Hk_nosoc(norb,norb))
      Hk_nosoc = wann%get_ham(kpt)

      Hk = zero
      Hk(1:norb,1:norb) = Hk_nosoc
      Hk(norb+1:nbnd,norb+1:nbnd) = Hk_nosoc

      allocate(H0_soc(nbnd,nbnd))
      call GetHam_SOC_onsite(soc,lam,H0_soc)

      Hk = Hk + H0_soc

      deallocate(Hk_nosoc,H0_soc)

   end subroutine GetHam_SOC
!--------------------------------------------------------------------------------------
   subroutine GetHam_NOSOC(kpt,wann,Hk)
      real(dp),intent(in)           :: kpt(3)
      type(wann90_tb_t),intent(in)  :: wann
      complex(dp),intent(inout)     :: Hk(:,:)
      complex(dp),allocatable :: Hk_nosoc(:,:)
      integer :: norb,nbnd

      norb = wann%num_wann
      nbnd = 2 * norb

      call assert_shape(Hk, [nbnd,nbnd], "GetHam_NOSOC", "Hk")

      allocate(Hk_nosoc(norb,norb))
      Hk_nosoc = wann%get_ham(kpt)

      Hk = zero
      Hk(1:norb,1:norb) = Hk_nosoc
      Hk(norb+1:nbnd,norb+1:nbnd) = Hk_nosoc

      deallocate(Hk_nosoc)

   end subroutine GetHam_NOSOC
!--------------------------------------------------------------------------------------
   subroutine GetHam_SOC_onsite(soc,lam,H0)
      type(ham_soc_t),intent(in)    :: soc
      real(dp),intent(in)           :: lam(:)
      complex(dp),intent(inout)     :: H0(:,:)
      integer :: norb,nbnd,ig,imin,imax
      integer,allocatable :: orb_min(:)   

      norb = soc%norb
      nbnd = 2 * norb

      call assert_shape(lam, [soc%ngroups], "GetHam_SOC_onsite", "lam")
      call assert_shape(H0, [nbnd,nbnd], "GetHam_SOC_onsite", "H0")

      allocate(orb_min(soc%ngroups)); orb_min = 1
      do ig=1,soc%ngroups-1
         orb_min(ig+1) = orb_min(ig) + soc%ndim(ig)
      end do

      H0 = zero
      do ig=1,soc%ngroups
         imin = orb_min(ig)
         imax = orb_min(ig) + soc%ndim(ig) - 1

         ! sigma_x
         H0(imin:imax,imin+norb:imax+norb) = H0(imin:imax,imin+norb:imax+norb) &
            + lam(ig) * soc%Lmat(imin:imax,imin:imax,1)
         H0(imin+norb:imax+norb,imin:imax) = H0(imin+norb:imax+norb,imin:imax) &
            + lam(ig) * soc%Lmat(imin:imax,imin:imax,1)

         ! sigma_y
         H0(imin:imax,imin+norb:imax+norb) = H0(imin:imax,imin+norb:imax+norb) &
            - iu * lam(ig) * soc%Lmat(imin:imax,imin:imax,2)
         H0(imin+norb:imax+norb,imin:imax) = H0(imin+norb:imax+norb,imin:imax) &
            + iu * lam(ig) * soc%Lmat(imin:imax,imin:imax,2)

         ! sigma_z
         H0(imin:imax,imin:imax) = H0(imin:imax,imin:imax) &
            + lam(ig) * soc%Lmat(imin:imax,imin:imax,3)
         H0(imin+norb:imax+norb,imin+norb:imax+norb) = H0(imin+norb:imax+norb,imin+norb:imax+norb) &
            - lam(ig) * soc%Lmat(imin:imax,imin:imax,3)
      end do

   end subroutine GetHam_SOC_onsite
!--------------------------------------------------------------------------------------
   subroutine GetObs_SOC(soc,vect,Lvec,Svec)
      type(ham_soc_t),intent(in) :: soc
      complex(dp),intent(in)     :: vect(:)
      real(dp),intent(out)       :: Lvec(3),Svec(3)
      integer :: norb,ig,imin,imax
      integer,allocatable :: orb_min(:)   

      call assert_shape(vect, [2*soc%norb], "GetObs_SOC", "vect")

      Lvec = 0.0_dp
      Svec = 0.0_dp

      norb = soc%norb
      allocate(orb_min(soc%ngroups)); orb_min = 1
      do ig=1,soc%ngroups-1
         orb_min(ig+1) = orb_min(ig) + soc%ndim(ig)
      end do

      do ig=1,soc%ngroups
         imin = orb_min(ig)
         imax = orb_min(ig) + soc%ndim(ig) - 1

         Lvec(1) = Lvec(1) + exp_val(soc%Lmat(:,:,1), vect, imin, imax, imin, imax) &
            + exp_val(soc%Lmat(:,:,1), vect, imin, imax, imin+norb, imax+norb)
         Lvec(2) = Lvec(2) + exp_val(soc%Lmat(:,:,2), vect, imin, imax, imin, imax) &
            + exp_val(soc%Lmat(:,:,2), vect, imin, imax, imin+norb, imax+norb)
         Lvec(3) = Lvec(3) + exp_val(soc%Lmat(:,:,3), vect, imin, imax, imin, imax) &
            + exp_val(soc%Lmat(:,:,3), vect, imin, imax, imin+norb, imax+norb)

         Svec(1) = Svec(1) + 2.0_dp * dble(dot_product(vect(imin:imax), vect(imin+norb:imax+norb)))
         Svec(2) = Svec(2) + 2.0_dp * aimag(dot_product(vect(imin:imax), vect(imin+norb:imax+norb)))
         Svec(3) = Svec(3) + sum(abs(vect(imin:imax))**2 - abs(vect(imin+norb:imax+norb))**2)

      end do

      deallocate(orb_min)

   end subroutine GetObs_SOC
!--------------------------------------------------------------------------------------
   real(dp) function exp_val(A,v,istart_A,iend_A,istart_v,iend_v)
      complex(dp),intent(in) :: A(:,:)
      complex(dp),intent(in) :: v(:)
      integer,intent(in) :: istart_A,iend_A,istart_v,iend_v

      exp_val = dble(dot_product(v(istart_v:iend_v), &
         matmul(A(istart_A:iend_A,istart_A:iend_A), v(istart_v:iend_v))))

   end function exp_val
!--------------------------------------------------------------------------------------

!======================================================================================    
end module Mwannier_calc
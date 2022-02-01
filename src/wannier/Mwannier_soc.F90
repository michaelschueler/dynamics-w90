module Mwannier_soc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero
   use Munits,only: DPi
   use Mham_w90,only: wann90_tb_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: ham_soc_t
   public :: GetHam_SOC, GetHam_NOSOC, GetHam_SOC_onsite, GetObs_SOC, Read_SOC_lambda
!--------------------------------------------------------------------------------------
   type ham_soc_t
      integer                                  :: ngroups,norb
      integer,allocatable,dimension(:)         :: ndim,Lorb
      complex(dp),allocatable,dimension(:,:,:) :: Lmat
   contains
      procedure,public  :: ReadFromHDF5 => ham_soc_ReadFromHDF5
   end type ham_soc_t
!--------------------------------------------------------------------------------------
contains
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
end module Mwannier_soc
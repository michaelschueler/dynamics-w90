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
      integer :: ik
      complex(dp),allocatable :: Hk(:,:)

      open(newunit=unit_inp,file=trim(file_inp),status='OLD',action='READ')
      read(unit_inp,nml=HAMILTONIAN) 
      close(unit_inp)

      ppos = scan(trim(file_ham),".", BACK= .true.)
      if(trim(file_ham(ppos+1:)) == "h5") then
#if WITHHDF5
      call me%Ham%ReadFromHDF5(file_ham)
#else
      write(error_unit,fmt900) "No HDF5 support. Can't read "//trim(file_ham)
      stop
#endif
      else
         call me%Ham%ReadFromW90(file_ham)
      end if

      me%soc_mode = w90_with_soc

      if(me%soc_mode) then
         if(mod(me%ham%num_wann,2) /= 0) then
            write(error_unit,fmt900) "Number of Wannier orbitals is odd!"
            stop
         end if
         me%nwan = int(me%ham%num_wann/2.0_dp)
      else
         me%nwan = me%ham%num_wann
      end if
      me%nbnd = me%ham%num_wann

      call Read_Kpoints(file_inp,me%kpts,print_info=.true.)

      me%Nk = size(me%kpts,1)
      allocate(me%epsk(me%nbnd,me%Nk),me%vectk(me%nbnd,me%nbnd,me%Nk))

      allocate(Hk(me%nbnd,me%nbnd))
      do ik=1,me%Nk
         Hk = me%Ham%get_ham(me%kpts(ik,:))
         call EigHE(Hk,me%epsk(:,ik),me%vectk(:,:,ik))
      end do
      deallocate(Hk)

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine GetOrbitalWeight(me,orb_weight)
      class(wannier_calc_t) :: me
      real(dp),allocatable,intent(out) :: orb_weight(:,:,:)
      integer :: iorb

      allocate(orb_weight(me%nwan,me%nbnd,me%Nk)); orb_weight = 0.0_dp
      if(me%soc_mode) then
         do iorb=1,me%nwan
            orb_weight(iorb,:,:) = abs(me%vectk(2*iorb-1,:,:))**2 & 
               + abs(me%vectk(2*iorb,:,:))**2 
         end do
      else
         orb_weight = abs(me%vectk)**2
      end if

   end subroutine GetOrbitalWeight
!--------------------------------------------------------------------------------------
   subroutine GetSpin(me,spin)
      class(wannier_calc_t) :: me
      real(dp),allocatable,intent(out) :: spin(:,:,:)
      integer :: ibnd,iorb,ik
      complex(dp),allocatable :: vectk_up(:),vectk_dn(:)

      allocate(spin(3,me%nbnd,me%Nk)); spin = 0.0_dp
      if(.not.me%soc_mode) return

      allocate(vectk_up(me%nwan),vectk_dn(me%nwan))
      do ik=1,me%Nk
         do ibnd=1,me%nbnd
            do iorb=1,me%nwan
               vectk_up(iorb) = me%vectk(2*iorb-1,ibnd,ik)
               vectk_dn(iorb) = me%vectk(2*iorb,ibnd,ik)
            end do
            spin(1,ibnd,ik) = 2.0_dp*dble(dot_product(vectk_up,vectk_dn))
            spin(2,ibnd,ik) = 2.0_dp*aimag(dot_product(vectk_up,vectk_dn))
            spin(3,ibnd,ik) = sum(abs(vectk_up)**2) - sum(abs(vectk_dn)**2) 
         end do
      end do
     
   end subroutine GetSpin
!--------------------------------------------------------------------------------------
   subroutine GetBerryCurvature(me,berry,gauge)
      class(wannier_calc_t) :: me
      real(dp),allocatable,intent(out) :: berry(:,:,:)
      integer,intent(in),optional :: gauge
      integer :: gauge_
      integer :: ik
      real(dp),allocatable :: berry_disp(:,:),berry_dip(:,:)

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
         end do
         deallocate(berry_disp,berry_dip)
      case default
         write(error_unit,fmt900) "GetBerryCurvature: unrecognized gauge!"
         return
      end select

   end subroutine GetBerryCurvature
!--------------------------------------------------------------------------------------
   subroutine GetOAM(me,oam,gauge)
      class(wannier_calc_t) :: me
      real(dp),allocatable,intent(out) :: oam(:,:,:)
      integer,intent(in),optional :: gauge
      integer :: gauge_
      integer :: ik

      gauge_ = 0
      if(present(gauge)) gauge_ = gauge

      allocate(oam(me%nbnd,3,me%Nk)); oam = 0.0_dp

      select case(gauge_)
      case(0) 
         do ik=1,me%Nk
            oam(:,:,ik) = me%Ham%get_oam(me%kpts(ik,:))
         end do
      case(1)
         do ik=1,me%Nk
            oam(:,:,ik) = me%Ham%get_oam_dip(me%kpts(ik,:))
         end do
      case default
         write(error_unit,fmt900) "GetOAM: unrecognized gauge!"
         return
      end select

   end subroutine GetOAM
!--------------------------------------------------------------------------------------


!======================================================================================    
end module Mwannier_calc
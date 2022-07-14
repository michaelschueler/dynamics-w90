module Mwannier_calc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero
   use Mutils,only: str
   use Mlinalg,only: EigHE
   use Mlatt_kpts,only: Read_Kpoints
   use Mham_w90,only: wann90_tb_t
   use Mwann_compress,only: PruneHoppings
   use Mwann_slab,only: Wannier_BulkToSlab
   use Mio_hamiltonian,only: ReadHamiltonian
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
      procedure,public  :: GetSpinBerryCurvature
      procedure,public  :: GetOAM
      procedure,public  :: GetMetric
   end type wannier_calc_t
!--------------------------------------------------------------------------------------
   integer,parameter :: velocity_gauge=0,dipole_gauge=1
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,file_inp,slab_mode)
      class(wannier_calc_t) :: me
      character(len=*),intent(in) :: file_inp
      logical,intent(in)          :: slab_mode
      character(len=255) :: file_ham
      logical  :: w90_with_soc=.false.,expert_params=.false.
      real(dp) :: energy_thresh=0.0_dp
      namelist/HAMILTONIAN/file_ham,w90_with_soc,energy_thresh,expert_params
      integer :: nlayer=0,max_zhop=10
      namelist/SLAB/nlayer,max_zhop
      !......................................
      integer :: unit_inp
      integer :: ik
      real(dp) :: comp_rate
      complex(dp),allocatable :: Hk(:,:)
      type(wann90_tb_t) :: ham_tmp
      !......................................

      open(newunit=unit_inp,file=trim(file_inp),status='OLD',action='READ')
      read(unit_inp,nml=HAMILTONIAN); rewind(unit_inp)
      if(slab_mode)  read(unit_inp,nml=SLAB); rewind(unit_inp)
      close(unit_inp)

      call ReadHamiltonian(file_ham,ham_tmp)
      if(energy_thresh > 0.0_dp) then
         call PruneHoppings(energy_thresh,ham_tmp,me%ham,comp_rate)
         write(output_unit,fmt_info) "compression rate: "//str(nint(100 * comp_rate)) // "%"
      else
         call me%ham%Set(ham_tmp)
      end if
      call ham_tmp%Clean()

      if(expert_params) call me%ham%ReadParams(file_inp) 

      me%soc_mode = w90_with_soc

      if(slab_mode .and. nlayer > 0) then
         write(output_unit,fmt_info) "building slab with "//str(nlayer)//" layers"
         call ham_tmp%Set(me%ham)
         call Wannier_BulkToSlab(ham_tmp,nlayer,me%ham,ijmax=max_zhop)
      end if
      call ham_tmp%Clean()

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
            spin(ibnd,1,ik) = 2.0_dp*dble(dot_product(vectk_up,vectk_dn))
            spin(ibnd,2,ik) = 2.0_dp*aimag(dot_product(vectk_up,vectk_dn))
            spin(ibnd,3,ik) = sum(abs(vectk_up)**2) - sum(abs(vectk_dn)**2) 
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
      case(velocity_gauge) 
         do ik=1,me%Nk
            berry(:,:,ik) = me%Ham%get_berrycurv(me%kpts(ik,:))
         end do
      case(dipole_gauge)
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
   subroutine GetSpinBerryCurvature(me,spin_berry,gauge)
      class(wannier_calc_t) :: me
      real(dp),allocatable,intent(out) :: spin_berry(:,:,:,:)
      integer,intent(in),optional :: gauge
      integer :: gauge_
      integer :: ik
      real(dp),allocatable :: berry_disp(:,:,:),berry_dip(:,:,:)

      gauge_ = 0
      if(present(gauge)) gauge_ = gauge

      allocate(spin_berry(me%nbnd,3,3,me%Nk)); spin_berry = 0.0_dp

      select case(gauge_)
      case(velocity_gauge) 
         do ik=1,me%Nk
            spin_berry(:,:,:,ik) = me%Ham%get_spin_berrycurv(me%kpts(ik,:))
         end do
      case(dipole_gauge)
         allocate(berry_disp(me%nbnd,3,3),berry_dip(me%nbnd,3,3))
         do ik=1,me%Nk
            call me%Ham%get_spin_berrycurv_dip(me%kpts(ik,:),berry_disp,berry_dip)
            spin_berry(:,:,:,ik) = berry_disp + berry_dip
         end do
         deallocate(berry_disp,berry_dip)
      case default
         write(error_unit,fmt900) "GetBerryCurvature: unrecognized gauge!"
         return
      end select

   end subroutine GetSpinBerryCurvature
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
      case(velocity_gauge) 
         do ik=1,me%Nk
            oam(:,:,ik) = me%Ham%get_oam(me%kpts(ik,:))
         end do
      case(dipole_gauge)
         do ik=1,me%Nk
            oam(:,:,ik) = me%Ham%get_oam_dip(me%kpts(ik,:))
         end do
      case default
         write(error_unit,fmt900) "GetOAM: unrecognized gauge!"
         return
      end select

   end subroutine GetOAM
!--------------------------------------------------------------------------------------
   subroutine GetMetric(me,metric)
      class(wannier_calc_t) :: me
      real(dp),allocatable,intent(out) :: metric(:,:,:,:)
      integer :: ik

      allocate(metric(me%nbnd,3,3,me%Nk))
      do ik=1,me%Nk
         metric(:,:,:,ik) = me%Ham%get_metric(me%kpts(ik,:))
      end do

   end subroutine GetMetric
!--------------------------------------------------------------------------------------

!======================================================================================    
end module Mwannier_calc

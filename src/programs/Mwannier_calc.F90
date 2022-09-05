module Mwannier_calc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero
   use Mutils,only: str
   use Mlinalg,only: EigHE
   use Mham_w90,only: wann90_tb_t
   use Mwann_compress,only: PruneHoppings
   use Mwann_slab,only: Wannier_BulkToSlab
   use Mio_params,only: WannierCalcParams_t, HamiltonianParams_t
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
   integer,parameter :: field_mode_positions=0,field_mode_dipole=1,field_mode_berry=2
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,par_ham,par_calc,kpts)
      class(wannier_calc_t) :: me
      type(HamiltonianParams_t),intent(in) :: par_ham
      type(WannierCalcParams_t),intent(in) :: par_calc
      real(dp),intent(in)                  :: kpts(:,:)
      !......................................
      integer :: ik
      real(dp) :: comp_rate
      complex(dp),allocatable :: Hk(:,:)
      type(wann90_tb_t) :: ham_tmp
      !......................................

      call ReadHamiltonian(par_ham%file_ham,ham_tmp,file_xyz=par_ham%file_xyz)
      if(par_ham%energy_thresh > 0.0_dp) then
         call PruneHoppings(par_ham%energy_thresh,ham_tmp,me%ham,comp_rate)
         write(output_unit,fmt_info) "compression rate: "//str(nint(100 * comp_rate)) // "%"
      else
         call me%ham%Set(ham_tmp)
      end if
      call ham_tmp%Clean()

      me%soc_mode = par_ham%w90_with_soc

      if(par_ham%slab_mode .and. par_ham%slab_nlayer > 0) then
         write(output_unit,fmt_info) "building slab with "//str(par_ham%slab_nlayer)//" layers"
         call ham_tmp%Set(me%ham)
         call Wannier_BulkToSlab(ham_tmp,par_ham%slab_nlayer,me%ham,ijmax=par_ham%slab_max_zhop)
      end if
      call ham_tmp%Clean()

      call me%ham%SetExpertParams(use_degen_pert=par_ham%use_degen_pert,&
         degen_thresh=par_ham%degen_thresh,&
         force_herm=par_ham%force_herm,&
         force_antiherm=par_ham%force_antiherm)

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

      allocate(me%kpts(size(kpts,1), size(kpts,2)))
      me%kpts = kpts

      me%Nk = size(me%kpts,1)
      allocate(me%epsk(me%nbnd,me%Nk),me%vectk(me%nbnd,me%nbnd,me%Nk))

      allocate(Hk(me%nbnd,me%nbnd))
      if(par_ham%apply_field) then
         write(output_unit,fmt_info) "E-field: ["//str(par_ham%Efield(1))//" , "//str(par_ham%Efield(2)) &
            // " , "//str(par_ham%Efield(3))//"]"
         select case(par_ham%field_mode)
         case(field_mode_positions)
            write(output_unit,fmt_info) "field couples to positions"
         case(field_mode_dipole)
            write(output_unit,fmt_info) "field couples to dipoles"
         case(field_mode_berry)
            write(output_unit,fmt_info) "field couples to Berry connection"
         end select 
         do ik=1,me%Nk
            Hk = me%Ham%get_ham_field(me%kpts(ik,:),par_ham%Efield,par_ham%field_mode)
            call EigHE(Hk,me%epsk(:,ik),me%vectk(:,:,ik))
         end do
      else
         do ik=1,me%Nk
            Hk = me%Ham%get_ham(me%kpts(ik,:))
            call EigHE(Hk,me%epsk(:,ik),me%vectk(:,:,ik))
         end do
      end if
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

module Marpes_calc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero
   use Mutils,only: str
   use Mlatt_utils,only: utility_Cart2Red_2D
   use Mham_w90,only: wann90_tb_t
   use Mwannier_orbitals,only: wannier_orbs_t
   use Mscattwf,only: scattwf_t
   use Mpes_intensity,only: PES_Intensity
   use Mio_params,only: HamiltonianParams_t, PESParams_t
   use Mio_hamiltonian,only: ReadHamiltonian
   use Mio_orbitals,only: ReadWannierOrbitals
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: arpes_calc_t
!--------------------------------------------------------------------------------------
   type arpes_calc_t
      integer     :: nbnd
      integer     :: gauge
      integer     :: Nepe
      real(dp)    :: wphot,MuChem,Eshift,lambda_esc,eta_smear
      complex(dp) :: polvec(3)
      integer     :: Nk
      real(dp),allocatable,dimension(:)   :: Epe
      real(dp),allocatable,dimension(:,:) :: kpts,spect
      type(wann90_tb_t)    :: ham
      type(wannier_orbs_t) :: orbs
      type(scattwf_t)      :: chi
   contains
      procedure,public  :: Init
      procedure,public  :: CalcPES
      procedure,public  :: WriteSpectrum
      procedure,private :: WriteSpectrum_txt
#ifdef WITHHDF5
      procedure,private :: WriteSpectrum_hdf5
#endif      
   end type arpes_calc_t
!--------------------------------------------------------------------------------------
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Init(me,par_ham,par_pes,kpts)
      use Mutils,only: linspace
      class(arpes_calc_t) :: me
      type(HamiltonianParams_t),intent(in) :: par_ham
      type(PESParams_t),intent(in)         :: par_pes
      real(dp),intent(in)                  :: kpts(:,:)
      integer :: ik
      real(dp) :: kvec(3),kred(3)

      call ReadHamiltonian(par_ham%file_ham,me%ham,file_xyz=par_ham%file_xyz)
      me%nbnd = me%ham%num_wann
      me%MuChem = par_ham%MuChem

      call ReadWannierOrbitals(par_pes%file_orbs,me%orbs)
      me%gauge = par_pes%gauge
      me%Nepe = par_pes%Nepe
      me%wphot = par_pes%wphot
      me%Eshift = par_pes%Eshift
      me%lambda_esc = par_pes%lambda_esc
      me%eta_smear = par_pes%eta_smear    
      me%polvec = par_pes%polvec

      allocate(me%Epe(me%Nepe))
      if(me%Nepe == 1) then
         allocate(me%Epe(1))
         me%Epe(1) = par_pes%Epe_min
      else
         me%Epe = linspace(par_pes%Epe_min, par_pes%Epe_max, me%Nepe)
      end if

      call me%chi%Init(par_pes%scatt_type,par_pes%Zeff)

      me%Nk = size(kpts, dim=1)
      allocate(me%kpts(me%Nk,2))
      if(par_pes%kpts_reduced) then
         do ik=1,me%Nk
            me%kpts(ik,1:2) = me%ham%recip_lattice(1,1:2) * kpts(ik,1) + &
               me%ham%recip_lattice(2,1:2) * kpts(ik,2) 
         end do
      else
         me%kpts(1:me%Nk,1:2) = kpts(1:me%Nk,1:2)
      end if

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine CalcPES(me)
      use Mlinalg,only: EigHE
      class(arpes_calc_t) :: me
      integer :: ik,iepe
      real(dp) :: kpar(2),kpt(3)
      real(dp),allocatable :: epsk(:)
      complex(dp),allocatable :: Hk(:,:),vectk(:,:)

      allocate(epsk(me%nbnd),Hk(me%nbnd,me%nbnd),vectk(me%nbnd,me%nbnd))


      allocate(me%spect(me%Nepe,me%Nk))

      kpt = 0.0_dp
      do ik=1,me%Nk
         kpar(1:2) = me%kpts(ik,1:2)
         kpt(1:2) = utility_Cart2Red_2D(me%ham%recip_reduced,kpar)

         print*, kpt

         Hk = me%ham%get_ham(kpt)
         call EigHE(Hk,epsk,vectk)
         epsk = epsk + me%Eshift

         !$OMP PARALLEL DO
         do iepe=1,me%Nepe
            me%spect(iepe,ik) = PES_Intensity(me%orbs,me%ham,me%chi,kpar,me%wphot,&
               me%polvec,me%Epe(iepe),epsk,vectk,me%MuChem,me%lambda_esc,me%eta_smear,me%gauge)
         end do
      end do

      deallocate(epsk,Hk,vectk)

   end subroutine CalcPES
!--------------------------------------------------------------------------------------
   subroutine WriteSpectrum(me,prefix)
      class(arpes_calc_t) :: me
      character(len=*),intent(in) :: prefix

#ifdef WITHHDF5
      call me%WriteSpectrum_hdf5(prefix)
#else
      call me%WriteSpectrum_txt(prefix)
#endif

   end subroutine WriteSpectrum 
!--------------------------------------------------------------------------------------
   subroutine WriteSpectrum_txt(me,prefix)
      use Mutils,only: savetxt
      class(arpes_calc_t) :: me
      character(len=*),intent(in) :: prefix

      call savetxt(trim(prefix)//'_pes.txt', me%spect, transp=.true.)

   end subroutine WriteSpectrum_txt
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine WriteSpectrum_hdf5(me,prefix)
      use Mhdf5_utils
      class(arpes_calc_t) :: me
      character(len=*),intent(in) :: prefix
      integer(HID_t) :: file_id

      call hdf_open_file(file_id, trim(prefix)//'_pes.h5', STATUS='NEW')      
      call hdf_write_dataset(file_id,'epe',me%Epe)
      call hdf_write_dataset(file_id,'spect',me%spect)     
      call hdf_close_file(file_id)

   end subroutine WriteSpectrum_hdf5
#endif
!--------------------------------------------------------------------------------------

!======================================================================================    
end module Marpes_calc

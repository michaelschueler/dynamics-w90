module wan_fft_ham
!! Provides tools for reading/writing the Wannier Hamiltonian and for computing 
!! band structures, observables, and Berry-phase properties.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use,intrinsic :: ISO_C_binding
   use omp_lib
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
   use scitools_utils,only: str,stop_error
   use wan_hamiltonian,only: wann90_tb_t
#ifdef WITHSPFFT
   use spfft
#endif
   implicit none
   include '../units_inc.f90'
   include '../formats.h'
   include 'fftw3.f03'
!--------------------------------------------------------------------------------------
   private
   public :: wann_fft_t
!--------------------------------------------------------------------------------------
   
!--------------------------------------------------------------------------------------
   type :: wann_fft_t
      integer                                    :: nwan,nwan2
      integer                                    :: nrpts
      integer                                    :: nkx,nky,nkz,nkpts,nkpts_loc,kdim
      integer,allocatable,dimension(:)           :: ndegen
      integer,allocatable,dimension(:,:)         :: irvec
      real(dp),allocatable,dimension(:,:)        :: crvec
      complex(dp),allocatable,dimension(:)       :: HA_r
      complex(dp),allocatable,dimension(:,:)     :: gradH_R
      complex(dp),allocatable,dimension(:,:,:)   :: ham_r
      complex(dp),allocatable,dimension(:,:,:,:) :: pos_r
      integer :: nthreads = 1
#ifdef WITHSPFFT
      ! .. SpFFT .. 
      integer :: processingUnit = 1
      integer :: maxNumLocalZColumns
      type(c_ptr) :: grid = c_null_ptr
      type(c_ptr) :: transform = c_null_ptr
      type(c_ptr) :: realValuesPtr
      integer,allocatable,dimension(:)           :: indices
      real(C_DOUBLE),allocatable,dimension(:)    :: spaceDomain
#else
      ! .. FFTW ..
      integer :: nx,ny,nz
      integer(kind=8) :: plan_fw,plan_bw
      complex(dp),allocatable,dimension(:)       :: work_r,work_k   
      ! .. parallelization ..
      logical :: fftw_omp
      integer :: nthreads,nthreads_fft,nthreads_orb
#endif
   contains
      procedure, public   :: InitFromW90
      procedure, public   :: Clean
      procedure, public   :: GetHam
      procedure, public   :: GetGradHam
      procedure, public   :: GetDipole
      procedure, private  :: ApplyPhaseFactor
   end type wann_fft_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
#ifdef WITHSPFFT
   include 'wan_spfft_ham_inc.F90'
#else
   include 'wan_fftw_ham_inc.F90'
#endif
!--------------------------------------------------------------------------------------
   subroutine GetGradiant(crvec,OO_R,vec_R)
      real(dp),intent(in) :: crvec(:,:)
      complex(dp),intent(in) :: OO_R(:)      
      complex(dp),intent(inout) :: vec_R(:,:)
      integer :: nrpts,ir

      nrpts = size(crvec, dim=1)
      vec_R(1:nrpts,1) = iu * crvec(1:nrpts,1) * OO_R(1:nrpts)
      vec_R(1:nrpts,2) = iu * crvec(1:nrpts,2) * OO_R(1:nrpts)
      vec_R(1:nrpts,3) = iu * crvec(1:nrpts,3) * OO_R(1:nrpts)
   
   end subroutine GetGradiant
!--------------------------------------------------------------------------------------


!======================================================================================

end module wan_fft_ham
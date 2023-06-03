module wan_fft_ham
!! Provides tools for reading/writing the Wannier Hamiltonian and for computing 
!! band structures, observables, and Berry-phase properties.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
   use wan_hamiltonian,only: wann90_tb_t
   implicit none
   include '../units_inc.f90'
   include '../formats.h'
   include '../fftw3.f03'
!--------------------------------------------------------------------------------------
   private
   public :: wann90_tb_t,ReadTB_from_w90,utility_recip_lattice
!--------------------------------------------------------------------------------------
   
!--------------------------------------------------------------------------------------
   type :: wann_fft_t
      integer                                    :: nwan
      integer                                    :: nx,ny,nz,nrpts
      integer,allocatable,dimension(:)           :: ndegen
      complex(dp),allocatable,dimension(:,:,:)   :: ham_r
      complex(dp),allocatable,dimension(:,:,:,:) :: pos_r
   contains
      procedure,public  :: InitFromW90
      procedure,public  :: Clean
   end type wann_fft_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine InitFromW90(me,w90)
      class(wann_fft_t) :: me
      type(wann90_tb_t),intent(in) :: w90
      integer :: ir,ix,iy,iz
      integer,allocatable :: iix(:),iiy(:),iiz(:)

      me%nrpts = w90%nrpts
      allocate(me%ndegen(me%nrpts))
      me%ndegen = w90%ndegen

      me%nwan = w90%num_wann

      allocate( iix(minval(w90%irvec(:,1)) : maxval(w90%irvec(:,1)) ) )
      allocate( iiy(minval(w90%irvec(:,2)) : maxval(w90%irvec(:,2)) ) )
      allocate( iiz(minval(w90%irvec(:,3)) : maxval(w90%irvec(:,3)) ) )
      me%nx = size(iix)
      me%ny = size(iiy)
      me%nz = size(iiz)

      do ir=1,me%nrpts
         w90%irvec(ir,1)
      end do


   end subroutine InitFromW90

!======================================================================================

end module wan_fft_ham
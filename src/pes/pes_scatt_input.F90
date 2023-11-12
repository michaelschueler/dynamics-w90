module pes_scatt_input
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
   use scitools_utils,only: stop_error
#ifdef WITHHDF5
   use scitools_hdf5_utils
#endif
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: scatt_input_t

   type :: scatt_input_t
      integer :: norb,lmax
      integer :: nE
      real(dp),allocatable,dimension(:)     :: Ex
      real(dp),allocatable,dimension(:,:,:) :: radint, phase
   contains
      procedure, public  :: ReadFromFile
      procedure, private :: ReadFromFile_txt
#ifdef WITHHDF5
      procedure, private :: ReadFromFile_hdf5
#endif
   end type scatt_input_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine ReadFromFile(me,fname)
      class(scatt_input_t) :: me
      character(len=*),intent(in) :: fname
      logical :: file_ok

      inquire(file=trim(fname),exist=file_ok)
      if(.not.file_ok) then
         call stop_error('Input file does not exist: '//trim(fname))
      end if

#ifdef WITHHDF5
      call me%ReadFromFile_hdf5(fname)
#else
      call me%ReadFromFile_txt(fname)
#endif

   end subroutine ReadFromFile
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine ReadFromFile_hdf5(me,fname)
      class(scatt_input_t) :: me
      character(len=*),intent(in) :: fname
      integer :: iorb
      integer(HID_t) :: file_id

      call hdf_open_file(file_id, trim(fname), STATUS='OLD', ACTION='READ')

      call hdf_read_attribute(file_id,'','ne', me%nE)
      call hdf_read_attribute(file_id,'','norb', me%norb)
      call hdf_read_attribute(file_id,'','lmax', me%lmax)

      allocate(me%Ex(me%nE))
      allocate(me%radint(me%nE,0:me%lmax,me%norb))
      allocate(me%Phase(me%nE,0:me%lmax,me%norb))

      call hdf_read_dataset(file_id,'energy',me%Ex)
      call hdf_read_dataset(file_id,'radint',me%radint)
      call hdf_read_dataset(file_id,'phase',me%Phase)

      call hdf_close_file(file_id) 

   end subroutine ReadFromFile_hdf5
#endif
!--------------------------------------------------------------------------------------
   subroutine ReadFromFile_txt(me,fname)
      class(scatt_input_t) :: me
      character(len=*),intent(in) :: fname
      integer :: iE, l, iorb
      integer :: unit_inp

      open(newunit=unit_inp, file=trim(fname), status='old', action='read')

      read(unit_inp,*) me%nE, me%norb, me%lmax

      allocate(me%Ex(me%nE))
      allocate(me%radint(me%nE,0:me%lmax,me%norb))
      allocate(me%Phase(me%nE,0:me%lmax,me%norb))      

      read(unit_inp, *) ! comment line: ## energy grid

      do iE=1,me%nE
         read(unit_inp,*) me%Ex(iE)
      end do

      read(unit_inp, *) ! comment line: ## integrals + phase

      do iorb=1,me%norb
         do iE=1,me%nE
            read(unit_inp,*) me%radint(iE,0:me%lmax,iorb), me%phase(iE,0:me%lmax,iorb)
         end do
      end do

      close(unit_inp)

   end subroutine ReadFromFile_txt
!--------------------------------------------------------------------------------------


!======================================================================================    
end module pes_scatt_input
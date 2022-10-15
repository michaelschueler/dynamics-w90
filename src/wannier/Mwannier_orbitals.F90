module Mwannier_orbitals
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
#ifdef WITHHDF5
   use scitools_hdf5_utils
#endif
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: wannier_orbs_t
!--------------------------------------------------------------------------------------
   integer,parameter :: wf_slater=0, wf_grid=1
!--------------------------------------------------------------------------------------
   type wannier_orbs_t
      logical :: real_lm=.true.
      integer :: wf_type
      integer :: natoms,norb,nr
      integer,allocatable,dimension(:)    :: atom_indx
      integer,allocatable,dimension(:)    :: L_indx,M_indx,N_indx
      real(dp),allocatable,dimension(:)   :: weight
      real(dp),allocatable,dimension(:)   :: rs,Zorb,Zscatt
      real(dp),allocatable,dimension(:,:) :: Rrad
   contains
      procedure,public  :: ReadFromTXT
#ifdef WITHHDF5
      procedure,public  :: ReadFromHDF5
#endif
   end type wannier_orbs_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine ReadFromTXT(me,fname)
      class(wannier_orbs_t) :: me
      character(len=*),intent(in) :: fname
      integer(HID_t) :: file_id
      integer :: real_flag
      integer :: unit_inp
      integer :: ir

      open(newunit=unit_inp, file=trim(fname), status='OLD', ACTION='READ')
      read(unit_inp,*) real_flag
      me%real_lm = (real_flag == 1)

      read(unit_inp,*) me%wf_type, me%natoms, me%norb

      allocate(me%atom_indx(me%norb))
      read(unit_inp,*) me%atom_indx

      allocate(me%L_indx(me%norb),me%M_indx(me%norb))
      read(unit_inp,*) me%L_indx
      read(unit_inp,*) me%M_indx           

      allocate(me%weight(me%norb))
      read(unit_inp,*) me%weight

      select case(me%wf_type)
      case(wf_slater)
         allocate(me%N_indx(me%norb),me%Zorb(me%norb),me%Zscatt(me%norb))
         read(unit_inp,*) me%N_indx
         read(unit_inp,*) me%Zorb
         read(unit_inp,*) me%Zscatt
      case(wf_grid) 
         read(unit_inp,*) me%nr         
         allocate(me%rs(me%nr),me%Rrad(me%nr,me%norb))
         do ir=1,me%nr
            read(unit_inp,*) me%rs(ir), me%Rrad(ir,1:me%norb)
         end do
      case default
         write(error_unit,fmt900) 'Invalid value for wf_type'
         stop
      end select      

      close(unit_inp)

   end subroutine ReadFromTXT
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine ReadFromHDF5(me,fname)
      class(wannier_orbs_t)       :: me
      character(len=*),intent(in) :: fname
      integer(HID_t) :: file_id
      integer :: real_flag

      call hdf_open_file(file_id, trim(fname), STATUS='OLD', ACTION='READ')

      call hdf_read_attribute(file_id,'','real_lm', real_flag)
      me%real_lm = (real_flag == 1)

      call hdf_read_attribute(file_id,'','wf_type', me%wf_type)
      call hdf_read_attribute(file_id,'','natoms', me%natoms)
      call hdf_read_attribute(file_id,'','norb', me%norb)
  

      allocate(me%atom_indx(me%norb))
      call hdf_read_dataset(file_id,'atom_indx',me%atom_indx)

      allocate(me%L_indx(me%norb),me%M_indx(me%norb))
      call hdf_read_dataset(file_id,'l_indx',me%L_indx)
      call hdf_read_dataset(file_id,'m_indx',me%M_indx)

      allocate(me%weight(me%norb))
      call hdf_read_dataset(file_id,'weight',me%weight)

      allocate(me%Zscatt(me%norb))
      call hdf_read_dataset(file_id,'zscatt',me%Zscatt)

      select case(me%wf_type)
      case(wf_slater)
         allocate(me%N_indx(me%norb),me%Zorb(me%norb))
         call hdf_read_dataset(file_id,'n_indx',me%N_indx)
         call hdf_read_dataset(file_id,'zorb',me%Zorb)
      case(wf_grid) 
         call hdf_read_attribute(file_id,'','nr', me%nr)
         allocate(me%rs(me%nr),me%Rrad(me%nr,me%norb))
         call hdf_read_dataset(file_id,'rs',me%rs)
         call hdf_read_dataset(file_id,'rrad',me%Rrad)
      case default
         write(error_unit,fmt900) 'Invalid value for wf_type'
         stop
      end select

      call hdf_close_file(file_id) 

   end subroutine ReadFromHDF5
#endif
!--------------------------------------------------------------------------------------


!======================================================================================    
end module Mwannier_orbitals
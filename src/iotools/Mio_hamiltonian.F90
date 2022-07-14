module Mio_hamiltonian
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero
   use Mutils,only: get_file_ext, check_file_ext
   use Mham_w90,only: wann90_tb_t
   use Mwann_soc,only: ham_soc_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------  
   private
   public :: ReadHamiltonian, WriteHamiltonian, Read_SOC_Hamiltonian, Read_SOC_lambda
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine ReadHamiltonian(file_wann,wann)
      character(len=*),intent(in)   :: file_wann
      type(wann90_tb_t),intent(out) :: wann

      if(check_file_ext(file_wann, "tb") .or. check_file_ext(file_wann, "dat")) then
         call wann%ReadFromW90(file_wann)
      elseif(check_file_ext(file_wann, "h5")) then
#ifdef WITHHDF5
         call wann%ReadFromHDF5(file_wann)
#else
      write(error_unit,fmt900) "No HDF5 support. Can't read "//trim(file_wann)
      stop         
#endif
      else
         write(error_unit,fmt900) "Unrecognized file extension: "//get_file_ext(file_wann)
         stop
      end if     

   end subroutine ReadHamiltonian
!--------------------------------------------------------------------------------------
   subroutine WriteHamiltonian(wann,file_wann)
      type(wann90_tb_t),intent(in)  :: wann
      character(len=*),intent(in)   :: file_wann

      if(check_file_ext(file_wann, "tb") .or. check_file_ext(file_wann, "dat")) then
         call wann%SaveToW90(file_wann)
      elseif(check_file_ext(file_wann, "h5")) then
#ifdef WITHHDF5
         call wann%SaveToHDF5(file_wann)
#else
      write(error_unit,fmt900) "No HDF5 support. Can't save "//trim(file_wann)
      stop         
#endif
      else
         write(error_unit,fmt900) "Unrecognized file extension: "//get_file_ext(file_wann)
         stop
      end if         

   end subroutine WriteHamiltonian
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
   subroutine Read_SOC_Hamiltonian(fname,soc)
      character(len=*),intent(in)   :: fname
      type(ham_soc_t),intent(out)   :: soc

      if(check_file_ext(fname, "h5")) then
#ifdef WITHHDF5
         call soc%ReadFromHDF5(fname)
#else
         write(error_unit,fmt900) "No HDF5 support. Can't read "//trim(fname)
         stop         
#endif
      else 
         call soc%ReadFromTXT(fname)
      end if
   end subroutine Read_SOC_Hamiltonian
!--------------------------------------------------------------------------------------



!======================================================================================
end module Mio_hamiltonian
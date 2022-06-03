module Mio_hamiltonian
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,iu,zero
   use Mutils,only: get_file_ext, check_file_ext
   use Mham_w90,only: wann90_tb_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------  
   private
   public :: ReadHamiltonian, WriteHamiltonian
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
!======================================================================================
end module Mio_hamiltonian
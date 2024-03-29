module io_orbitals
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero
   use scitools_utils,only: get_file_ext, check_file_ext
   use wan_orbitals,only: wannier_orbs_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------  
   private
   public :: ReadWannierOrbitals
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine ReadWannierOrbitals(fname,orbs)
      character(len=*),intent(in)            :: fname
      type(wannier_orbs_t),intent(out)       :: orbs

      if(check_file_ext(fname, "txt") .or. check_file_ext(fname, "dat")) then
         call orbs%ReadFromTXT(fname)
      elseif(check_file_ext(fname, "h5")) then
#ifdef WITHHDF5
         call orbs%ReadFromHDF5(fname)
#else
      write(error_unit,fmt900) "No HDF5 support. Can't read "//trim(fname)
      stop         
#endif
      else
         write(error_unit,fmt900) "Unrecognized file extension: "//get_file_ext(fname)
         stop
      end if     

   end subroutine ReadWannierOrbitals
!--------------------------------------------------------------------------------------



!======================================================================================
end module io_orbitals
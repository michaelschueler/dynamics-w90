module io_hamiltonian
!! Provides generic interfaces for in/output the Wannier Hamiltonian in `Wannier90`
!! or hdf5 format. Provides also tools for reading spin-orbit coupling parameters from file.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero
   use scitools_utils,only: get_file_ext, check_file_ext
   use wan_hamiltonian,only: wann90_tb_t
   use wan_soc,only: ham_soc_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------  
   private
   public :: ReadHamiltonian, WriteHamiltonian, Read_SOC_Hamiltonian, Read_SOC_lambda
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine ReadHamiltonian(file_wann,wann,file_xyz)
      !! Reads the Wannier Hamiltonian from file. If the file extension is `tb` or `dat`,
      !! we assume the Hamiltonian has been produced by `Wannier90` (`*_tb.dat` file).
      !! If the file extension is `h5` we read from hdf5 if `dynamics-w90` has been
      !! build with hdf5 support.
      !! Optionally the coordinates of the Wannier centeres can be read from file,
      !! following the format of the `*_centres.xyz` of `Wannier90`.
      character(len=*),intent(in)            :: file_wann !! file name for the Hamiltonian
      type(wann90_tb_t),intent(out)          :: wann !! Wannier class to be initialized from the input
      character(len=*),intent(in),optional   :: file_xyz !! file name for the coordinates

      if(check_file_ext(file_wann, "tb") .or. check_file_ext(file_wann, "dat")) then
         if(present(file_xyz)) then
            call wann%ReadFromW90(file_wann,file_xyz=file_xyz)
         else
            call wann%ReadFromW90(file_wann)
         end if
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
      !! Writes the Hamiltonian stored in the Wannier class `wann` to file.
      !! If the file extension of `file_wann` is `tb` or `dat` a text file in 
      !! `Wannier90` format will be produced. If the file extension is `h5`, the 
      !! Hamiltonian will be written in hdf5 format if `dynamics-w90` has been
      !! build with hdf5 support.
      type(wann90_tb_t),intent(in)  :: wann !! Wannier class
      character(len=*),intent(in)   :: file_wann !! file name for the output

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
      !! Reads the spin-orbit coupling (SOC) constants from a simple text file. We assume
      !! complete shells of orbitals (s, p, d, ...) for the orbitals of the Wannier Hamiltonian.
      !! For each group and for each atom we assign a SOC constant.
      character(len=*),intent(in) :: fname !! file name with the constants
      real(dp),allocatable,intent(out) :: lam(:) !! SOC constants for each group of orbitals
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
      !! Reads the spin-orbit coupling (SOC) Hamiltonioan from file. The SOC Hamiltonian
      !! is given by the matrix representation of the angular momentum operator for each
      !! shell. The script `SOCInput_txt.py` / `SOCInput.py` prepares such a file in 
      !! plain text / hdf5 format.
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
end module io_hamiltonian
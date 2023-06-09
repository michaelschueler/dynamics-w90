module io_density
!! Provides generic interfaces for in/output the Wannier Hamiltonian in `Wannier90`
!! or hdf5 format. Provides also tools for reading spin-orbit coupling parameters from file.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero
   use scitools_utils,only: get_file_ext, check_file_ext, stop_error
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------  
   private
   public :: ReadDensityMatrix, WriteDensityMatrix
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine ReadDensityMatrix(file_dens,Rhok,tstart)
      !! Reads the density matrix from file. If the file extension is `bin``,
      !! we assume plain binary format.
      !! If the file extension is `h5` we read from hdf5 if `dynamics-w90` has been
      !! build with hdf5 support.
      character(len=*),intent(in)            :: file_dens !! file name for the density matrix
      complex(dp),allocatable,intent(out)    :: Rhok(:,:,:)
      real(dp),intent(out)                   :: tstart

      if(check_file_ext(file_dens, "bin")) then
         call ReadDensityMatrix_binary(file_dens,Rhok,tstart)
      elseif(check_file_ext(file_dens, "h5")) then
#ifdef WITHHDF5
         call ReadDensityMatrix_hdf5(file_dens,Rhok,tstart)
#else
      call stop_error("No HDF5 support. Can't read "//trim(file_dens))
#endif
      else
         call stop_error("Unrecognized file extension: "//get_file_ext(file_dens))
      end if     

   end subroutine ReadDensityMatrix
!--------------------------------------------------------------------------------------
   subroutine WriteDensityMatrix(prefix_dens,Rhok,tstop)
      !! Writes the density matrix to file.
      !! If the file extension of `file_dens` is `bin` a binary file in will be produced. 
      !! If the file extension is `h5`, the 
      !! density matrix will be written in hdf5 format if `dynamics-w90` has been
      !! build with hdf5 support.
      character(len=*),intent(in)            :: prefix_dens !! file name for the density matrix
      complex(dp),intent(in)                 :: Rhok(:,:,:)
      real(dp),intent(in)                    :: tstop
      character(len=256) :: file_dens

      file_dens = trim(prefix_dens) // "_dens."
#ifdef WITHHDF5
      file_dens = trim(file_dens) // "h5"
      call WriteDensityMatrix_hdf5(file_dens,Rhok,tstop)
#else
      file_dens = trim(file_dens) // "bin"
      call WriteDensityMatrix_binary(file_dens,Rhok,tstop)
#endif

   end subroutine WriteDensityMatrix
!--------------------------------------------------------------------------------------
   subroutine ReadDensityMatrix_binary(file_dens,Rhok,tstart)
      character(len=*),intent(in)            :: file_dens !! file name for the density matrix
      complex(dp),allocatable,intent(out)    :: Rhok(:,:,:)
      real(dp),intent(out)                   :: tstart
      integer :: nbnd,nk
      integer :: unit_inp

      open(newunit=unit_inp, file=trim(file_dens), form='unformatted' , status='old', action='read')

      read(unit_inp) nbnd , nk
      read(unit_inp) tstart
      allocate(Rhok(nbnd,nbnd,nk))
      read(unit_inp) Rhok

      close(unit_inp)

   end subroutine ReadDensityMatrix_binary
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine ReadDensityMatrix_hdf5(file_dens,Rhok,tstart)
      use scitools_hdf5_utils
      character(len=*),intent(in)            :: file_dens !! file name for the density matrix
      complex(dp),allocatable,intent(out)    :: Rhok(:,:,:)
      real(dp),intent(out)                   :: tstart
      integer(HID_t) :: file_id
      integer :: nbnd,nk,nbnd_,nk_
      integer :: unit_inp
      real(dp),allocatable :: rdata(:,:,:)

      call hdf_open_file(file_id, trim(file_dens), STATUS='OLD', ACTION='READ')

      call hdf_read_attribute(file_id,'','nbnd', nbnd)
      call hdf_read_attribute(file_id,'','nk', nk)
      call hdf_read_attribute(file_id,'','tstart', tstart)

      allocate(Rhok(nbnd,nbnd,nk))
      allocate(rdata(nbnd,nbnd,nk))

      call hdf_read_dataset(file_id,'real',rdata)
      Rhok = rdata
      call hdf_read_dataset(file_id,'imag',rdata)
      Rhok = Rhok + iu * rdata

      deallocate(rdata)

      call hdf_close_file(file_id)

   end subroutine ReadDensityMatrix_hdf5
#endif
!--------------------------------------------------------------------------------------
   subroutine WriteDensityMatrix_binary(file_dens,Rhok,tstop)
      character(len=*),intent(in)            :: file_dens !! file name for the density matrix
      complex(dp),intent(in)                 :: Rhok(:,:,:)
      real(dp),intent(in)                    :: tstop
      integer :: nbnd,nk
      integer :: unit_out

      nbnd = size(Rhok, dim=1)
      nk = size(Rhok, dim=3)

      open(newunit=unit_out, file=trim(file_dens), form='unformatted' , status='replace')

      write(unit_out) nbnd, nk
      write(unit_out) tstop
      write(unit_out) Rhok

      close(unit_out)

   end subroutine WriteDensityMatrix_binary
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine WriteDensityMatrix_hdf5(file_dens,Rhok,tstop)
      use scitools_hdf5_utils
      character(len=*),intent(in)            :: file_dens !! file name for the density matrix
      complex(dp),intent(in)                 :: Rhok(:,:,:)
      real(dp),intent(in)                    :: tstop
      integer :: nbnd,nk
      integer(HID_t) :: file_id
      real(dp),allocatable :: rdata(:,:,:)

      nbnd = size(Rhok, dim=1)
      nk = size(Rhok, dim=3)

      call hdf_open_file(file_id, trim(file_dens), STATUS='NEW')

      call hdf_write_attribute(file_id,'','nbnd', nbnd)
      call hdf_write_attribute(file_id,'','nk', nk)
      call hdf_write_attribute(file_id,'','tstart', tstop)

      allocate(rdata(nbnd,nbnd,nk))
      rdata = dble(Rhok)
      call hdf_write_dataset(file_id,'real',rdata)
      rdata = aimag(Rhok)
      call hdf_write_dataset(file_id,'imag',rdata)
      deallocate(rdata)

      call hdf_close_file(file_id)

   end subroutine WriteDensityMatrix_hdf5
#endif
!--------------------------------------------------------------------------------------

!======================================================================================
end module io_density
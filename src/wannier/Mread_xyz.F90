module Mread_xyz
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use scitools_def,only: dp
   use scitools_utils,only: checkoption,loadtxt,str
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: ReadXYZ                                                                                                                                                                                                               
!--------------------------------------------------------------------------------------
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine ReadXYZ(fname,coords)
   !! Reads the coordinates of the Wannier centeres from the `*_centres.xyz` file
   !! produced by `Wannier90`.
      character(len=*),intent(in) :: fname !! file name of the `*_centres.xyz` file
      real(dp),intent(inout)      :: coords(:,:) !! coordinates of the Wannier centers
      integer :: num_wann,num_pos,i,iorb
      character(len=2) :: sx
      integer :: unit_inp
      real(dp) :: rvec(3)

      open(newunit=unit_inp, file=trim(fname), status='old', action='read')
      read(unit_inp,*) num_pos

      read(unit_inp,*) ! comment line

      iorb = 0
      do i=1,num_pos
         read(unit_inp,*) sx,rvec
         if(trim(sx) == "X") then
            iorb = iorb + 1
         end if
      end do

      if(size(coords,1) /= iorb) then
         write(error_unit,fmt900) "wrong number of coordinates"
         close(unit_inp)
         stop
      end if

      rewind(unit_inp)

      read(unit_inp,*) num_pos
      read(unit_inp,*) ! comment line

      iorb = 0
      do i=1,num_pos
         read(unit_inp,*) sx,rvec
         if(trim(sx) == "X") then
            iorb = iorb + 1
            coords(iorb,1:3) = rvec
         end if
      end do

      close(unit_inp)

   end subroutine ReadXYZ
!--------------------------------------------------------------------------------------

!======================================================================================   
end module Mread_xyz
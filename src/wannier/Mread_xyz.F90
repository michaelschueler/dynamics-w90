module Mread_xyz
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp
   use Mutils,only: checkoption,loadtxt,str
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
      character(len=*),intent(in) :: fname
      real(dp),intent(inout)      :: coords(:,:)
      integer :: num_wann,i
      character(len=1) :: sx
      integer :: unit_inp

      open(newunit=unit_inp, file=trim(fname), status='old', action='read')
      read(unit_inp,*) num_wann

      if(size(coords,1) /= num_wann) then
         write(error_unit,fmt900) "wrong number of coordinates"
         close(unit_inp)
         stop
      end if

      read(unit_inp,*) ! comment line

      do i=1,num_wann
         read(unit_inp,*) sx,coords(i,:)
      end do

      close(unit_inp)

   end subroutine ReadXYZ
!--------------------------------------------------------------------------------------

!======================================================================================   
end module Mread_xyz
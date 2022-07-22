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
module Mlatt_kpts
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp
   use Mutils,only: checkoption,loadtxt,str
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: Read_Kpoints,GenKgrid                                                                                                                                                                                                                    
!--------------------------------------------------------------------------------------
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Read_Kpoints(fname,kpts,print_info)
      character(len=*),intent(in) :: fname
      real(dp),allocatable,intent(out) :: kpts(:,:)
      logical,intent(in),optional :: print_info
      logical :: info_
      character(len=32)  :: kpoints_type="grid"
      character(len=255) :: file_kpts=""
      integer :: unit_inp, unit_k
      integer :: nk1=1,nk2=1,nk3=1
      integer :: nk
      namelist/KPOINTS/kpoints_type,file_kpts,nk1,nk2,nk3

      info_ = .true.
      if(present(print_info)) info_ = print_info

      open(newunit=unit_inp, file=trim(fname), status='old', action='read')
      read(unit_inp, nml=KPOINTS)
      close(unit_inp)

      if(checkoption(kpoints_type, "list")) then
         call loadtxt(file_kpts,kpts)
      elseif(checkoption(kpoints_type, "path")) then
         call ReadKpath(file_kpts, kpts)
      elseif(checkoption(kpoints_type, "grid")) then
         call GenKgrid(kpts,nk1,nk2,nk3)
      else
         write(error_unit,fmt900) "invalid k-points format"
         stop
      end if

      nk = size(kpts,1)

      if(info_) then
         if(nk == 1) then
            write(output_unit,fmt_info) str(size(kpts,1))//" k-point was read from file."
         else
            write(output_unit,fmt_info) str(size(kpts,1))//" k-points were read from file."
         end if
         write(output_unit,*)
      end if
         
   end subroutine Read_Kpoints
!--------------------------------------------------------------------------------------
   subroutine ReadKpath(Flname,kpts)
      character(len=*),intent(in)      :: Flname
      real(dp),allocatable,intent(out) :: kpts(:,:)
      logical  :: file_there
      integer  :: iunit
      integer  :: num_points,Nk,i,j,k
      real(dp) :: x
      integer,allocatable  :: num_seg_points(:)
      real(dp),allocatable :: points(:,:)

      inquire(file=trim(Flname),exist=file_there)

      if(.not.file_there) then
         write(output_unit,fmt900) 'Input file does not exist: '//trim(Flname)
         stop
      end if

      open(newunit=iunit,file=trim(Flname),status='OLD',action='READ')
      read(iunit,*) num_points
      allocate(num_seg_points(num_points-1))
      read(iunit,*) num_seg_points
      allocate(points(3,num_points))
      do i=1,num_points
         read(iunit,*) points(:,i)
      end do

      Nk = sum(num_seg_points)+1
      allocate(kpts(Nk,3))

      k = 0
      do i=1,num_points-1
         do j=1,num_seg_points(i)
            k = k + 1
            x = (j-1)/dble(num_seg_points(i))
            kpts(k,:) = (1.0_dp-x) * points(:,i) + x * points(:,i+1)
         end do
      end do

      kpts(Nk,:) = points(:,num_points)

      deallocate(num_seg_points,points)

      close(iunit)

   end subroutine ReadKpath
!--------------------------------------------------------------------------------------
   subroutine GenKgrid(kpts,nk1,nk2,nk3)
      real(dp),allocatable,intent(out) :: kpts(:,:)
      integer,intent(in) :: nk1,nk2
      integer,intent(in),optional :: nk3
      integer :: nk3_
      integer :: nk,i1,i2,i3,ik

      nk3_ = 1
      if(present(nk3)) nk3_ = nk3

      nk = nk1 * nk2 * nk3_
      allocate(kpts(nk,3)); kpts = 0.0_dp

      if(nk3_ == 1) then
         ik = 0
         do i1=1,nk1
            do i2=1,nk2
               ik = ik + 1
               kpts(ik,1) = -0.5_dp + (i1-1)/dble(nk1)
               kpts(ik,2) = -0.5_dp + (i2-1)/dble(nk2)
            end do
         end do
      else 
         ik = 0
         do i1=1,nk1
            do i2=1,nk2
               do i3=1,nk3
                  ik = ik + 1
                  kpts(ik,1) = -0.5_dp + (i1-1)/dble(nk1)
                  kpts(ik,2) = -0.5_dp + (i2-1)/dble(nk2)
                  kpts(ik,3) = -0.5_dp + (i3-1)/dble(nk3)
               end do
            end do
         end do
      end if

   end subroutine GenKgrid
!--------------------------------------------------------------------------------------

!======================================================================================   
end module Mlatt_kpts
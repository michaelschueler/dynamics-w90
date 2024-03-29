module wan_latt_kpts
!! Provides tools to read k-points from file.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use scitools_def,only: dp
   use scitools_utils,only: checkoption,loadtxt,str
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
   subroutine Read_Kpoints(fname,kpts,print_info,root_tag)
   !! Reads k-points from file. There are three formats of specifying k-points:
   !!
   !! * A row-by-row list of points. Triggered by `kpoints_type="list"`.
   !! * A path specified by the special points and the number of points in between.
   !!   Triggerd by `kpoints_type="path"`.
   !! * A grid covering the full Brillouin zone. Triggerd by `kpoints_type="grid"`.
   !! 
   !! If `kpoints_type="list"` the code expects a list of k-points with kx (1st column),
   !! ky (2nd column), and kz (3d column). Input in reduced coordinates is expected.
   !! If `kpoints_type="path"` a path in k-space will be constructed from the input. 
   !! The `file_kpts` file has the following format:
   !! ```
   !!    npoints ndim
   !!    nseg1 nseg2 ...
   !!    point1 
   !!    point2
   !!    ..
   !! ```
   !! Here, `npoints` is the number of points to pass through, while `nseg1` is the number 
   !! of segments between `point1` and `point2` and so on. Below the special points are listed.
   !! 
   !! 
   !! For `kpoints_type="grid"`, the number of points along the reciprocal lattice directions 
   !! `nk1`, `nk2`, `nk3` is read from the `KPOINTS` name list.
      character(len=*),intent(in) :: fname !! name of the input file containing the `KPOINTS` name list
      real(dp),allocatable,intent(out) :: kpts(:,:) !! the k-points read from file
      logical,intent(in),optional :: print_info !! if `.true.`, some basic info about the k-points is written
      logical,intent(in),optional :: root_tag !! prints the k-point info only if `root_tag=.true.`. 
                                              !! Useful for calls in an MPI program.
      logical :: info_,root_
      character(len=32)  :: kpoints_type="grid"
      character(len=255) :: file_kpts=""
      integer :: unit_inp, unit_k
      integer :: nk1=1,nk2=1,nk3=1
      integer :: nk
      namelist/KPOINTS/kpoints_type,file_kpts,nk1,nk2,nk3

      info_ = .true.
      if(present(print_info)) info_ = print_info

      root_ = .true.
      if(present(root_tag)) root_ = root_tag

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

      if(info_ .and. root_) then
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
   !! Generates a grid covering either the 2D or the 3D Brillouin zone by the
   !! Monkhorst-Pack scheme. 
      real(dp),allocatable,intent(out) :: kpts(:,:) !! the k-points
      integer,intent(in) :: nk1,nk2 !! number of points along the first two reciprocal lattice directions
      integer,intent(in),optional :: nk3 !! Number of points along the third reciprocal lattice direction.
                                         !! If not specified, we assume a 2D Brillouin zone.
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
end module wan_latt_kpts
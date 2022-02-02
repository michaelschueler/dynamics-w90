module Mlatt
!====================================================================================== 
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdef,only: dp
   implicit none
!--------------------------------------------------------------------------------------
   private
   public &
      latt3d_t
!--------------------------------------------------------------------------------------
   type latt3d_t
      integer    :: Nk1,Nk2,Nk3,Nk !! number of k-points in each direction
      integer,allocatable,dimension(:,:)  :: kindex
      real(dp),allocatable,dimension(:,:) :: kcoord(:,:),kpts(:,:)
   contains
      procedure,public :: Init => latt3d_Init
      procedure,public :: Clean => latt3d_Clean
   end type latt3d_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine latt3d_Init(me,Nk1,Nk2,Nk3,recip_lat)
      class(latt3d_t)     :: me
      integer,intent(in)  :: Nk1
      integer,intent(in)  :: Nk2
      integer,intent(in)  :: Nk3
      real(dp),intent(in) :: recip_lat(3,3)
      integer  :: i1,i2,i3,ik

      me%Nk1 = Nk1
      me%Nk2 = Nk2
      me%Nk3 = Nk3
      me%Nk = me%Nk1 * me%Nk2 * me%Nk3

      allocate(me%kindex(me%Nk,3),me%kcoord(me%Nk,3),me%kpts(me%Nk,3))

      if(me%Nk3 > 1) then
         ik = 0
         do i1=1,me%Nk1
            do i2=1,me%Nk2
               do i3=1,me%Nk3
                  ik = ik + 1
                  me%kcoord(ik,1) = (i1-1)/dble(me%Nk1) - 0.5_dp
                  me%kcoord(ik,2) = (i2-1)/dble(me%Nk2) - 0.5_dp
                  me%kcoord(ik,3) = (i3-1)/dble(me%Nk3) - 0.5_dp
                  me%kindex(ik,1:3) = [i1,i2,i3]
               end do
            end do
         end do
      else
         me%kcoord = 0.0_dp
         ik = 0
         do i1=1,me%Nk1
            do i2=1,me%Nk2
               ik = ik + 1
               me%kcoord(ik,1) = (i1-1)/dble(me%Nk1) - 0.5_dp
               me%kcoord(ik,2) = (i2-1)/dble(me%Nk2) - 0.5_dp
               me%kindex(ik,1:3) = [i1,i2,1]
            end do
         end do
      end if

      do ik=1,me%Nk
         me%kpts(ik,1) = me%kcoord(ik,1) * recip_lat(1,1) &
            + me%kcoord(ik,2) * recip_lat(2,1) &
            + me%kcoord(ik,3) * recip_lat(3,1)
         me%kpts(ik,2) = me%kcoord(ik,1) * recip_lat(1,2) &
            + me%kcoord(ik,2) * recip_lat(2,2) &
            + me%kcoord(ik,3) * recip_lat(3,2)
         me%kpts(ik,3) = me%kcoord(ik,1) * recip_lat(1,3) &
            + me%kcoord(ik,2) * recip_lat(2,3) &
            + me%kcoord(ik,3) * recip_lat(3,3)
      end do

   end subroutine latt3d_Init
!--------------------------------------------------------------------------------------
   subroutine latt3d_Clean(me)
      class(latt3d_t) :: me

      deallocate(me%kcoord)
      deallocate(me%kindex)
      deallocate(me%kpts)

   end subroutine latt3d_Clean
!--------------------------------------------------------------------------------------


!======================================================================================   
end module Mlatt

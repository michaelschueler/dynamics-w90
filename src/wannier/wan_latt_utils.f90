module wan_latt_utils
!======================================================================================
   use scitools_def,only: dp
   use scitools_linalg,only: Inv
   implicit none
   include '../units_inc.f90'
!--------------------------------------------------------------------------------------   
   private
   public :: utility_recip_lattice, utility_recip_reduced
   public :: utility_Cart2Red_2D, utility_Cart2Red_3D
!-------------------------------------------------------------------------------------- 
contains
!-------------------------------------------------------------------------------------- 
   subroutine utility_recip_lattice(real_lat, recip_lat, volume)  
      !==================================================================!
      !                                                                  !
      !!  Calculates the reciprical lattice vectors and the cell volume
      !                                                                  !
      !===================================================================
      implicit none
      real(kind=dp), intent(in)  :: real_lat(3, 3)
      real(kind=dp), intent(out) :: recip_lat(3, 3)
      real(kind=dp), intent(out),optional :: volume
      real(kind=dp) :: vol

      recip_lat(1, 1) = real_lat(2, 2)*real_lat(3, 3) - real_lat(3, 2)*real_lat(2, 3)
      recip_lat(1, 2) = real_lat(2, 3)*real_lat(3, 1) - real_lat(3, 3)*real_lat(2, 1)
      recip_lat(1, 3) = real_lat(2, 1)*real_lat(3, 2) - real_lat(3, 1)*real_lat(2, 2)
      recip_lat(2, 1) = real_lat(3, 2)*real_lat(1, 3) - real_lat(1, 2)*real_lat(3, 3)
      recip_lat(2, 2) = real_lat(3, 3)*real_lat(1, 1) - real_lat(1, 3)*real_lat(3, 1)
      recip_lat(2, 3) = real_lat(3, 1)*real_lat(1, 2) - real_lat(1, 1)*real_lat(3, 2)
      recip_lat(3, 1) = real_lat(1, 2)*real_lat(2, 3) - real_lat(2, 2)*real_lat(1, 3)
      recip_lat(3, 2) = real_lat(1, 3)*real_lat(2, 1) - real_lat(2, 3)*real_lat(1, 1)
      recip_lat(3, 3) = real_lat(1, 1)*real_lat(2, 2) - real_lat(2, 1)*real_lat(1, 2)

      vol = real_lat(1, 1)*recip_lat(1, 1) + &
         real_lat(1, 2)*recip_lat(1, 2) + &
         real_lat(1, 3)*recip_lat(1, 3)

      recip_lat = DPI*recip_lat/vol
      vol = abs(vol)

      if(present(volume)) volume = vol

   end subroutine utility_recip_lattice
!--------------------------------------------------------------------------------------
   subroutine utility_recip_reduced(recip_lat, recip_red)
      !==================================================================!
      !                                                                  !
      !!  Calculates the matrix to convert cartesian k coordinates to
      !!  reduced coordinates
      !                                                                  !
      !===================================================================
      real(dp),intent(in)  :: recip_lat(3,3)
      real(dp),intent(out) :: recip_red(3,3)

      recip_red = transpose(Inv(recip_lat))

   end subroutine utility_recip_reduced
!--------------------------------------------------------------------------------------
   function utility_Cart2Red_2D(recip_red,kvec) result(kred)
      real(dp),intent(in) :: recip_red(:,:)  
      real(dp),intent(in) :: kvec(2)
      real(dp)            :: kred(2)

      kred = matmul(recip_red(1:2,1:2), kvec)

   end function utility_Cart2Red_2D
!--------------------------------------------------------------------------------------
   function utility_Cart2Red_3D(recip_red,kvec) result(kred)
      real(dp),intent(in) :: recip_red(:,:)  
      real(dp),intent(in) :: kvec(3)
      real(dp)            :: kred(3)

      kred = matmul(recip_red(1:3,1:3), kvec)

   end function utility_Cart2Red_3D
!--------------------------------------------------------------------------------------

!======================================================================================    
end module wan_latt_utils
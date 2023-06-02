module wan_utils
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp, zero
   use scitools_linalg,only: Inv, util_zgemm
   implicit none
   include '../units_inc.f90'
!--------------------------------------------------------------------------------------   
   private
   public :: utility_recip_lattice, utility_recip_reduced
   public :: utility_Cart2Red_2D, utility_Cart2Red_3D
   public :: utility_diagonalize, utility_rotate_diag, utility_matmul_diag
   public :: Batch_Diagonalize_t
!--------------------------------------------------------------------------------------   
   type :: Batch_Diagonalize_t
      integer :: nbnd
      real(dp),allocatable,dimension(:)             :: eps
      complex(dp),allocatable,dimension(:,:)        :: vect
      integer,allocatable,dimension(:)              :: iwork,ifail
      real(dp),allocatable,dimension(:)             :: rwork
      complex(dp),allocatable,dimension(:)          :: mat_pack, cwork
   contains
      procedure, public :: Init => Batch_Diag_Init
      procedure, public :: Clean => Batch_Diag_Clean
      procedure, public :: Diagonalize => Batch_Diag_Diagonalize
   end type Batch_Diagonalize_t
!-------------------------------------------------------------------------------------- 
contains
!-------------------------------------------------------------------------------------- 
   subroutine Batch_Diag_Init(me,nbnd)
      class(Batch_Diagonalize_t) :: me
      integer,intent(in) :: nbnd

      me%nbnd = nbnd
      allocate(me%eps(me%nbnd))
      allocate(me%vect(me%nbnd,me%nbnd))
      allocate(me%iwork(5*me%nbnd))
      allocate(me%ifail(me%nbnd))
      allocate(me%rwork(7*me%nbnd))
      allocate(me%cwork(5*me%nbnd))   
      allocate(me%mat_pack((me%nbnd*(me%nbnd + 1))/2))   

   end subroutine Batch_Diag_Init
!-------------------------------------------------------------------------------------- 
   subroutine Batch_Diag_Clean(me)
      class(Batch_Diagonalize_t) :: me

      if(allocated(me%eps)) deallocate(me%eps)
      if(allocated(me%vect)) deallocate(me%vect)
      if(allocated(me%iwork)) deallocate(me%iwork)
      if(allocated(me%ifail)) deallocate(me%ifail)
      if(allocated(me%rwork)) deallocate(me%rwork)
      if(allocated(me%cwork)) deallocate(me%cwork)
      if(allocated(me%mat_pack)) deallocate(me%mat_pack)

   end subroutine Batch_Diag_Clean
!-------------------------------------------------------------------------------------- 
   subroutine Batch_Diag_Diagonalize(me,Hk,epsk,vectk,info)
      class(Batch_Diagonalize_t) :: me
      complex(dp),intent(in)     :: Hk(:,:)
      real(dp),intent(inout),optional :: epsk(:)
      complex(dp),intent(inout),optional :: vectk(:,:)
      integer,intent(out),optional :: info
      integer :: i,j
      integer :: info_,nfound

      do j = 1, me%nbnd
         do i = 1, j
            me%mat_pack(i + ((j - 1)*j)/2) = Hk(i, j)
         end do
      end do
      me%vect = zero; me%eps = 0.0_dp; me%cwork = zero; me%rwork = 0.0_dp; me%iwork = 0
      call ZHPEVX('V', 'A', 'U', me%nbnd, me%mat_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
         nfound, me%eps(1), me%vect, me%nbnd, me%cwork, me%rwork, me%iwork, me%ifail, info_)
      if (info_ < 0) then
         write(output_unit, '(a,i3,a)') 'THE ', -info_, &
         ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
         write(error_unit,*) 'Error in utility_diagonalize'
      end if
      if (info_ > 0) then
         write(output_unit, '(i3,a)') info_, ' EIGENVECTORS FAILED TO CONVERGE'
         write(error_unit,*) 'Error in utility_diagonalize'
      end if

      if(present(epsk)) epsk = me%eps
      if(present(vectk)) vectk = me%vect
      if(present(info)) info = info_

   end subroutine Batch_Diag_Diagonalize
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
!--------------------------------------------------------------------------------------
   subroutine utility_diagonalize(mat, dim, eig, rot)
      !============================================================!
      !                                                            !
      !! Diagonalize the dim x dim  hermitian matrix 'mat' and
      !! return the eigenvalues 'eig' and the unitary rotation 'rot'
      !                                                            !
      !============================================================!
      use,intrinsic::iso_fortran_env,only: output_unit, error_unit

      integer, intent(in)           :: dim
      complex(kind=dp), intent(in)  :: mat(dim, dim)
      real(kind=dp), intent(out)    :: eig(dim)
      complex(kind=dp), intent(out) :: rot(dim, dim)

      complex(kind=dp)   :: mat_pack((dim*(dim + 1))/2), cwork(2*dim)
      real(kind=dp)      :: rwork(7*dim)
      integer            :: i, j, info, nfound, iwork(5*dim), ifail(dim)

      do j = 1, dim
         do i = 1, j
            mat_pack(i + ((j - 1)*j)/2) = mat(i, j)
         end do
      end do
      rot = zero; eig = 0.0_dp; cwork = zero; rwork = 0.0_dp; iwork = 0
      call ZHPEVX('V', 'A', 'U', dim, mat_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
         nfound, eig(1), rot, dim, cwork, rwork, iwork, ifail, info)
      if (info < 0) then
         write (output_unit, '(a,i3,a)') 'THE ', -info, &
         ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
         write(error_unit,*) 'Error in utility_diagonalize'
      end if
      if (info > 0) then
         write (output_unit, '(i3,a)') info, ' EIGENVECTORS FAILED TO CONVERGE'
         write(error_unit,*) 'Error in utility_diagonalize'
      end if

  end subroutine utility_diagonalize
!--------------------------------------------------------------------------------------
   function utility_rotate_diag(mat, rot, dim)
      !===========================================================!
      !                                                           !
      !! Rotates the dim x dim matrix 'mat' according to
      !! (rot)^dagger.mat.rot, where 'rot' is a unitary matrix.
      !! Computes only the diagonal elements of rotated matrix.
      !                                                           !
      !===========================================================!
      integer          :: dim
      complex(kind=dp) :: utility_rotate_diag(dim)
      complex(kind=dp) :: mat(dim, dim)
      complex(kind=dp) :: rot(dim, dim)
      complex(kind=dp) :: tmp(dim, dim)

      call util_zgemm(rot, mat, tmp, transa_opt='C')
      utility_rotate_diag = utility_matmul_diag(tmp, rot, dim)

   end function utility_rotate_diag
!--------------------------------------------------------------------------------------
   function utility_matmul_diag(mat1, mat2, dim)
      !===========================================================!
      !                                                           !
      !! Computes the diagonal elements of the matrix mat1.mat2
      !                                                           !
      !===========================================================!
      integer          :: dim
      complex(kind=dp) :: utility_matmul_diag(dim)
      complex(kind=dp) :: mat1(dim, dim)
      complex(kind=dp) :: mat2(dim, dim)

      integer i, j

      utility_matmul_diag = zero
      do i = 1, dim
         do j = 1, dim
            utility_matmul_diag(i) = utility_matmul_diag(i) + mat1(i, j)*mat2(j, i)
         end do
      end do

   end function utility_matmul_diag
!--------------------------------------------------------------------------------------

!======================================================================================    
end module wan_utils
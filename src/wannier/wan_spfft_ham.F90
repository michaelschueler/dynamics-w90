module wan_spfft_ham
!! Provides tools for reading/writing the Wannier Hamiltonian and for computing 
!! band structures, observables, and Berry-phase properties.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use,intrinsic :: ISO_C_binding
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
   use scitools_utils,only: str,stop_error
   use wan_hamiltonian,only: wann90_tb_t
   use spfft
   implicit none
   include '../units_inc.f90'
   include '../formats.h'
   include 'fftw3.f03'
!--------------------------------------------------------------------------------------
   private
   public :: wann_spfft_t
!--------------------------------------------------------------------------------------
   type :: wann_spfft_t
      integer                                    :: nwan
      integer                                    :: nx,ny,nz,nrpts
      integer                                    :: nkx,nky,nkz,nkpts,nkpts_loc,kdim
      integer,allocatable,dimension(:)           :: ndegen
      real(dp),allocatable,dimension(:,:)        :: crvec
      complex(dp),allocatable,dimension(:)       :: HA_r
      complex(dp),allocatable,dimension(:,:)     :: gradH_R
      complex(dp),allocatable,dimension(:,:,:)   :: ham_r
      complex(dp),allocatable,dimension(:,:,:,:) :: pos_r
      ! .. SpFFT .. 
      integer :: processingUnit = 1
      integer :: maxNumThreads = -1
      integer :: maxNumLocalZColumns
      type(c_ptr) :: grid = c_null_ptr
      type(c_ptr) :: transform = c_null_ptr
      type(c_ptr) :: realValuesPtr
      integer,allocatable,dimension(:)           :: indices
      complex(C_DOUBLE_COMPLEX),allocatable,dimension(:) :: frequencyElements
      real(C_DOUBLE),allocatable,dimension(:)    :: spaceDomain
   contains
      procedure, public   :: InitFromW90
      procedure, public   :: Clean
      procedure, public   :: GetHam
      procedure, public   :: GetGradHam
      procedure, public   :: GetDipole
   end type wann_spfft_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine InitFromW90(me,w90,nk)
      class(wann_spfft_t),target :: me
      type(wann90_tb_t),intent(in) :: w90
      integer,intent(in) :: nk(:)
      integer :: ir,irw,ix,iy,iz,i,j,counter
      integer :: kx,ky
      integer :: ix_bound(2),iy_bound(2),iz_bound(2)
      integer :: ierr,tid
      integer :: errorCode = 0
      integer,allocatable :: w90_rindx(:)
      integer,allocatable :: indx_2d(:,:), indx_3d(:,:,:)

      me%nwan = w90%num_wann
      me%kdim = size(nk)

      ix_bound(1) = minval(w90%irvec(:,1))
      ix_bound(2) = maxval(w90%irvec(:,1))
      iy_bound(1) = minval(w90%irvec(:,2))
      iy_bound(2) = maxval(w90%irvec(:,2))
      iz_bound(1) = minval(w90%irvec(:,3))
      iz_bound(2) = maxval(w90%irvec(:,3))     

      me%nkx = 1
      me%nky = 1
      me%nkz = 1

      select case(me%kdim)
      case(1)
         me%nx = ix_bound(2) - ix_bound(1) + 1
         if(mod(me%nx,2) /= 0) me%nx = me%nx + 1
         me%ny = 1
         me%nz = 1
         me%nkx = nk(1)

         if(me%nkx < me%nx) then
            write(error_unit,fmt900) "Less k-points than R-points."
            write(output_unit,'(A)') " ---> real-space grid: "//str(me%nx)
            write(output_unit,'(A)') " ---> k-space grid: "//str(me%nkx)
         end if

         me%maxNumLocalZColumns = me%nkx
         allocate(me%indices(me%nx * 3)); me%indices = 0
         allocate(me%spaceDomain(2 * me%nkx))

         counter = 0
         do ix=1,me%nx
            me%indices(counter * 3 + 1) = ix - 1
            counter = counter + 1
         end do

      case(2)
         me%nx = ix_bound(2) - ix_bound(1) + 1
         me%ny = iy_bound(2) - iy_bound(1) + 1
         if(mod(me%nx,2) /= 0) me%nx = me%nx + 1
         if(mod(me%ny,2) /= 0) me%ny = me%ny + 1
         me%nz = 1

         me%nkx = nk(1)
         me%nky = nk(2)

         if(me%nkx < me%nx .or. me%nky < me%ny) then
            write(error_unit,fmt900) "Less k-points than R-points."
            write(output_unit,'(A)') " ---> real-space grid: "//str(me%nx)//"x"//str(me%ny)
            write(output_unit,'(A)') " ---> k-space grid: "//str(me%nkx)//"x"//str(me%nky)
         end if

         me%maxNumLocalZColumns = me%nkx * me%nky
         allocate(me%indices(me%nx * me%ny * 3)); me%indices = 0
         allocate(me%spaceDomain(2 * me%nkx * me%nky))

         counter = 0
         do iy=1,me%ny
            do ix=1,me%nx
               me%indices(counter * 3 + 1) = ix - 1
               me%indices(counter * 3 + 2) = iy - 1
               counter = counter + 1
            end do
         end do        

      case(3)
         me%nx = ix_bound(2) - ix_bound(1) + 1
         me%ny = iy_bound(2) - iy_bound(1) + 1
         me%nz = iz_bound(2) - iz_bound(1) + 1
         if(mod(me%nx,2) /= 0) me%nx = me%nx + 1
         if(mod(me%ny,2) /= 0) me%ny = me%ny + 1
         if(mod(me%nz,2) /= 0) me%nz = me%nz + 1

         me%nkx = nk(1)
         me%nky = nk(2)
         me%nkz = nk(3)

         if(me%nkx < me%nx .or. me%nky < me%ny .or. me%nkz < me%nz) then
            write(error_unit,fmt900) "Less k-points than R-points."
            write(output_unit,'(A)') " ---> real-space grid: "//str(me%nx)//"x"//str(me%ny)//"x"//str(me%nz)
            write(output_unit,'(A)') " ---> k-space grid: "//str(me%nkx)//"x"//str(me%nky)//"x"//str(me%nkz)
         end if

         me%maxNumLocalZColumns = me%nx * me%ny
         allocate(me%indices(me%nx * me%ny * me%nz * 3)); me%indices = 0
         allocate(me%spaceDomain(2 * me%nkx * me%nky * me%nkz))

         counter = 0
         do iz=1,me%nz
            do iy=1,me%ny
               do ix=1,me%nx
                  me%indices(counter * 3 + 1) = ix - 1
                  me%indices(counter * 3 + 2) = iy - 1
                  me%indices(counter * 3 + 3) = iz - 1
                  counter = counter + 1
               end do
            end do    
         end do    

      end select

      me%nrpts = me%nx * me%ny * me%nz
      me%nkpts = me%nkx * me%nky * me%nkz

      allocate(me%HA_r(me%nrpts))
      allocate(me%gradH_R(me%nrpts,3))
      allocate(me%frequencyElements(me%nkpts))

      allocate(w90_rindx(me%nrpts)); w90_rindx = -999

      select case(me%kdim)
      case(1)
         do ir=1,w90%nrpts
            ix = GetFFTIndex(me%nx, w90%irvec(ir,1))
            w90_rindx(ix) = ir
         end do         
      case(2)
         allocate(indx_2d(me%nx,me%ny)); indx_2d = -999
         do ir=1,w90%nrpts
            ix = GetFFTIndex(me%nx, w90%irvec(ir,1))
            iy = GetFFTIndex(me%ny, w90%irvec(ir,2))
            indx_2d(ix,iy) = ir
         end do         
         w90_rindx = reshape(indx_2d, [me%nrpts])  
         deallocate(indx_2d)
      case(3)
         allocate(indx_3d(me%nx,me%ny,me%nz)); indx_3d = -999
         do ir=1,w90%nrpts
            ix = GetFFTIndex(me%nx, w90%irvec(ir,1))
            iy = GetFFTIndex(me%ny, w90%irvec(ir,2))
            iz = GetFFTIndex(me%nz, w90%irvec(ir,3))
            indx_3d(ix,iy,iz) = ir
         end do  
         w90_rindx = reshape(indx_3d, [me%nrpts])   
         deallocate(indx_3d)              
      end select

      allocate(me%ndegen(me%nrpts)); me%ndegen = 1000000
      allocate(me%ham_r(me%nrpts,me%nwan,me%nwan)); me%ham_r = zero
      allocate(me%pos_r(me%nrpts,me%nwan,me%nwan,3)); me%pos_r = zero
      allocate(me%crvec(me%nrpts,3)); me%crvec = 0.0_dp

      do ir=1,me%nrpts
         irw = w90_rindx(ir)
         if(irw == -999) cycle
         me%ndegen(ir) = w90%ndegen(irw)
         me%crvec(ir,:) = w90%crvec(irw,:)
         do j=1,me%nwan
            do i=1,me%nwan
               me%ham_r(ir,i,j) = w90%ham_r(i,j,irw) / me%ndegen(ir)
               me%pos_r(ir,i,j,:) = w90%pos_r(i,j,:,irw) / me%ndegen(ir)
            end do
         end do
      end do

      deallocate(w90_rindx)

      ! create grid
      errorCode = spfft_grid_create(me%grid, me%nkx, me%nky, me%nkz, me%maxNumLocalZColumns, &
         me%processingUnit, me%maxNumThreads)
      if (errorCode /= SPFFT_SUCCESS) error stop

      ! create transform
      ! Note: A transform handle can be created without a grid if no resource sharing is desired.
      errorCode = spfft_transform_create(me%transform, me%grid, me%processingUnit, 0, me%nkx, me%nky, me%nkz, me%nkz,&
        me%nrpts, SPFFT_INDEX_TRIPLETS, me%indices)
      if (errorCode /= SPFFT_SUCCESS) error stop

      ! grid can be safely destroyed after creating all required transforms
      errorCode = spfft_grid_destroy(me%grid)
      if (errorCode /= SPFFT_SUCCESS) error stop

   end subroutine InitFromW90
!--------------------------------------------------------------------------------------
   subroutine Clean(me)
      class(wann_spfft_t) :: me

      if(allocated(me%ham_r)) deallocate(me%ham_r)
      if(allocated(me%pos_r)) deallocate(me%pos_r)
      if(allocated(me%ndegen)) deallocate(me%ndegen)
      if(allocated(me%crvec)) deallocate(me%crvec)
      if(allocated(me%HA_r)) deallocate(me%HA_r)
      if(allocated(me%gradH_R)) deallocate(me%gradH_R)
      if(allocated(me%indices)) deallocate(me%indices)
      if(allocated(me%frequencyElements)) deallocate(me%frequencyElements)
      if(allocated(me%spaceDomain)) deallocate(me%spaceDomain)

   end subroutine Clean
!--------------------------------------------------------------------------------------
   subroutine GetHam(me,Hk,Ar)
      class(wann_spfft_t) :: me
      complex(dp),intent(inout) :: Hk(:,:,:)
      real(dp),intent(in),optional :: Ar(3)
      integer :: i,j,ik,ir
      integer,allocatable :: rdims(:), kdims(:)
      integer :: errorCode = 0
      complex(C_DOUBLE_COMPLEX), pointer :: spaceDomainPtr(:)

      call assert_shape(Hk, [me%nwan,me%nwan,me%nkpts_loc], "GetHam", "Hk")

      select case(me%kdim)
      case(1)
         allocate(rdims(1),kdims(1))
         rdims(1) = me%nx
         kdims(1) = me%nkx 
      case(2)      
         allocate(rdims(2),kdims(2))
         rdims(:) = [me%nx, me%ny]
         kdims(:) = [me%nkx, me%nky]
      case(3)
         allocate(rdims(3),kdims(3))
         rdims(:) = [me%nx, me%ny, me%nz]
         kdims(:) = [me%nkx, me%nky, me%nkz]
      end select

      do j=1,me%nwan
         do i=1,j
            me%HA_r = me%ham_r(:,i,j)
            if(present(Ar)) then
               call ApplyPhaseFactor(rdims, Ar, me%HA_r)
            end if

            ! set space domain array to use memory allocted by the library
            errorCode = spfft_transform_get_space_domain(me%transform, me%processingUnit, me%realValuesPtr)
            if (errorCode /= SPFFT_SUCCESS) error stop

            ! transform backward
            errorCode = spfft_transform_backward(me%transform, me%HA_r, me%processingUnit)
            if (errorCode /= SPFFT_SUCCESS) error stop

            call c_f_pointer(me%realValuesPtr, spaceDomainPtr, [me%nkpts])

            do ik=1,me%nkpts
               Hk(i,j,ik) = spaceDomainPtr(ik)
               if(i < j) Hk(j,i,ik) = conjg(Hk(i,j,ik))
            end do
         end do
      end do

      deallocate(rdims,kdims)

   end subroutine GetHam
!--------------------------------------------------------------------------------------
   subroutine GetGradHam(me,GradHk,Ar)
      class(wann_spfft_t) :: me
      complex(dp),intent(inout) :: GradHk(:,:,:,:)
      real(dp),intent(in),optional :: Ar(3)
      integer :: i,j,idir,ik,ik_glob
      integer,allocatable :: rdims(:), kdims(:)
      integer :: errorCode = 0
      complex(C_DOUBLE_COMPLEX), pointer :: spaceDomainPtr(:)

      call assert_shape(GradHk, [me%nwan,me%nwan,3,me%nkpts_loc], "GetGradHam", "GradHk")

      select case(me%kdim)
      case(1)
         allocate(rdims(1),kdims(1))
         rdims(1) = me%nx
         kdims(1) = me%nkx 
      case(2)
    
         allocate(rdims(2),kdims(2))
         rdims(:) = [me%nx, me%ny]
         kdims(:) = [me%nkx, me%nky]
      case(3)
         allocate(rdims(3),kdims(3))
         rdims(:) = [me%nx, me%ny, me%nz]
         kdims(:) = [me%nkx, me%nky, me%nkz]
      end select

      do j=1,me%nwan
         do i=1,j
            call GetGradiant(me%crvec, me%ham_r(:,i,j), me%gradH_R)
            do idir=1,3
               if(present(Ar)) then
                  call ApplyPhaseFactor(rdims, Ar, me%gradH_R(:,idir))
               end if

               ! set space domain array to use memory allocted by the library
               errorCode = spfft_transform_get_space_domain(me%transform, me%processingUnit, me%realValuesPtr)
               if (errorCode /= SPFFT_SUCCESS) error stop

               ! transform backward
               errorCode = spfft_transform_backward(me%transform, me%gradH_R(:,idir), me%processingUnit)
               if (errorCode /= SPFFT_SUCCESS) error stop

               call c_f_pointer(me%realValuesPtr, spaceDomainPtr, [me%nkpts])

               do ik=1,me%nkpts
                  GradHk(i,j,idir,ik) = spaceDomainPtr(ik)
                  if(i < j) GradHk(j,i,idir,ik) = conjg(GradHk(i,j,idir,ik))
               end do
            end do
         end do
      end do

      deallocate(rdims,kdims)

   end subroutine GetGradHam
!--------------------------------------------------------------------------------------
     subroutine GetDipole(me,Dk,Ar)
      class(wann_spfft_t) :: me
      complex(dp),intent(inout) :: Dk(:,:,:,:)
      real(dp),intent(in),optional :: Ar(3)
      integer :: i,j,idir,ik,ik_glob
      integer,allocatable :: rdims(:), kdims(:)
      integer :: errorCode = 0
      complex(C_DOUBLE_COMPLEX), pointer :: spaceDomainPtr(:)

      call assert_shape(Dk, [me%nwan,me%nwan,3,me%nkpts_loc], "GetDipole", "Dk")

      select case(me%kdim)
      case(1)
         allocate(rdims(1),kdims(1))
         rdims(1) = me%nx
         kdims(1) = me%nkx 
      case(2)
    
         allocate(rdims(2),kdims(2))
         rdims(:) = [me%nx, me%ny]
         kdims(:) = [me%nkx, me%nky]
      case(3)
         allocate(rdims(3),kdims(3))
         rdims(:) = [me%nx, me%ny, me%nz]
         kdims(:) = [me%nkx, me%nky, me%nkz]
      end select


      do j=1,me%nwan
         do i=1,j
            do idir=1,3
               me%HA_r = me%pos_r(:,i,j,idir)
               if(present(Ar)) then
                  call ApplyPhaseFactor(rdims, Ar, me%HA_r)
               end if

               ! set space domain array to use memory allocted by the library
               errorCode = spfft_transform_get_space_domain(me%transform, me%processingUnit, me%realValuesPtr)
               if (errorCode /= SPFFT_SUCCESS) error stop

               ! transform backward
               errorCode = spfft_transform_backward(me%transform, me%HA_r, me%processingUnit)
               if (errorCode /= SPFFT_SUCCESS) error stop

               call c_f_pointer(me%realValuesPtr, spaceDomainPtr, [me%nkpts])

               do ik=1,me%nkpts_loc
                  Dk(i,j,idir,ik) = spaceDomainPtr(ik)
                  if(i < j) Dk(j,i,idir,ik) = conjg(Dk(i,j,idir,ik))
               end do
            end do
         end do
      end do

      deallocate(rdims,kdims)

   end subroutine GetDipole
!--------------------------------------------------------------------------------------
   pure elemental function GetFFTIndex(Nx,ilat) result(ix)
      integer,intent(in) :: Nx
      integer,intent(in) :: ilat
      integer :: ix

      if(ilat >= 0) then
         ix = ilat + 1 
      else
         ix = Nx + ilat + 1
      end if

   end function GetFFTIndex   
!--------------------------------------------------------------------------------------
   pure elemental function FFT_Freq(n,i) result(k)
      integer,intent(in) :: n, i
      integer :: k

      if(i <= n/2) then
         k = i - 1
      else
         k = -(n - i + 1)
      end if

   end function FFT_Freq
!--------------------------------------------------------------------------------------
   subroutine ApplyPhaseFactor(nk,Ar,OO_R)
      integer,intent(in)  :: nk(:)
      real(dp),intent(in) :: Ar(3)
      complex(dp),target,intent(inout) :: OO_R(:)
      integer :: kdim,nx,ny,nz,ix,iy,iz,kx,ky,kz
      real(dp) :: adot,c,s
      complex(dp),pointer :: OO_1d(:)
      complex(dp),pointer :: OO_2d(:,:)
      complex(dp),pointer :: OO_3d(:,:,:)

      kdim = size(nk)

      select case(kdim)
      case(1)
         nx = nk(1)
         OO_1d(1:nx) => OO_R
         do ix=1,nx
            kx = FFT_Freq(nx, ix)
            adot = Ar(1) * kx
            c = cos(DPI * adot)
            s = sin(DPI * adot)
            OO_1d(ix) = cmplx(c,-s,kind=dp) * OO_1d(ix)
         end do
      case(2)
         nx = nk(1); ny = nk(2)
         OO_2d(1:nx, 1:ny) => OO_R
         do concurrent(iy=1:ny, ix=1:nx)
            kx = FFT_Freq(nx, ix)
            ky = FFT_Freq(ny, iy)
            adot = Ar(1) * kx + Ar(2) * ky
            c = cos(DPI * adot)
            s = sin(DPI * adot)
            OO_2d(ix,iy) = cmplx(c,-s,kind=dp) * OO_2d(ix,iy)
         end do         
      case(3)
         nx = nk(1); ny = nk(2); nz = nk(3)
         OO_3d(1:nx, 1:ny, 1:nz) => OO_R
         do concurrent(iz=1:nz, iy=1:ny, ix=1:nx)
            kx = FFT_Freq(nx, ix)
            ky = FFT_Freq(ny, iy)
            kz = FFT_Freq(nz, iz)
            adot = Ar(1) * kx + Ar(2) * ky  + Ar(3) * kz
            c = cos(DPI * adot)
            s = sin(DPI * adot)
            OO_3d(ix,iy,iz) = cmplx(c,-s,kind=dp) * OO_3d(ix,iy,iz)
         end do   
      end select

   end subroutine ApplyPhaseFactor
!--------------------------------------------------------------------------------------
   subroutine GetGradiant(crvec,OO_R,vec_R)
      real(dp),intent(in) :: crvec(:,:)
      complex(dp),intent(in) :: OO_R(:)      
      complex(dp),intent(inout) :: vec_R(:,:)
      integer :: nrpts,ir

      nrpts = size(crvec, dim=1)
      vec_R(1:nrpts,1) = iu * crvec(1:nrpts,1) * OO_R(1:nrpts)
      vec_R(1:nrpts,2) = iu * crvec(1:nrpts,2) * OO_R(1:nrpts)
      vec_R(1:nrpts,3) = iu * crvec(1:nrpts,3) * OO_R(1:nrpts)
   
   end subroutine GetGradiant
!--------------------------------------------------------------------------------------


!======================================================================================

end module wan_spfft_ham
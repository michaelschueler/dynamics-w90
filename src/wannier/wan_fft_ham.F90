module wan_fft_ham
!! Provides tools for reading/writing the Wannier Hamiltonian and for computing 
!! band structures, observables, and Berry-phase properties.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
   use scitools_utils,only: str
   use wan_hamiltonian,only: wann90_tb_t
   implicit none
   include '../units_inc.f90'
   include '../formats.h'
   include '../fftw3.fh'
!--------------------------------------------------------------------------------------
   private
   public :: wann_fft_t
!--------------------------------------------------------------------------------------
   
!--------------------------------------------------------------------------------------
   type :: wann_fft_t
      integer                                    :: nwan
      integer                                    :: nx,ny,nz,nrpts
      integer                                    :: nkx,nky,nkz,nkpts,kdim
      integer,allocatable,dimension(:)           :: ndegen
      real(dp),allocatable,dimension(:,:)        :: crvec
      complex(dp),allocatable,dimension(:,:,:)   :: ham_r
      complex(dp),allocatable,dimension(:,:,:,:) :: pos_r
      ! .. FFTW ..
      integer(kind=8) :: plan_fw,plan_bw
   contains
      procedure, public   :: InitFromW90
      procedure, public   :: Clean
      procedure, public   :: GetHam
      procedure, public   :: GetGradHam
      procedure, public   :: GetHam_Dressed
      procedure, public   :: GetGradHam_Dressed
      procedure, private  :: GetHam_1d
      procedure, private  :: GetHam_2d
      procedure, private  :: GetHam_3d
      procedure, private  :: GetGradHam_1d
      procedure, private  :: GetGradHam_2d
      procedure, private  :: GetGradHam_3d
      procedure, private  :: GetHam_Dressed_1d
      procedure, private  :: GetHam_Dressed_2d
      procedure, private  :: GetHam_Dressed_3d
      procedure, private  :: GetGradHam_Dressed_1d
      procedure, private  :: GetGradHam_Dressed_2d
      procedure, private  :: GetGradHam_Dressed_3d
      procedure, public   :: GetDipole
      procedure, public   :: GetDipole_Dressed
      procedure, private  :: GetDipole_1d
      procedure, private  :: GetDipole_2d
      procedure, private  :: GetDipole_3d
      procedure, private  :: GetDipole_Dressed_1d
      procedure, private  :: GetDipole_Dressed_2d
      procedure, private  :: GetDipole_Dressed_3d
   end type wann_fft_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine InitFromW90(me,w90,nk)
      class(wann_fft_t) :: me
      type(wann90_tb_t),intent(in) :: w90
      integer,intent(in) :: nk(:)
      integer :: ir,irw,ix,iy,iz,i,j
      integer :: kx,ky
      integer :: ix_bound(2),iy_bound(2),iz_bound(2)
      integer,allocatable :: w90_rindx(:)
      integer,allocatable :: indx_2d(:,:), indx_3d(:,:,:)
      complex(dp),allocatable :: zr_1d(:),zk_1d(:)
      complex(dp),allocatable :: zr_2d(:,:),zk_2d(:,:)
      complex(dp),allocatable :: zr_3d(:,:,:),zk_3d(:,:,:)
  
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
            stop
         end if

         allocate(zr_1d(me%nkx),zk_1d(me%nkx))
         call DFFTW_PLAN_DFT_1D(me%plan_fw,me%nkx,zk_1d,zr_1d,FFTW_FORWARD,FFTW_MEASURE)
         call DFFTW_PLAN_DFT_1D(me%plan_bw,me%nkx,zr_1d,zk_1d,FFTW_BACKWARD,FFTW_MEASURE)
         deallocate(zr_1d,zk_1d)
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
            stop
         end if

         allocate(zr_2d(me%nkx,me%nky),zk_2d(me%nkx,me%nky))
         call DFFTW_PLAN_DFT_2D(me%plan_fw,me%nkx,me%nky,zk_2d,zr_2d,FFTW_FORWARD,FFTW_MEASURE)
         call DFFTW_PLAN_DFT_2D(me%plan_bw,me%nkx,me%nky,zr_2d,zk_2d,FFTW_BACKWARD,FFTW_MEASURE)
         deallocate(zr_2d,zk_2d)
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
            stop
         end if

         allocate(zr_3d(me%nkx,me%nky,me%nkz),zk_3d(me%nkx,me%nky,me%nkz))
         call DFFTW_PLAN_DFT_3D(me%plan_fw,me%nkx,me%nky,me%nkz,zk_3d,zr_3d,FFTW_FORWARD,FFTW_MEASURE)
         call DFFTW_PLAN_DFT_3D(me%plan_bw,me%nkx,me%nky,me%nkz,zr_3d,zk_3d,FFTW_BACKWARD,FFTW_MEASURE)
         deallocate(zr_3d,zk_3d)

      end select

      me%nrpts = me%nx * me%ny * me%nz
      me%nkpts = me%nkx * me%nky * me%nkz

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

   end subroutine InitFromW90
!--------------------------------------------------------------------------------------
   subroutine Clean(me)
      class(wann_fft_t) :: me

      if(allocated(me%ham_r)) deallocate(me%ham_r)
      if(allocated(me%pos_r)) deallocate(me%pos_r)
      if(allocated(me%ndegen)) deallocate(me%ndegen)
      if(allocated(me%crvec)) deallocate(me%crvec)

   end subroutine Clean
!--------------------------------------------------------------------------------------
   subroutine GetHam(me,Hk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: Hk(:,:,:)

      call assert_shape(Hk, [me%nwan,me%nwan,me%nkpts], "GetHam", "Hk")

      select case(me%kdim)
      case(1)
         call me%GetHam_1d(Hk)
      case(2)
         call me%GetHam_2d(Hk)
      case(3)
         call me%GetHam_3d(Hk)
      end select

   end subroutine GetHam
!--------------------------------------------------------------------------------------
   subroutine GetGradHam(me,GradHk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: GradHk(:,:,:,:)

      call assert_shape(GradHk, [me%nwan,me%nwan,3,me%nkpts], "GetGradHam", "GradHk")

      select case(me%kdim)
      case(1)
         call me%GetGradHam_1d(GradHk)
      case(2)
         call me%GetGradHam_2d(GradHk)
      case(3)
         call me%GetGradHam_3d(GradHk)
      end select

   end subroutine GetGradHam
!--------------------------------------------------------------------------------------
   subroutine GetHam_Dressed(me,Ar,Hk)
      class(wann_fft_t) :: me
      real(dp),intent(in) :: Ar(3)
      complex(dp),intent(inout) :: Hk(:,:,:)

      call assert_shape(Hk, [me%nwan,me%nwan,me%nkpts], "GetHam", "Hk")

      select case(me%kdim)
      case(1)
         call me%GetHam_Dressed_1d(Ar,Hk)
      case(2)
         call me%GetHam_Dressed_2d(Ar,Hk)
      case(3)
         call me%GetHam_Dressed_3d(Ar,Hk)
      end select

   end subroutine GetHam_Dressed
!--------------------------------------------------------------------------------------
   subroutine GetGradHam_Dressed(me,Ar,GradHk)
      class(wann_fft_t) :: me
      real(dp),intent(in) :: Ar(3)
      complex(dp),intent(inout) :: GradHk(:,:,:,:)

      call assert_shape(GradHk, [me%nwan,me%nwan,3,me%nkpts], "GetGradHam_Dressed", "GradHk")

      select case(me%kdim)
      case(1)
         call me%GetGradHam_Dressed_1d(Ar,GradHk)
      case(2)
         call me%GetGradHam_Dressed_2d(Ar,GradHk)
      case(3)
         call me%GetGradHam_Dressed_3d(Ar,GradHk)
      end select

   end subroutine GetGradHam_Dressed
!--------------------------------------------------------------------------------------
   subroutine GetDipole(me,Dk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: Dk(:,:,:,:)

      call assert_shape(Dk, [me%nwan,me%nwan,3,me%nkpts], "GetDipole", "Dk")

      select case(me%kdim)
      case(1)
         call me%GetDipole_1d(Dk)
      case(2)
         call me%GetDipole_2d(Dk)
      case(3)
         call me%GetDipole_3d(Dk)
      end select

   end subroutine GetDipole
!--------------------------------------------------------------------------------------
   subroutine GetDipole_Dressed(me,Ar,Dk)
      class(wann_fft_t) :: me
      real(dp),intent(in) :: Ar(3)
      complex(dp),intent(inout) :: Dk(:,:,:,:)

      call assert_shape(Dk, [me%nwan,me%nwan,3,me%nkpts], "GetDipole", "Dk")

      select case(me%kdim)
      case(1)
         call me%GetDipole_Dressed_1d(Ar,Dk)
      case(2)
         call me%GetDipole_Dressed_2d(Ar,Dk)
      case(3)
         call me%GetDipole_Dressed_3d(Ar,Dk)
      end select

   end subroutine GetDipole_Dressed
!--------------------------------------------------------------------------------------
   subroutine GetHam_1d(me,Hk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: Hk(:,:,:)
      integer :: i,j,ik
      complex(dp),allocatable :: work_r(:),work_k(:)

      !$OMP PARALLEL PRIVATE(i,j,ik,work_r,work_k)
      allocate(work_r(me%nkx),work_k(me%nkx))
      do j=1,me%nwan
         do i=1,me%nwan
            call Smooth2Dense_1d(me%nx, me%nkx, me%ham_r(:,i,j), work_r)
            call dfftw_execute_dft(me%plan_bw,work_r,work_k)
            do ik=1,me%nkx
               Hk(i,j,ik) = work_k(ik)
            end do
         end do
      end do
      deallocate(work_r,work_k)
      !$OMP END PARALLEL

   end subroutine GetHam_1d
!--------------------------------------------------------------------------------------
   subroutine GetHam_2d(me,Hk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: Hk(:,:,:)
      integer :: i,j,ik
      complex(dp),allocatable :: work_r(:,:),work_k(:,:),work_1d(:)

      !$OMP PARALLEL PRIVATE(i,j,ik,work_r,work_k,work_1d)
      allocate(work_r(me%nkx,me%nky),work_k(me%nkx,me%nky),work_1d(me%nkpts))
      !$OMP DO COLLAPSE(2)
      do j=1,me%nwan
         do i=1,me%nwan
            call Smooth2Dense_2d(me%nx, me%ny, me%nkx, me%nky, me%ham_r(:,i,j), work_r)
            call dfftw_execute_dft(me%plan_bw,work_r,work_k)
            work_1d = reshape(work_k, [me%nkpts])
            do ik=1,me%nkpts
               Hk(i,j,ik) = work_1d(ik)
            end do
         end do
      end do
      !$OMP END DO
      deallocate(work_r,work_k,work_1d)
      !$OMP END PARALLEL

   end subroutine GetHam_2d
!--------------------------------------------------------------------------------------
   subroutine GetHam_3d(me,Hk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: Hk(:,:,:)
      integer :: i,j,ik
      complex(dp),allocatable :: work_r(:,:,:),work_k(:,:,:),work_1d(:)

      !$OMP PARALLEL PRIVATE(i,j,work_r,work_k,work_1d)
      allocate(work_r(me%nkx,me%nky,me%nkz),work_k(me%nkx,me%nky,me%nkz),work_1d(me%nkpts))
      !$OMP DO COLLAPSE(2)
      do j=1,me%nwan
         do i=1,me%nwan
            call Smooth2Dense_3d(me%nx, me%ny, me%nz, me%nkx, me%nky, me%nkz, me%ham_r(:,i,j), work_r)
            call dfftw_execute_dft(me%plan_bw,work_r,work_k)
            work_1d = reshape(work_k, [me%nkpts])
            do ik=1,me%nkpts
               Hk(i,j,ik) = work_1d(ik)
            end do
         end do
      end do
      !$OMP END DO
      deallocate(work_r,work_k,work_1d)
      !$OMP END PARALLEL

   end subroutine GetHam_3d
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
   subroutine GetGradHam_1d(me,GradHk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: GradHk(:,:,:,:)
      integer :: i,j,ik,idir
      complex(dp),allocatable :: work_r(:),work_k(:)
      complex(dp),allocatable :: gradH_R(:,:)

      !$OMP PARALLEL PRIVATE(i,j,ik,idir,work_r,work_k,gradH_R)
      allocate(work_r(me%nkx),work_k(me%nkx))
      allocate(gradH_R(me%nrpts,3))
      do j=1,me%nwan
         do i=1,me%nwan
            ! call me%DressCrvec(me%ham_r(:,i,j), gradH_R)
            call GetGradiant(me%crvec, me%ham_r(:,i,j), gradH_R)
            do idir=1,3
               call Smooth2Dense_1d(me%nx, me%nkx, gradH_R(:,idir), work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               do ik=1,me%nkx
                  GradHk(i,j,idir,ik) = work_k(ik)
               end do
            end do
         end do
      end do
      deallocate(work_r,work_k)
      deallocate(gradH_R)
      !$OMP END PARALLEL

   end subroutine GetGradHam_1d
!--------------------------------------------------------------------------------------
   subroutine GetGradHam_2d(me,GradHk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: GradHk(:,:,:,:)
      integer :: i,j,ik,idir
      complex(dp),allocatable :: work_r(:,:),work_k(:,:),work_1d(:)
      complex(dp),allocatable :: gradH_R(:,:)

      !$OMP PARALLEL PRIVATE(i,j,ik,work_r,work_k,work_1d,gradH_R)
      allocate(work_r(me%nkx,me%nky),work_k(me%nkx,me%nky),work_1d(me%nkpts))
      allocate(gradH_R(me%nrpts,3))
      !$OMP DO COLLAPSE(2)
      do j=1,me%nwan
         do i=1,me%nwan
            ! call me%DressCrvec(me%ham_r(:,i,j), gradH_R)
            call GetGradiant(me%crvec, me%ham_r(:,i,j), gradH_R)
            do idir=1,3
               call Smooth2Dense_2d(me%nx, me%ny, me%nkx, me%nky, gradH_R(:,idir), work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               work_1d = reshape(work_k, [me%nkpts])
               do ik=1,me%nkpts
                  GradHk(i,j,idir,ik) = work_1d(ik)
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(work_r,work_k,work_1d)
      deallocate(gradH_R)
      !$OMP END PARALLEL

   end subroutine GetGradHam_2d
!--------------------------------------------------------------------------------------
   subroutine GetGradHam_3d(me,GradHk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: GradHk(:,:,:,:)
      integer :: i,j,ik,idir
      complex(dp),allocatable :: work_r(:,:,:),work_k(:,:,:),work_1d(:)
      complex(dp),allocatable :: gradH_R(:,:)

      !$OMP PARALLEL PRIVATE(i,j,work_r,work_k,work_1d,gradH_R)
      allocate(work_r(me%nkx,me%nky,me%nkz),work_k(me%nkx,me%nky,me%nkz),work_1d(me%nkpts))
      allocate(gradH_R(me%nrpts,3))
      !$OMP DO COLLAPSE(2)
      do j=1,me%nwan
         do i=1,me%nwan
            ! call me%DressCrvec(me%ham_r(:,i,j), gradH_R)
            call GetGradiant(me%crvec, me%ham_r(:,i,j), gradH_R)
            do idir=1,3
               call Smooth2Dense_3d(me%nx, me%ny, me%nz, me%nkx, me%nky, me%nkz, gradH_R(:,idir), work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               work_1d = reshape(work_k, [me%nkpts])
               do ik=1,me%nkpts
                  GradHk(i,j,idir,ik) = work_1d(ik)
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(work_r,work_k,work_1d)
      deallocate(gradH_R)
      !$OMP END PARALLEL

   end subroutine GetGradHam_3d
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
   subroutine GetGradHam_Dressed_1d(me,Ar,GradHk)
      class(wann_fft_t) :: me
      real(dp),intent(in)       :: Ar(3)
      complex(dp),intent(inout) :: GradHk(:,:,:,:)
      integer :: i,j,ik,idir
      complex(dp),allocatable :: work_r(:),work_k(:)
      complex(dp),allocatable :: gradH_R(:,:)

      !$OMP PARALLEL PRIVATE(i,j,ik,idir,work_r,work_k,gradH_R)
      allocate(work_r(me%nkx),work_k(me%nkx))
      allocate(gradH_R(me%nrpts,3))
      do j=1,me%nwan
         do i=1,me%nwan
            ! call me%DressCrvec(me%ham_r(:,i,j), gradH_R)
            call GetGradiant(me%crvec, me%ham_r(:,i,j), gradH_R)
            do idir=1,3
               ! call me%DressPhase(Ar, gradH_R(:,idir))
               call ApplyPhaseFactor([me%nx], Ar, gradH_R(:,idir))
               call Smooth2Dense_1d(me%nx, me%nkx, gradH_R(:,idir), work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               do ik=1,me%nkx
                  GradHk(i,j,idir,ik) = work_k(ik)
               end do
            end do
         end do
      end do
      deallocate(work_r,work_k)
      deallocate(gradH_R)
      !$OMP END PARALLEL

   end subroutine GetGradHam_Dressed_1d
!--------------------------------------------------------------------------------------
   subroutine GetGradHam_Dressed_2d(me,Ar,GradHk)
      class(wann_fft_t) :: me
      real(dp),intent(in)       :: Ar(3)
      complex(dp),intent(inout) :: GradHk(:,:,:,:)
      integer :: i,j,ik,idir
      complex(dp),allocatable :: work_r(:,:),work_k(:,:),work_1d(:)
      complex(dp),allocatable :: gradH_R(:,:)

      !$OMP PARALLEL PRIVATE(i,j,ik,work_r,work_k,work_1d,gradH_R)
      allocate(work_r(me%nkx,me%nky),work_k(me%nkx,me%nky),work_1d(me%nkpts))
      allocate(gradH_R(me%nrpts,3))
      !$OMP DO COLLAPSE(2)
      do j=1,me%nwan
         do i=1,me%nwan
            ! call me%DressCrvec(me%ham_r(:,i,j), gradH_R)
            call GetGradiant(me%crvec, me%ham_r(:,i,j), gradH_R)
            do idir=1,3
               ! call me%DressPhase(Ar, gradH_R(:,idir))
               call ApplyPhaseFactor([me%nx,me%ny], Ar, gradH_R(:,idir))
               call Smooth2Dense_2d(me%nx, me%ny, me%nkx, me%nky, gradH_R(:,idir), work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               work_1d = reshape(work_k, [me%nkpts])
               do ik=1,me%nkpts
                  GradHk(i,j,idir,ik) = work_1d(ik)
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(work_r,work_k,work_1d)
      deallocate(gradH_R)
      !$OMP END PARALLEL

   end subroutine GetGradHam_Dressed_2d
!--------------------------------------------------------------------------------------
   subroutine GetGradHam_Dressed_3d(me,Ar,GradHk)
      class(wann_fft_t) :: me
      real(dp),intent(in)       :: Ar(3)
      complex(dp),intent(inout) :: GradHk(:,:,:,:)
      integer :: i,j,ik,idir
      complex(dp),allocatable :: work_r(:,:,:),work_k(:,:,:),work_1d(:)
      complex(dp),allocatable :: gradH_R(:,:)

      !$OMP PARALLEL PRIVATE(i,j,work_r,work_k,work_1d,gradH_R)
      allocate(work_r(me%nkx,me%nky,me%nkz),work_k(me%nkx,me%nky,me%nkz),work_1d(me%nkpts))
      allocate(gradH_R(me%nrpts,3))
      !$OMP DO COLLAPSE(2)
      do j=1,me%nwan
         do i=1,me%nwan
            ! call me%DressCrvec(me%ham_r(:,i,j), gradH_R)
            call GetGradiant(me%crvec, me%ham_r(:,i,j), gradH_R)
            do idir=1,3
               ! call me%DressPhase(Ar, gradH_R(:,idir))
               call ApplyPhaseFactor([me%nx,me%ny,me%nz], Ar, gradH_R(:,idir))
               call Smooth2Dense_3d(me%nx, me%ny, me%nz, me%nkx, me%nky, me%nkz, gradH_R(:,idir), work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               work_1d = reshape(work_k, [me%nkpts])
               do ik=1,me%nkpts
                  GradHk(i,j,idir,ik) = work_1d(ik)
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(work_r,work_k,work_1d)
      deallocate(gradH_R)
      !$OMP END PARALLEL

   end subroutine GetGradHam_Dressed_3d
!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------
   subroutine GetHam_Dressed_1d(me,Ar,Hk)
      class(wann_fft_t) :: me
      real(dp),intent(in)       :: Ar(3)
      complex(dp),intent(inout) :: Hk(:,:,:)
      integer :: i,j,ik
      complex(dp),allocatable :: HA_r(:),work_r(:),work_k(:)

      !$OMP PARALLEL PRIVATE(i,j,ik,HA_r,work_r,work_k)
      allocate(HA_r(me%nx),work_r(me%nkx),work_k(me%nkx))
      do j=1,me%nwan
         do i=1,me%nwan
            HA_r = me%ham_r(:,i,j)
            ! call me%DressPhase(Ar, HA_r)
            call ApplyPhaseFactor([me%nx], Ar, HA_r)
            call Smooth2Dense_1d(me%nx, me%nkx, HA_r, work_r)
            call dfftw_execute_dft(me%plan_bw,work_r,work_k)
            do ik=1,me%nkx
               Hk(i,j,ik) = work_k(ik)
            end do
         end do
      end do
      deallocate(HA_r,work_r,work_k)
      !$OMP END PARALLEL

   end subroutine GetHam_Dressed_1d
!--------------------------------------------------------------------------------------
   subroutine GetHam_Dressed_2d(me,Ar,Hk)
      class(wann_fft_t) :: me
      real(dp),intent(in)       :: Ar(3)
      complex(dp),intent(inout) :: Hk(:,:,:)
      integer :: i,j,ik
      complex(dp),allocatable :: HA_r(:)
      complex(dp),allocatable :: work_r(:,:),work_k(:,:),work_1d(:)

      !$OMP PARALLEL PRIVATE(i,j,ik,HA_r,work_r,work_k,work_1d)
      allocate(HA_r(me%nrpts))
      allocate(work_r(me%nkx,me%nky),work_k(me%nkx,me%nky),work_1d(me%nkpts))
      !$OMP DO COLLAPSE(2)
      do j=1,me%nwan
         do i=1,me%nwan
            HA_r = me%ham_r(:,i,j)
            ! call me%DressPhase(Ar, HA_r)
            call ApplyPhaseFactor([me%nx,me%ny], Ar, HA_r)
            call Smooth2Dense_2d(me%nx, me%ny, me%nkx, me%nky, HA_r, work_r)
            call dfftw_execute_dft(me%plan_bw,work_r,work_k)
            work_1d = reshape(work_k, [me%nkpts])
            do ik=1,me%nkpts
               Hk(i,j,ik) = work_1d(ik)
            end do
         end do
      end do
      !$OMP END DO
      deallocate(HA_r)
      deallocate(work_r,work_k,work_1d)
      !$OMP END PARALLEL

   end subroutine GetHam_Dressed_2d
!--------------------------------------------------------------------------------------
   subroutine GetHam_Dressed_3d(me,Ar,Hk)
      class(wann_fft_t) :: me
      real(dp),intent(in)       :: Ar(3)
      complex(dp),intent(inout) :: Hk(:,:,:)
      integer :: i,j,ik
      complex(dp),allocatable :: HA_r(:)
      complex(dp),allocatable :: work_r(:,:,:),work_k(:,:,:),work_1d(:)

      !$OMP PARALLEL PRIVATE(i,j,HA_r,work_r,work_k,work_1d)
      allocate(HA_r(me%nrpts))
      allocate(work_r(me%nkx,me%nky,me%nkz),work_k(me%nkx,me%nky,me%nkz),work_1d(me%nkpts))
      !$OMP DO COLLAPSE(2)
      do j=1,me%nwan
         do i=1,me%nwan
            HA_r = me%ham_r(:,i,j)
            call ApplyPhaseFactor([me%nx,me%ny,me%nz], Ar, HA_r)
            call Smooth2Dense_3d(me%nx, me%ny, me%nz, me%nkx, me%nky, me%nkz, HA_r, work_r)
            call dfftw_execute_dft(me%plan_bw,work_r,work_k)
            work_1d = reshape(work_k, [me%nkpts])
            do ik=1,me%nkpts
               Hk(i,j,ik) = work_1d(ik)
            end do
         end do
      end do
      !$OMP END DO
      deallocate(HA_r)
      deallocate(work_r,work_k,work_1d)
      !$OMP END PARALLEL

   end subroutine GetHam_Dressed_3d
!--------------------------------------------------------------------------------------
   subroutine GetDipole_1d(me,Dk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: Dk(:,:,:,:)
      integer :: i,j,idir,ik
      complex(dp),allocatable :: work_r(:),work_k(:)

      !$OMP PARALLEL PRIVATE(i,j,idir,ik,work_r,work_k)
      allocate(work_r(me%nkx),work_k(me%nkx))
      !$OMP DO COLLAPSE(3)
      do idir=1,3
         do j=1,me%nwan
            do i=1,me%nwan
               call Smooth2Dense_1d(me%nx, me%nkx, me%pos_r(:,i,j,idir), work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               do ik=1,me%nkx
                  Dk(i,j,idir,ik) = work_k(ik)
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(work_r,work_k)
      !$OMP END PARALLEL

   end subroutine GetDipole_1d
!--------------------------------------------------------------------------------------
   subroutine GetDipole_2d(me,Dk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: Dk(:,:,:,:)
      integer :: i,j,idir,ik
      complex(dp),allocatable :: work_r(:,:),work_k(:,:),work_1d(:)

      !$OMP PARALLEL PRIVATE(i,j,ik,work_r,work_k,work_1d)
      allocate(work_r(me%nkx,me%nky),work_k(me%nkx,me%nky),work_1d(me%nkpts))
      !$OMP DO COLLAPSE(3)
      do idir=1,3
         do j=1,me%nwan
            do i=1,me%nwan
               call Smooth2Dense_2d(me%nx, me%ny, me%nkx, me%nky, me%pos_r(:,i,j,idir), work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               work_1d = reshape(work_k, [me%nkpts])
               do ik=1,me%nkpts
                  Dk(i,j,idir,ik) = work_1d(ik)
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(work_r,work_k,work_1d)
      !$OMP END PARALLEL

   end subroutine GetDipole_2d
!--------------------------------------------------------------------------------------
   subroutine GetDipole_3d(me,Dk)
      class(wann_fft_t) :: me
      complex(dp),intent(inout) :: Dk(:,:,:,:)
      integer :: i,j,idir,ik
      complex(dp),allocatable :: work_r(:,:,:),work_k(:,:,:),work_1d(:)

      !$OMP PARALLEL PRIVATE(i,j,work_r,work_k,work_1d)
      allocate(work_r(me%nkx,me%nky,me%nkz),work_k(me%nkx,me%nky,me%nkz),work_1d(me%nkpts))
      !$OMP DO COLLAPSE(3)
      do idir=1,3
         do j=1,me%nwan
            do i=1,me%nwan
               call Smooth2Dense_3d(me%nx, me%ny, me%nz, me%nkx, me%nky, me%nkz, me%pos_r(:,i,j,idir), work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               work_1d = reshape(work_k, [me%nkpts])
               do ik=1,me%nkpts
                  Dk(i,j,idir,ik) = work_1d(ik)
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(work_r,work_k,work_1d)
      !$OMP END PARALLEL

   end subroutine GetDipole_3d
!--------------------------------------------------------------------------------------
   subroutine GetDipole_Dressed_1d(me,Ar,Dk)
      class(wann_fft_t) :: me
      real(dp),intent(in)       :: Ar(3)
      complex(dp),intent(inout) :: Dk(:,:,:,:)
      integer :: i,j,idir,ik
      complex(dp),allocatable :: DA_r(:)
      complex(dp),allocatable :: work_r(:),work_k(:)

      !$OMP PARALLEL PRIVATE(i,j,idir,ik,DA_r,work_r,work_k)
      allocate(DA_r(me%nrpts))
      allocate(work_r(me%nkx),work_k(me%nkx))
      !$OMP DO COLLAPSE(3)
      do idir=1,3
         do j=1,me%nwan
            do i=1,me%nwan
               DA_r = me%pos_r(:,i,j,idir)
               call ApplyPhaseFactor([me%nx], Ar, DA_r)
               ! call me%DressPhase(Ar,DA_r)
               call Smooth2Dense_1d(me%nx, me%nkx, DA_r, work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               do ik=1,me%nkx
                  Dk(i,j,idir,ik) = work_k(ik)
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(DA_r)
      deallocate(work_r,work_k)
      !$OMP END PARALLEL

   end subroutine GetDipole_Dressed_1d
!--------------------------------------------------------------------------------------
   subroutine GetDipole_Dressed_2d(me,Ar,Dk)
      class(wann_fft_t) :: me
      real(dp),intent(in)       :: Ar(3)
      complex(dp),intent(inout) :: Dk(:,:,:,:)
      integer :: i,j,idir,ik
      complex(dp),allocatable :: DA_r(:)
      complex(dp),allocatable :: work_r(:,:),work_k(:,:),work_1d(:)

      !$OMP PARALLEL PRIVATE(i,j,idir,ik,DA_r,work_r,work_k,work_1d)
      allocate(DA_r(me%nrpts))
      allocate(work_r(me%nkx,me%nky),work_k(me%nkx,me%nky),work_1d(me%nkpts))
      !$OMP DO COLLAPSE(3)
      do idir=1,3
         do j=1,me%nwan
            do i=1,me%nwan
               DA_r = me%pos_r(:,i,j,idir)
               call ApplyPhaseFactor([me%nx,me%ny], Ar, DA_r)
               ! call me%DressPhase(Ar,DA_r)
               call Smooth2Dense_2d(me%nx, me%ny, me%nkx, me%nky, DA_r, work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               work_1d = reshape(work_k, [me%nkpts])
               do ik=1,me%nkpts
                  Dk(i,j,idir,ik) = work_1d(ik)
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(DA_r)
      deallocate(work_r,work_k,work_1d)
      !$OMP END PARALLEL

   end subroutine GetDipole_Dressed_2d
!--------------------------------------------------------------------------------------
   subroutine GetDipole_Dressed_3d(me,Ar,Dk)
      class(wann_fft_t) :: me
      real(dp),intent(in)       :: Ar(3)
      complex(dp),intent(inout) :: Dk(:,:,:,:)
      integer :: i,j,idir,ik
      complex(dp),allocatable :: DA_r(:)
      complex(dp),allocatable :: work_r(:,:,:),work_k(:,:,:),work_1d(:)

      !$OMP PARALLEL PRIVATE(i,j,idir,ik,DA_r,work_r,work_k,work_1d)
      allocate(DA_r(me%nrpts))
      allocate(work_r(me%nkx,me%nky,me%nkz),work_k(me%nkx,me%nky,me%nkz),work_1d(me%nkpts))
      !$OMP DO COLLAPSE(3)
      do idir=1,3
         do j=1,me%nwan
            do i=1,me%nwan
               DA_r = me%pos_r(:,i,j,idir)
               call ApplyPhaseFactor([me%nx,me%ny,me%nz], Ar, DA_r)
               ! call me%DressPhase(Ar,DA_r)
               call Smooth2Dense_3d(me%nx, me%ny, me%nz, me%nkx, me%nky, me%nkz, DA_r, work_r)
               call dfftw_execute_dft(me%plan_bw,work_r,work_k)
               work_1d = reshape(work_k, [me%nkpts])
               do ik=1,me%nkpts
                  Dk(i,j,idir,ik) = work_1d(ik)
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(DA_r)
      deallocate(work_r,work_k,work_1d)
      !$OMP END PARALLEL

   end subroutine GetDipole_Dressed_3d

!--------------------------------------------------------------------------------------
   subroutine Smooth2Dense_1d(nx,nkx,psi_s,psi_d)
      integer,intent(in) :: nx
      integer,intent(in) :: nkx
      complex(dp),intent(in) :: psi_s(:)
      complex(dp),intent(inout) :: psi_d(:)

      psi_d = zero
      psi_d(1:nx/2) = psi_s(1:nx/2)
      psi_d(nkx - nx/2 + 1:nkx) = psi_s(nx/2 + 1:nx)

   end subroutine Smooth2Dense_1d
!--------------------------------------------------------------------------------------
   subroutine Smooth2Dense_2d(nx,ny,nkx,nky,psi_s,psi_d)
      integer,intent(in) :: nx,ny
      integer,intent(in) :: nkx,nky
      complex(dp),target,intent(in) :: psi_s(:)
      complex(dp),intent(inout) :: psi_d(:,:)
      complex(dp),pointer :: work(:,:)

      work(1:nx, 1:ny) => psi_s

      psi_d = zero

      psi_d( 1:nx/2 , 1:ny/2 ) = work( 1:nx/2 , 1:ny/2 )

      psi_d( nkx - nx/2 + 1:nkx, 1:ny/2 ) = &
         work( nx/2 + 1:nx, 1:ny/2  )

      psi_d( 1:nx/2, nky - ny/2 + 1:nky) = &
         work( 1:nx/2, ny/2 + 1:ny  )

      psi_d( nkx - nx/2 + 1:nkx, nky - ny/2 + 1:nky) = &
         work( nx/2 + 1:nx, ny/2 + 1:ny   )         

   end subroutine Smooth2Dense_2d
!--------------------------------------------------------------------------------------
   subroutine Smooth2Dense_3d(nx,ny,nz,nkx,nky,nkz,psi_s,psi_d)
      integer,intent(in) :: nx,ny,nz
      integer,intent(in) :: nkx,nky,nkz
      complex(dp),target,intent(in) :: psi_s(:)
      complex(dp),intent(inout) :: psi_d(:,:,:)
      complex(dp),pointer :: work(:,:,:)

      work(1:nx, 1:ny, 1:nz) => psi_s

      psi_d = zero

      psi_d( 1:nx/2 , 1:ny/2, 1:nz/2) = work( 1:nx/2 , 1:ny/2, 1:nz/2)

      psi_d( nkx - nx/2 + 1:nkx, 1:ny/2, 1:nz/2 ) = &
         work( nx/2 + 1:nx, 1:ny/2, 1:nz/2  )

      psi_d( 1:nx/2, nky - ny/2 + 1:nky, 1:nz/2) = &
         work( 1:nx/2, ny/2 + 1:ny, 1:nz/2  )

      psi_d( nkx - nx/2 + 1:nkx, nky - ny/2 + 1:nky, 1:nz/2) = &
         work( nx/2 + 1:nx, ny/2 + 1:ny, 1:nz/2   )         

      psi_d( 1:nx/2 , 1:ny/2, nkz - nz/2 + 1:nkz) = work( 1:nx/2 , 1:ny/2, nz/2 + 1 : nz)

      psi_d( nkx - nx/2 + 1:nkx, 1:ny/2, nkz - nz/2 + 1:nkz ) = &
         work( nx/2 + 1:nx, 1:ny/2, nz/2 + 1 : nz  )

      psi_d( 1:nx/2, nky - ny/2 + 1:nky, nkz - nz/2 + 1:nkz) = &
         work( 1:nx/2, ny/2 + 1:ny, nz/2 + 1 : nz  )

      psi_d( nkx - nx/2 + 1:nkx, nky - ny/2 + 1:nky, nkz - nz/2 + 1:nkz) = &
         work( nx/2 + 1:nx, ny/2 + 1:ny, nz/2 + 1 : nz   )   

   end subroutine Smooth2Dense_3d
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

end module wan_fft_ham
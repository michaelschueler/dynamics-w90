module wan_fft_ham
!! Provides tools for reading/writing the Wannier Hamiltonian and for computing 
!! band structures, observables, and Berry-phase properties.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
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
      complex(dp),allocatable,dimension(:,:,:)   :: ham_r
      complex(dp),allocatable,dimension(:,:,:,:) :: pos_r
      ! .. FFTW ..
      integer(kind=8) :: plan_fw,plan_bw
   contains
      procedure, public   :: InitFromW90
      procedure, public   :: Clean
      procedure, public   :: GetHam
      procedure, private  :: GetHam_1d
      procedure, private  :: GetHam_2d
      procedure, private  :: GetHam_3d
   end type wann_fft_t

   type :: latt_t
      integer :: ndim
      integer :: n1,n2,n3
      integer :: npts
      integer,allocatable :: itab(:,:)
   contains
      procedure,public :: Init => latt_Init      
      procedure,public :: Clean => latt_Clean
   end type latt_t
!--------------------------------------------------------------------------------------
contains

!--------------------------------------------------------------------------------------
   subroutine latt_Init(me,nr)
      class(latt_t) :: me
      integer,intent(in) :: nr(:)
      integer :: i1,i2,i3
      integer,allocatable :: indx_2d(:,:,:),indx_3d(:,:,:,:)

      me%ndim = size(nr)
      select case(me%ndim)
      case(1)
         me%n1 = nr(1)
         me%n2 = 1
         me%n3 = 1
         me%npts = me%n1 * me%n2 * me%n3
         allocate(me%itab(me%npts,1))
         do i1=1,nr(1)
            me%itab(i1,1) = FFT_index(me%n1, i1)
         end do
      case(2)
         me%n1 = nr(1)
         me%n2 = nr(2)
         me%n3 = 1
         me%npts = me%n1 * me%n2 * me%n3
         allocate(me%itab(me%npts,2))     
         allocate(indx_2d(me%n1,me%n2,2))
         do concurrent(i2=1:me%n2, i1=1:me%n1)
            indx_2d(i1,i2,1) = FFT_index(me%n1, i1)
            indx_2d(i1,i2,2) = FFT_index(me%n2, i2)
         end do     
         me%itab(:,1) = reshape(indx_2d(:,:,1), [me%npts])
         me%itab(:,2) = reshape(indx_2d(:,:,2), [me%npts])
         deallocate(indx_2d)
      case(3)
         me%n1 = nr(1)
         me%n2 = nr(2)
         me%n3 = nr(3)
         me%npts = me%n1 * me%n2 * me%n3
         allocate(me%itab(me%npts,3))     
         allocate(indx_3d(me%n1,me%n2,me%n3,3))
         do concurrent(i3=1:me%n3, i2=1:me%n2, i1=1:me%n1)
            indx_3d(i1,i2,i3,1) = FFT_index(me%n1, i1)
            indx_3d(i1,i2,i3,2) = FFT_index(me%n2, i2)
            indx_3d(i1,i2,i3,3) = FFT_index(me%n3, i3)
         end do     
         me%itab(:,1) = reshape(indx_3d(:,:,:,1), [me%npts])
         me%itab(:,2) = reshape(indx_3d(:,:,:,2), [me%npts])
         me%itab(:,3) = reshape(indx_3d(:,:,:,3), [me%npts])
         deallocate(indx_3d)
      end select

   end subroutine latt_Init
!--------------------------------------------------------------------------------------
   subroutine latt_Clean(me)
      class(latt_t) :: me

      deallocate(me%itab)

   end subroutine latt_Clean
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
      type(latt_t) :: lat

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
         allocate(zr_1d(me%nkx),zk_1d(me%nkx))
         call DFFTW_PLAN_DFT_1D(me%plan_fw,me%nkx,zk_1d,zr_1d,FFTW_FORWARD,FFTW_MEASURE)
         call DFFTW_PLAN_DFT_1D(me%plan_bw,me%nkx,zr_1d,zk_1d,FFTW_BACKWARD,FFTW_MEASURE)
         deallocate(zr_1d,zk_1d)

         call lat%Init([me%nx])
      case(2)
         me%nx = ix_bound(2) - ix_bound(1) + 1
         me%ny = iy_bound(2) - iy_bound(1) + 1
         if(mod(me%nx,2) /= 0) me%nx = me%nx + 1
         if(mod(me%ny,2) /= 0) me%ny = me%ny + 1
         me%nz = 1

         me%nkx = nk(1)
         me%nky = nk(2)
         allocate(zr_2d(me%nkx,me%nky),zk_2d(me%nkx,me%nky))
         call DFFTW_PLAN_DFT_2D(me%plan_fw,me%nkx,me%nky,zk_2d,zr_2d,FFTW_FORWARD,FFTW_MEASURE)
         call DFFTW_PLAN_DFT_2D(me%plan_bw,me%nkx,me%nky,zr_2d,zk_2d,FFTW_BACKWARD,FFTW_MEASURE)
         deallocate(zr_2d,zk_2d)

         call lat%Init([me%nx,me%ny])
      case(3)
         me%nx = ix_bound(2) - ix_bound(1) + 1
         me%ny = iy_bound(2) - iy_bound(1) + 1
         me%nz = iz_bound(2) - iz_bound(1) + 1
         if(mod(me%nx,2) /= 0) me%nx = me%nx + 1
         if(mod(me%nx,2) /= 0) me%nx = me%nx + 1
         if(mod(me%ny,2) /= 0) me%ny = me%ny + 1

         me%nkx = nk(1)
         me%nky = nk(2)
         me%nkz = nk(3)
         allocate(zr_3d(me%nkx,me%nky,me%nz),zk_3d(me%nkx,me%nky,me%nz))
         call DFFTW_PLAN_DFT_3D(me%plan_fw,me%nkx,me%nky,me%nkz,zk_3d,zr_3d,FFTW_FORWARD,FFTW_MEASURE)
         call DFFTW_PLAN_DFT_3D(me%plan_bw,me%nkx,me%nky,me%nkz,zr_3d,zk_3d,FFTW_BACKWARD,FFTW_MEASURE)
         deallocate(zr_3d,zk_3d)

         call lat%Init([me%nx,me%ny,me%nz])
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
         do ir=1,lat%Npts
            ix = lat%itab(ir,1)
            iy = lat%itab(ir,2)
            irw = minloc( abs(w90%irvec(:,1) - ix) + abs(w90%irvec(:,2) - iy), dim=1)
            if((w90%irvec(irw,1) == ix) .and. (w90%irvec(irw,2) == iy)) then
               w90_rindx(ir) = irw
            end if
         end do

         ! allocate(indx_2d(me%nx,me%ny)); indx_2d = -999
         ! do ir=1,w90%nrpts
         !    ix = GetFFTIndex(me%nx, w90%irvec(ir,1))
         !    iy = GetFFTIndex(me%ny, w90%irvec(ir,2))
         !    indx_2d(ix,iy) = ir
         ! end do         
         ! w90_rindx = reshape(indx_2d, [me%nrpts])  
         ! ! deallocate(indx_2d)
         ! do ix=1,me%nx
         !    if(ix <= me%nx/2) then
         !       kx = ix
         !    else
         !       kx = -(me%nx - ix + 1)

         !    do iy=1,me%nky

         !    end do
         ! end do
      case(3)
         allocate(indx_3d(me%nx,me%ny,me%nz)); indx_2d = -999
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

      do ir=1,me%nrpts
         irw = w90_rindx(ir)
         if(irw == -999) cycle
         me%ndegen(ir) = w90%ndegen(irw)
         do j=1,me%nwan
            do i=1,me%nwan
               ! me%ham_r(ir,i,j) = w90%ham_r(i,j,irw) / me%ndegen(ir)
               ! me%pos_r(ir,i,j,:) = w90%pos_r(i,j,:,irw) / me%ndegen(ir)
               me%ham_r(ir,i,j) = w90%ham_r(i,j,irw)  / me%ndegen(ir)
            end do
         end do
      end do

      deallocate(w90_rindx)
      call lat%Clean()

   end subroutine InitFromW90
!--------------------------------------------------------------------------------------
   subroutine Clean(me)
      class(wann_fft_t) :: me

      if(allocated(me%ham_r)) deallocate(me%ham_r)
      if(allocated(me%pos_r)) deallocate(me%pos_r)

   end subroutine Clean
!--------------------------------------------------------------------------------------
   subroutine GetHam(me,Hk)
      class(wann_fft_t),intent(in) :: me
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
   subroutine GetHam_1d(me,Hk)
      class(wann_fft_t),intent(in) :: me
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
      class(wann_fft_t),intent(in) :: me
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
      class(wann_fft_t),intent(in) :: me
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

      psi_d( 1:nx/2 , 1:ny/2, 1:nz/2) = work( 1:nx/2 , 1:ny/2, nz/2 + 1 : nz)

      psi_d( nkx - nx/2 + 1:nkx, 1:ny/2, nz/2 + 1 : nz ) = &
         work( nx/2 + 1:nx, 1:ny/2, nz/2 + 1 : nz  )

      psi_d( 1:nx/2, nky - ny/2 + 1:nky, 1:nz/2) = &
         work( 1:nx/2, ny/2 + 1:ny, nz/2 + 1 : nz  )

      psi_d( nkx - nx/2 + 1:nkx, nky - ny/2 + 1:nky, 1:nz/2) = &
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
   pure elemental function FFT_index(n,i) result(k)
      integer,intent(in) :: n, i
      integer :: k

      if(i <= n/2) then
         k = i - 1
      else
         k = -(n - i + 1)
      end if

   end function FFT_index
!--------------------------------------------------------------------------------------
   



!======================================================================================

end module wan_fft_ham
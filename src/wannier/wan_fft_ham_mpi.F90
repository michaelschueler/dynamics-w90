module wan_fft_ham_mpi
!! Provides tools for reading/writing the Wannier Hamiltonian and for computing 
!! band structures, observables, and Berry-phase properties.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use,intrinsic :: ISO_C_binding
   use mpi
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
   use scitools_utils,only: str,stop_error
   use scitools_array1d_dist,only: dist_array1d_t
   use wan_hamiltonian,only: wann90_tb_t
   implicit none
   include '../units_inc.f90'
   include '../formats.h'
   include 'fftw3.f03'
!--------------------------------------------------------------------------------------
   private
   public :: wann_fft_t
!--------------------------------------------------------------------------------------
   integer :: ntasks,taskid,ierr
!--------------------------------------------------------------------------------------
   type :: wann_fft_t
      integer                                    :: nwan,nwan2
      integer                                    :: nx,ny,nz,nrpts
      integer                                    :: nkx,nky,nkz,nkpts,nkpts_loc,kdim
      integer,allocatable,dimension(:)           :: ndegen
      real(dp),allocatable,dimension(:,:)        :: crvec
      complex(dp),allocatable,dimension(:)       :: work_r,work_k,HA_r
      complex(dp),allocatable,dimension(:,:)     :: gradH_R
      complex(dp),allocatable,dimension(:,:,:)   :: ham_r
      complex(dp),allocatable,dimension(:,:,:,:) :: pos_r
      ! .. FFTW ..
      integer(kind=8) :: plan_fw,plan_bw
   contains
      procedure, public   :: InitFromW90
      procedure, public   :: Clean
      procedure, public   :: GetHam
      procedure, public   :: GetGradHam
      procedure, public   :: GetDipole
   end type wann_fft_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine InitFromW90(me,w90,nk,kdist)
      class(wann_fft_t),target :: me
      type(wann90_tb_t),intent(in) :: w90
      integer,intent(in) :: nk(:)
      type(dist_array1d_t),intent(in) :: kdist
      integer :: ir,irw,ix,iy,iz,i,j
      integer :: kx,ky
      integer :: ix_bound(2),iy_bound(2),iz_bound(2)
      integer :: ierr,tid
      integer,allocatable :: w90_rindx(:)
      integer,allocatable :: indx_2d(:,:), indx_3d(:,:,:)
      
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

      me%nwan = w90%num_wann
      me%nwan2 = (me%nwan * (me%nwan + 1)) / 2
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
            call MPI_Finalize(ierr); stop
         end if

         allocate(me%work_r(me%nkx),me%work_k(me%nkx))

         call DFFTW_PLAN_DFT_1D(me%plan_fw,me%nkx,me%work_k,me%work_r,FFTW_FORWARD,FFTW_MEASURE)
         call DFFTW_PLAN_DFT_1D(me%plan_bw,me%nkx,me%work_r,me%work_k,FFTW_BACKWARD,FFTW_MEASURE)
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
            call MPI_Finalize(ierr); stop
         end if

         allocate(me%work_r(me%nkx*me%nky),me%work_k(me%nkx*me%nky))

         call DFFTW_PLAN_DFT_2D(me%plan_fw,me%nkx,me%nky,me%work_k,me%work_r,FFTW_FORWARD,FFTW_MEASURE)
         call DFFTW_PLAN_DFT_2D(me%plan_bw,me%nkx,me%nky,me%work_r,me%work_k,FFTW_BACKWARD,FFTW_MEASURE)
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
            call MPI_Finalize(ierr); stop
         end if

         allocate(me%work_r(me%nkx*me%nky*me%nkz),me%work_k(me%nkx*me%nky*me%nkz))
         call DFFTW_PLAN_DFT_3D(me%plan_fw,me%nkx,me%nky,me%nkz,me%work_k,me%work_r,FFTW_FORWARD,FFTW_MEASURE)
         call DFFTW_PLAN_DFT_3D(me%plan_bw,me%nkx,me%nky,me%nkz,me%work_r,me%work_k,FFTW_BACKWARD,FFTW_MEASURE)

      end select

      me%nrpts = me%nx * me%ny * me%nz
      me%nkpts = me%nkx * me%nky * me%nkz
      me%nkpts_loc = kdist%N_loc(taskid)

      allocate(me%HA_r(me%nrpts))
      allocate(me%gradH_R(me%nrpts,3))

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
      if(allocated(me%work_r)) deallocate(me%work_r)
      if(allocated(me%work_k)) deallocate(me%work_k)
      if(allocated(me%HA_r)) deallocate(me%HA_r)
      if(allocated(me%gradH_R)) deallocate(me%gradH_R)

   end subroutine Clean
!--------------------------------------------------------------------------------------
   subroutine GetHam(me,kdist,Hk,Ar)
      class(wann_fft_t) :: me
      type(dist_array1d_t),intent(in) :: kdist
      complex(dp),intent(inout) :: Hk(:,:,:)
      real(dp),intent(in),optional :: Ar(3)
      integer :: i,j,ik,ik_glob
      integer,allocatable :: rdims(:), kdims(:)

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

      if(present(Ar)) then
         do j=1,me%nwan
            do i=1,j
               me%HA_r = me%ham_r(:,i,j)
               call ApplyPhaseFactor(rdims, Ar, me%HA_r)
               call Smooth2Dense(rdims, kdims, me%HA_r, me%work_r)
               call dfftw_execute_dft(me%plan_bw,me%work_r,me%work_k)
               do ik=1,me%nkpts_loc
                  ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
                  Hk(i,j,ik) = me%work_k(ik_glob)
                  if(i < j) Hk(j,i,ik) = conjg(Hk(i,j,ik))
               end do
            end do
         end do
      else
         do j=1,me%nwan
            do i=1,j
               call Smooth2Dense(rdims, kdims, me%ham_r(:,i,j), me%work_r)
               call dfftw_execute_dft(me%plan_bw,me%work_r,me%work_k)
               do ik=1,me%nkpts_loc
                  ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
                  Hk(i,j,ik) = me%work_k(ik_glob)
                  if(i < j) Hk(j,i,ik) = conjg(Hk(i,j,ik))
               end do
            end do
         end do
      end if

      deallocate(rdims,kdims)

   end subroutine GetHam
!--------------------------------------------------------------------------------------
   subroutine GetGradHam(me,kdist,GradHk,Ar)
      class(wann_fft_t) :: me
      type(dist_array1d_t),intent(in) :: kdist
      complex(dp),intent(inout) :: GradHk(:,:,:,:)
      real(dp),intent(in),optional :: Ar(3)
      integer :: i,j,idir,ik,ik_glob
      integer,allocatable :: rdims(:), kdims(:)

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

      if(present(Ar)) then
         do j=1,me%nwan
            do i=1,j
               call GetGradiant(me%crvec, me%ham_r(:,i,j), me%gradH_R)
               do idir=1,3
                  call ApplyPhaseFactor(rdims, Ar, me%gradH_R(:,idir))
                  call Smooth2Dense(rdims, kdims, me%gradH_R(:,idir), me%work_r)
                  call dfftw_execute_dft(me%plan_bw,me%work_r,me%work_k)
                  do ik=1,me%nkpts_loc
                     ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
                     GradHk(i,j,idir,ik) = me%work_k(ik_glob)
                     if(i < j) GradHk(j,i,idir,ik) = conjg(GradHk(i,j,idir,ik))
                  end do
               end do
            end do
         end do
      else
         do j=1,me%nwan
            do i=1,j
               call GetGradiant(me%crvec, me%ham_r(:,i,j), me%gradH_R)
               do idir=1,3
                  call Smooth2Dense(rdims, kdims, me%gradH_R(:,idir), me%work_r)
                  call dfftw_execute_dft(me%plan_bw,me%work_r,me%work_k)
                  do ik=1,me%nkpts_loc
                     ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
                     GradHk(i,j,idir,ik) = me%work_k(ik_glob)
                     if(i < j) GradHk(j,i,idir,ik) = conjg(GradHk(i,j,idir,ik))
                  end do
               end do
            end do
         end do
      end if

      deallocate(rdims,kdims)

   end subroutine GetGradHam
!--------------------------------------------------------------------------------------
     subroutine GetDipole(me,kdist,Dk,Ar)
      class(wann_fft_t) :: me
      type(dist_array1d_t),intent(in) :: kdist
      complex(dp),intent(inout) :: Dk(:,:,:,:)
      real(dp),intent(in),optional :: Ar(3)
      integer :: i,j,idir,ik,ik_glob
      integer,allocatable :: rdims(:), kdims(:)


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

      if(present(Ar)) then
         do j=1,me%nwan
            do i=1,j
               do idir=1,3
                  me%HA_r = me%pos_r(:,i,j,idir)
                  call ApplyPhaseFactor(rdims, Ar, me%HA_r)
                  call Smooth2Dense(rdims, kdims, me%HA_r, me%work_r)
                  call dfftw_execute_dft(me%plan_bw,me%work_r,me%work_k)
                  do ik=1,me%nkpts_loc
                     ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
                     Dk(i,j,idir,ik) = me%work_k(ik_glob)
                     if(i < j) Dk(j,i,idir,ik) = conjg(Dk(i,j,idir,ik))
                  end do
               end do
            end do
         end do
      else
         do j=1,me%nwan
            do i=1,j
               do idir=1,3
                  call Smooth2Dense(rdims, kdims, me%pos_r(:,i,j,idir), me%work_r)
                  call dfftw_execute_dft(me%plan_bw,me%work_r,me%work_k)
                  do ik=1,me%nkpts_loc
                     ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
                     Dk(i,j,idir,ik) = me%work_k(ik_glob)
                     if(i < j) Dk(j,i,idir,ik) = conjg(Dk(i,j,idir,ik))
                  end do
               end do
            end do
         end do
      end if

      deallocate(rdims,kdims)

   end subroutine GetDipole
!--------------------------------------------------------------------------------------
   subroutine Smooth2Dense_1d(nx,nkx,psi_s,psi_d)
      integer,intent(in) :: nx
      integer,intent(in) :: nkx
      complex(dp),intent(in) :: psi_s(:)
      complex(dp),target,intent(inout) :: psi_d(:)

      psi_d = zero
      psi_d(1:nx/2) = psi_s(1:nx/2)
      psi_d(nkx - nx/2 + 1:nkx) = psi_s(nx/2 + 1:nx)

   end subroutine Smooth2Dense_1d
!--------------------------------------------------------------------------------------
   subroutine Smooth2Dense_2d(nx,ny,nkx,nky,psi_s,phi_d)
      integer,intent(in) :: nx,ny
      integer,intent(in) :: nkx,nky
      complex(dp),target,intent(in) :: psi_s(:)
      complex(dp),target,intent(inout) :: phi_d(:)
      complex(dp),pointer :: work(:,:),psi_d(:,:)

      work(1:nx, 1:ny) => psi_s
      psi_d(1:nkx, 1:nky) => phi_d

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
   subroutine Smooth2Dense_3d(nx,ny,nz,nkx,nky,nkz,psi_s,phi_d)
      integer,intent(in) :: nx,ny,nz
      integer,intent(in) :: nkx,nky,nkz
      complex(dp),target,intent(in) :: psi_s(:)
      complex(dp),target,intent(inout) :: phi_d(:)
      complex(dp),pointer :: work(:,:,:),psi_d(:,:,:)


      work(1:nx, 1:ny, 1:nz) => psi_s
      psi_d(1:nkx, 1:nky, 1:nkz) => phi_d

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
   subroutine Smooth2Dense(rdims,kdims,psi_s,psi_d)
      integer,intent(in) :: rdims(:)
      integer,intent(in) :: kdims(:)
      complex(dp),intent(in) :: psi_s(:)
      complex(dp),intent(inout) :: psi_d(:)
      integer :: kdim

      kdim = size(rdims)

      select case(kdim)
      case(1)
         call Smooth2Dense_1d(rdims(1), kdims(1), psi_s, psi_d)
      case(2)
         call Smooth2Dense_2d(rdims(1), rdims(2), kdims(1), kdims(2), psi_s, psi_d)
      case(3)
         call Smooth2Dense_3d(rdims(1), rdims(2), rdims(3), kdims(1), kdims(2), kdims(3),&
          psi_s, psi_d)
      end select

   end subroutine Smooth2Dense
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

end module wan_fft_ham_mpi
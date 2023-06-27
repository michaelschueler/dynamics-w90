!--------------------------------------------------------------------------------------
   subroutine InitFromW90(me,w90,nk,kdist,nthreads)
      class(wann_fft_t),target :: me
      type(wann90_tb_t),intent(in) :: w90
      integer,intent(in) :: nk(:)
      type(dist_array1d_t),intent(in) :: kdist
      integer,intent(in),optional :: nthreads
      integer :: ir,irw,ix,iy,iz,i,j
      integer :: kx,ky
      integer :: ix_bound(2),iy_bound(2),iz_bound(2)
      integer :: ierr,tid
      integer,allocatable :: w90_rindx(:)
      integer,allocatable :: indx_2d(:,:), indx_3d(:,:,:)
      integer(C_INTPTR_T) :: nk1,nk2,nk3
      integer(C_INTPTR_T) :: alloc_local
      logical :: on_root
      
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

      on_root = taskid == 0

      call fftw_mpi_init()

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

         nk1 = me%nkx
         nk2 = me%nky
         nk3 = 1

         alloc_local = fftw_mpi_local_size_2d(nk2, nk1, MPI_COMM_WORLD, me%local_n, me%local_offset)
         me%cdata = fftw_alloc_complex(alloc_local)
         call c_f_pointer(me%cdata, me%data, [nk1 * me%local_n])

         me%plan_fw = fftw_mpi_plan_dft_2d(nk2, nk1, me%data, me%data, MPI_COMM_WORLD, &
                                    FFTW_FORWARD, FFTW_MEASURE)
         me%plan_bw = fftw_mpi_plan_dft_2d(nk2, nk1, me%data, me%data, MPI_COMM_WORLD, &
                                    FFTW_BACKWARD, FFTW_MEASURE)

         me%nkpts_loc = nk1 * me%local_n

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

         nk1 = me%nkx
         nk2 = me%nky
         nk3 = me%nkz

         alloc_local = fftw_mpi_local_size_3d(nk3, nk2, nk1, MPI_COMM_WORLD, me%local_n, me%local_offset)
         me%cdata = fftw_alloc_complex(alloc_local)
         call c_f_pointer(me%cdata, me%data, [nk1 * nk2 * me%local_n])

         me%plan_fw = fftw_mpi_plan_dft_3d(nk3, nk2, nk1, me%data, me%data, MPI_COMM_WORLD, &
                                    FFTW_FORWARD, FFTW_MEASURE)
         me%plan_bw = fftw_mpi_plan_dft_3d(nk3, nk2, nk1, me%data, me%data, MPI_COMM_WORLD, &
                                    FFTW_BACKWARD, FFTW_MEASURE)

         me%nkpts_loc = nk1 * nk2 * me%local_n

      case default
         call stop_error("FFTW-MPI not implemented for dimensions other than 2 or 3", root_flag=taskid==0)
      end select

      if(me%nkpts_loc /= kdist%N_loc(taskid)) then
         call stop_error("Incompatible distributed grid", root_flag=on_root)
      end if

      me%nrpts = me%nx * me%ny * me%nz
      me%nkpts = me%nkx * me%nky * me%nkz

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
      if(allocated(me%HA_r)) deallocate(me%HA_r)
      if(allocated(me%gradH_R)) deallocate(me%gradH_R)

      call fftw_destroy_plan(me%plan_fw)
      call fftw_destroy_plan(me%plan_bw)
      call fftw_free(me%cdata)

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
               call me%ApplyPhaseFactor(Ar, me%HA_r)
               call Smooth2Dense(rdims, kdims, me%local_n, me%local_offset, me%HA_r, me%data)
               call fftw_mpi_execute_dft(me%plan_bw, me%data, me%data)
               do ik=1,me%nkpts_loc
                  Hk(i,j,ik) = me%data(ik)
                  if(i < j) Hk(j,i,ik) = conjg(Hk(i,j,ik))
               end do
            end do
         end do
      else
         do j=1,me%nwan
            do i=1,j
               call Smooth2Dense(rdims, kdims, me%local_n, me%local_offset, me%ham_r(:,i,j), me%data)
               call fftw_mpi_execute_dft(me%plan_bw, me%data, me%data)
               do ik=1,me%nkpts_loc
                  Hk(i,j,ik) = me%data(ik)
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
                  call me%ApplyPhaseFactor(Ar, me%gradH_R(:,idir))
                  call Smooth2Dense(rdims, kdims, me%local_n, me%local_offset, me%gradH_R(:,idir), me%data)
                  call fftw_mpi_execute_dft(me%plan_bw, me%data, me%data)
                  do ik=1,me%nkpts_loc
                     GradHk(i,j,idir,ik) = me%data(ik)
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
                  call Smooth2Dense(rdims, kdims, me%local_n, me%local_offset, me%gradH_R(:,idir), me%data)
                  call fftw_mpi_execute_dft(me%plan_bw, me%data, me%data)
                  do ik=1,me%nkpts_loc
                     GradHk(i,j,idir,ik) = me%data(ik)
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
                  call me%ApplyPhaseFactor(Ar, me%HA_r)
                  call Smooth2Dense(rdims, kdims, me%local_n, me%local_offset, me%HA_r, me%data)
                  call fftw_mpi_execute_dft(me%plan_bw, me%data, me%data)
                  do ik=1,me%nkpts_loc
                     Dk(i,j,idir,ik) = me%data(ik)
                     if(i < j) Dk(j,i,idir,ik) = conjg(Dk(i,j,idir,ik))
                  end do
               end do
            end do
         end do
      else
         do j=1,me%nwan
            do i=1,j
               do idir=1,3
                  call Smooth2Dense(rdims, kdims, me%local_n, me%local_offset, me%pos_r(:,i,j,idir), me%data)
                  call fftw_mpi_execute_dft(me%plan_bw, me%data, me%data)
                  do ik=1,me%nkpts_loc
                     Dk(i,j,idir,ik) = me%data(ik)
                     if(i < j) Dk(j,i,idir,ik) = conjg(Dk(i,j,idir,ik))
                  end do
               end do
            end do
         end do
      end if

      deallocate(rdims,kdims)

   end subroutine GetDipole
!--------------------------------------------------------------------------------------
   subroutine ApplyPhaseFactor(me,Ar,OO_R)
      class(wann_fft_t) :: me
      real(dp),intent(in) :: Ar(3)
      complex(dp),target,intent(inout) :: OO_R(:)
      integer :: ix,iy,iz,kx,ky,kz
      real(dp) :: adot,c,s
      complex(dp),pointer :: OO_1d(:)
      complex(dp),pointer :: OO_2d(:,:)
      complex(dp),pointer :: OO_3d(:,:,:)

      select case(me%kdim)
      case(2)
         OO_2d(1:me%nx, 1:me%ny) => OO_R
         do concurrent(iy=1:me%ny, ix=1:me%nx)
            kx = FFT_Freq(me%nx, ix)
            ky = FFT_Freq(me%ny, iy)
            adot = Ar(1) * kx + Ar(2) * ky
            c = cos(DPI * adot)
            s = sin(DPI * adot)
            OO_2d(ix,iy) = cmplx(c,-s,kind=dp) * OO_2d(ix,iy)
         end do         
      case(3)
         OO_3d(1:me%nx, 1:me%ny, 1:me%nz) => OO_R
         do concurrent(iz=1:me%nz, iy=1:me%ny, ix=1:me%nx)
            kx = FFT_Freq(me%nx, ix)
            ky = FFT_Freq(me%ny, iy)
            kz = FFT_Freq(me%nz, iz)
            adot = Ar(1) * kx + Ar(2) * ky  + Ar(3) * kz
            c = cos(DPI * adot)
            s = sin(DPI * adot)
            OO_3d(ix,iy,iz) = cmplx(c,-s,kind=dp) * OO_3d(ix,iy,iz)
         end do   
      end select
      
   end subroutine ApplyPhaseFactor
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
   subroutine Smooth2Dense_2d(nx,ny,nkx,nky,nloc,offset,psi_s,phi_d)
      integer,intent(in)  :: nx,ny
      integer,intent(in)  :: nkx,nky
      integer(C_INTPTR_T) :: nloc,offset
      integer :: iy,jy,i,j
      complex(dp),target,intent(in) :: psi_s(:)
      complex(dp),target,intent(inout) :: phi_d(:)
      complex(dp),pointer :: work(:,:),psi_d(:,:)

      work(1:nx, 1:ny) => psi_s
      psi_d(1:nkx, 1:nloc) => phi_d

      psi_d = zero

      do iy=1,ny/2
         i = iy - offset
         if(i <= 0 .or. i > nloc) cycle
         psi_d( 1:nx/2 , i ) = work( 1:nx/2 , iy )

         psi_d( nkx - nx/2 + 1:nkx, i ) = work( nx/2 + 1:nx, iy  )
      end do

      do iy=ny/2+1,ny
         jy = nky - ny + iy
         j = jy - offset
         if(j <= 0 .or. j > nloc) cycle

         psi_d( 1:nx/2, j) = work( 1:nx/2, iy  )
         psi_d( nkx - nx/2 + 1:nkx, j) = work( nx/2 + 1:nx, iy )
      end do

   end subroutine Smooth2Dense_2d
!--------------------------------------------------------------------------------------
   subroutine Smooth2Dense_3d(nx,ny,nz,nkx,nky,nkz,nloc,offset,psi_s,phi_d)
      integer,intent(in)  :: nx,ny,nz
      integer,intent(in)  :: nkx,nky,nkz
      integer(C_INTPTR_T) :: nloc,offset
      integer :: iz,jz,i,j
      complex(dp),target,intent(in) :: psi_s(:)
      complex(dp),target,intent(inout) :: phi_d(:)
      complex(dp),pointer :: work(:,:,:),psi_d(:,:,:)


      work(1:nx, 1:ny, 1:nz) => psi_s
      psi_d(1:nkx, 1:nky, 1:nloc) => phi_d

      psi_d = zero

      do iz=1,nz/2
         i = iz - offset
         if(i <= 0 .or. i > nloc) cycle
 
         psi_d( 1:nx/2 , 1:ny/2, i) = work( 1:nx/2 , 1:ny/2, iz)

         psi_d( nkx - nx/2 + 1:nkx, 1:ny/2, i ) = &
            work( nx/2 + 1:nx, 1:ny/2, iz  )

         psi_d( 1:nx/2, nky - ny/2 + 1:nky, i) = &
            work( 1:nx/2, ny/2 + 1:ny, iz  )

         psi_d( nkx - nx/2 + 1:nkx, nky - ny/2 + 1:nky, i) = &
            work( nx/2 + 1:nx, ny/2 + 1:ny, iz   )         

      end do

      do iz=nz/2+1,nz
         jz = nkz - nz + iz
         j = jz - offset
         if(j <= 0 .or. j > nloc) cycle

         psi_d( 1:nx/2 , 1:ny/2, j) = work( 1:nx/2 , 1:ny/2, iz)

         psi_d( nkx - nx/2 + 1:nkx, 1:ny/2, j ) = &
            work( nx/2 + 1:nx, 1:ny/2, iz  )

         psi_d( 1:nx/2, nky - ny/2 + 1:nky, j) = &
            work( 1:nx/2, ny/2 + 1:ny, iz  )

         psi_d( nkx - nx/2 + 1:nkx, nky - ny/2 + 1:nky, j) = &
            work( nx/2 + 1:nx, ny/2 + 1:ny, iz   )   

      end do

   end subroutine Smooth2Dense_3d
!--------------------------------------------------------------------------------------
   subroutine Smooth2Dense(rdims,kdims,nloc,offset,psi_s,psi_d)
      integer,intent(in)  :: rdims(:)
      integer,intent(in)  :: kdims(:)
      integer(C_INTPTR_T) :: nloc,offset
      complex(dp),intent(in) :: psi_s(:)
      complex(dp),intent(inout) :: psi_d(:)
      integer :: kdim

      kdim = size(rdims)

      select case(kdim)
      case(2)
         call Smooth2Dense_2d(rdims(1), rdims(2), kdims(1), kdims(2), nloc, offset, psi_s, psi_d)
      case(3)
         call Smooth2Dense_3d(rdims(1), rdims(2), rdims(3), kdims(1), kdims(2), kdims(3),&
          nloc, offset, psi_s, psi_d)
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

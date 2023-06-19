!--------------------------------------------------------------------------------------
   subroutine InitFromW90(me,w90,nk,kdist,nthreads)
      class(wann_spfft_t),target :: me
      type(wann90_tb_t),intent(in) :: w90
      type(dist_array1d_t),intent(in) :: kdist
      integer,intent(in),optional :: nthreads
      integer,intent(in) :: nk(:)
      integer :: ir,irw,ix,iy,iz,i,j,counter
      integer :: ierr,tid
      integer :: errorCode = 0
      integer,allocatable :: w90_rindx(:)
      integer,allocatable :: indx_2d(:,:), indx_3d(:,:,:)

      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

      if(present(nthreads)) me%maxNumThreads = nthreads

      me%nwan = w90%num_wann
      me%kdim = size(nk)
      me%nrpts = w90%nrpts
      me%nkpts_loc = kdist%N_loc(taskid)

      me%nkx = 1
      me%nky = 1
      me%nkz = 1

      select case(me%kdim)
      case(2)
         me%nkx = nk(1)
         me%nky = nk(2)
      case(3)
         me%nkx = nk(1)
         me%nky = nk(2)
         me%nkz = nk(3)
      case default
         call stop_error("SpFFT not implemented for dimensions other than 2 or 3", root_flag=taskid==0)
      end select       

      me%nkpts = me%nkx * me%nky * me%nkz  

      me%maxNumLocalZColumns = me%nkx * me%nky
      allocate(me%irvec(me%nrpts,3))
      allocate(me%indices(me%nrpts * 3)); me%indices = 0
      allocate(me%spaceDomain(2 * me%nkpts))

      me%irvec = w90%irvec

      do ir=0,w90%nrpts-1
         me%indices(ir * 3 + 1) = w90%irvec(ir+1,1)
         me%indices(ir * 3 + 2) = w90%irvec(ir+1,2)
         me%indices(ir * 3 + 3) = w90%irvec(ir+1,3)
      end do  

      allocate(me%ndegen(me%nrpts)); me%ndegen = 1000000
      allocate(me%ham_r(me%nrpts,me%nwan,me%nwan)); me%ham_r = zero
      allocate(me%pos_r(me%nrpts,me%nwan,me%nwan,3)); me%pos_r = zero
      allocate(me%crvec(me%nrpts,3)); me%crvec = 0.0_dp

      do ir=1,me%nrpts
         me%ndegen(ir) = w90%ndegen(ir)
         me%crvec(ir,:) = w90%crvec(ir,:)
         do j=1,me%nwan
            do i=1,me%nwan
               me%ham_r(ir,i,j) = w90%ham_r(i,j,ir) / me%ndegen(ir)
               me%pos_r(ir,i,j,:) = w90%pos_r(i,j,:,ir) / me%ndegen(ir)
            end do
         end do
      end do


      allocate(me%HA_r(me%nrpts))
      allocate(me%gradH_R(me%nrpts,3))

      ! create grid
      errorCode = spfft_grid_create(me%grid, me%nkx, me%nky, me%nkz, me%maxNumLocalZColumns, &
         me%processingUnit, me%maxNumThreads)
      if (errorCode /= SPFFT_SUCCESS) error stop

      ! create transform
      ! Note: A transform handle can be created without a grid if no resource sharing is desired.
      ! errorCode = spfft_transform_create(me%transform, me%grid, me%processingUnit, 0, &
      !    me%nkx, me%nky, me%nkz, me%nkz, me%nrpts, SPFFT_INDEX_TRIPLETS, me%indices)
      ! if (errorCode /= SPFFT_SUCCESS) error stop

      errorCode = spfft_transform_create(me%transform, me%grid, me%processingUnit, 0, &
         me%nkx, me%nky, me%nkz, me%nkz, me%nrpts, SPFFT_INDEX_TRIPLETS, me%indices)
      if (errorCode /= SPFFT_SUCCESS) error stop

      ! grid can be safely destroyed after creating all required transforms
      errorCode = spfft_grid_destroy(me%grid)
      if (errorCode /= SPFFT_SUCCESS) error stop

      ! set space domain array to use memory allocted by the library
      errorCode = spfft_transform_get_space_domain(me%transform, me%processingUnit, me%realValuesPtr)
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
      if(allocated(me%spaceDomain)) deallocate(me%spaceDomain)

   end subroutine Clean
!--------------------------------------------------------------------------------------
   subroutine GetHam(me,kdist,Hk,Ar)
      class(wann_spfft_t) :: me
      type(dist_array1d_t),intent(in) :: kdist
      complex(dp),intent(inout) :: Hk(:,:,:)
      real(dp),intent(in),optional :: Ar(3)
      integer :: i,j,ik,ik_glob
      integer :: errorCode = 0
      complex(C_DOUBLE_COMPLEX), pointer :: spaceDomainPtr(:)

      call assert_shape(Hk, [me%nwan,me%nwan,me%nkpts_loc], "GetHam", "Hk")

      do j=1,me%nwan
         do i=1,j
            me%HA_r = me%ham_r(:,i,j)
            if(present(Ar)) then
               call me%ApplyPhaseFactor(Ar, me%HA_r)
            end if

            ! ! set space domain array to use memory allocted by the library
            ! errorCode = spfft_transform_get_space_domain(me%transform, me%processingUnit, me%realValuesPtr)
            ! if (errorCode /= SPFFT_SUCCESS) error stop

            ! transform backward
            errorCode = spfft_transform_backward(me%transform, me%HA_r, me%processingUnit)
            if (errorCode /= SPFFT_SUCCESS) error stop

            call c_f_pointer(me%realValuesPtr, spaceDomainPtr, [me%nkpts])

            do ik=1,me%nkpts_loc
               ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
               Hk(i,j,ik) = spaceDomainPtr(ik_glob)
               if(i < j) Hk(j,i,ik) = conjg(Hk(i,j,ik))
            end do
         end do
      end do

   end subroutine GetHam
!--------------------------------------------------------------------------------------
   subroutine GetGradHam(me,kdist,GradHk,Ar)
      class(wann_spfft_t) :: me
      type(dist_array1d_t),intent(in) :: kdist
      complex(dp),intent(inout) :: GradHk(:,:,:,:)
      real(dp),intent(in),optional :: Ar(3)
      integer :: i,j,idir,ik,ik_glob
      integer :: errorCode = 0
      complex(C_DOUBLE_COMPLEX), pointer :: spaceDomainPtr(:)

      call assert_shape(GradHk, [me%nwan,me%nwan,3,me%nkpts_loc], "GetGradHam", "GradHk")

      do j=1,me%nwan
         do i=1,j
            call GetGradiant(me%crvec, me%ham_r(:,i,j), me%gradH_R)
            do idir=1,3
               if(present(Ar)) then
                  call me%ApplyPhaseFactor(Ar, me%gradH_R(:,idir))
               end if

               ! transform backward
               errorCode = spfft_transform_backward(me%transform, me%gradH_R(:,idir), me%processingUnit)
               if (errorCode /= SPFFT_SUCCESS) error stop

               call c_f_pointer(me%realValuesPtr, spaceDomainPtr, [me%nkpts])

               do ik=1,me%nkpts_loc
                  ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
                  GradHk(i,j,idir,ik) = spaceDomainPtr(ik_glob)
                  if(i < j) GradHk(j,i,idir,ik) = conjg(GradHk(i,j,idir,ik))
               end do
            end do
         end do
      end do

   end subroutine GetGradHam
!--------------------------------------------------------------------------------------
     subroutine GetDipole(me,kdist,Dk,Ar)
      class(wann_spfft_t) :: me
      type(dist_array1d_t),intent(in) :: kdist
      complex(dp),intent(inout) :: Dk(:,:,:,:)
      real(dp),intent(in),optional :: Ar(3)
      integer :: i,j,idir,ik,ik_glob
      integer :: errorCode = 0
      complex(C_DOUBLE_COMPLEX), pointer :: spaceDomainPtr(:)

      call assert_shape(Dk, [me%nwan,me%nwan,3,me%nkpts_loc], "GetDipole", "Dk")

      do j=1,me%nwan
         do i=1,j
            do idir=1,3
               me%HA_r = me%pos_r(:,i,j,idir)
               if(present(Ar)) then
                  call me%ApplyPhaseFactor(Ar, me%HA_r)
               end if

               ! ! set space domain array to use memory allocted by the library
               ! errorCode = spfft_transform_get_space_domain(me%transform, me%processingUnit, me%realValuesPtr)
               ! if (errorCode /= SPFFT_SUCCESS) error stop

               ! transform backward
               errorCode = spfft_transform_backward(me%transform, me%HA_r, me%processingUnit)
               if (errorCode /= SPFFT_SUCCESS) error stop

               call c_f_pointer(me%realValuesPtr, spaceDomainPtr, [me%nkpts])

               do ik=1,me%nkpts_loc
                  ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
                  Dk(i,j,idir,ik) = spaceDomainPtr(ik_glob)
                  if(i < j) Dk(j,i,idir,ik) = conjg(Dk(i,j,idir,ik))
               end do
            end do
         end do
      end do

   end subroutine GetDipole
!--------------------------------------------------------------------------------------
   pure subroutine ApplyPhaseFactor(me,Ar,OO_R)
      class(wann_spfft_t),intent(in) :: me
      real(dp),intent(in) :: Ar(3)
      complex(dp),target,intent(inout) :: OO_R(:)
      integer :: ir
      real(dp) :: adot,c,s

      select case(me%kdim)
      case(2)
         do ir=1,me%nrpts
            adot = Ar(1) * me%irvec(ir,1) + Ar(2) * me%irvec(ir,2)
            c = cos(DPI * adot)
            s = sin(DPI * adot)
            OO_R(ir) = cmplx(c,-s,kind=dp) * OO_R(ir)
         end do       
      case(3)
         do ir=1,me%nrpts
            adot = Ar(1) * me%irvec(ir,1) + Ar(2) * me%irvec(ir,2) + Ar(3) * me%irvec(ir,3)
            c = cos(DPI * adot)
            s = sin(DPI * adot)
            OO_R(ir) = cmplx(c,-s,kind=dp) * OO_R(ir)
         end do   
      end select

   end subroutine ApplyPhaseFactor
!--------------------------------------------------------------------------------------



!======================================================================================

end module wan_spfft_ham_mpi
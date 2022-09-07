module Marray1d_dist
!======================================================================================
!!  Provides tools for distributing 1D arrays over mpi tasks
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use mpi
   use Mdef,only: dp
   implicit none
!-------------------------------------------------------------------------------------- 
   private
   public :: GetDisplSize1D, dist_array1d_t
!-------------------------------------------------------------------------------------- 
   type dist_array1d_t
      !! This class stores the distribution scheme
      !! We try to split the array into pieces of equal length if possible;
      !! otherwise, we try to distribute as evenly as possible
      integer :: ntasks,N
      integer,dimension(:),allocatable :: N_loc
      integer,dimension(:,:),allocatable :: I_glob
   contains
      procedure,public  :: Init
      procedure,public  :: Clean
      procedure,public  :: Indx_loc2glob
   end type dist_array1d_t
!--------------------------------------------------------------------------------------
   integer,parameter :: master=0,from_master=1,from_worker=2
!-------------------------------------------------------------------------------------- 
contains
!--------------------------------------------------------------------------------------  
   subroutine GetDisplSize1D(N_loc,elem_size,displ,nsize,size_loc)
      !! returns local size and displacement for MPI_gatherv and MPI_allgatherv
      !! we assume an array A(1:elem_size,1:N) is distributed over the last index
      integer,intent(in)    :: N_loc(0:)  !! number of points for each MPI rank
      integer,intent(in)    :: elem_size  !! size of each element
      integer,intent(inout) :: displ(0:)  !! corresponds to displ in MPI_Gatherv
      integer,intent(out)   :: nsize      !! corresponds to sendcount in MPI_Gatherv
      integer,intent(inout) :: size_loc(0:) !! corresponds to recvcounts in MPI_Gatherv
      integer :: ntasks,taskid,iwork,ierr

      call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)

      displ = 0
      do iwork=1,ntasks-1
         displ(iwork) = displ(iwork-1) + elem_size * N_loc(iwork-1)
      end do
      nsize = elem_size * N_loc(taskid)

      do iwork=0,ntasks-1
         size_loc(iwork) = elem_size * N_loc(iwork)
      end do

   end subroutine GetDisplSize1D
!-------------------------------------------------------------------------------------- 
   subroutine Init(me,ntasks,taskid,N,PrintInfo,serial_mode)
      !! Creates layout scheme to distribute N elements over ntasks MPI ranks.
      !! The binning algorithm distributes the elements as evenly as possible.
      class(dist_array1d_t)  :: me
      integer,intent(in)  :: ntasks !! number of MPI ranks
      integer,intent(in)  :: taskid !! MPI rank index
      integer,intent(in)  :: N !! number of elements
      logical,intent(in),optional :: PrintInfo !! if .true., the distribution scheme will be printed
      logical,intent(in),optional :: serial_mode !! if .true. ntasks=1 is assumed
      logical :: info,serial
      integer :: itask
      integer :: ik,binmax
      integer :: nralloc,remainder,buckets
      integer :: bins(0:ntasks-1),offset(0:ntasks-1)

      info = .false.
      if(present(PrintInfo)) info = PrintInfo

      serial = .false.
      if(present(serial_mode)) serial = serial_mode

      me%ntasks = ntasks
      me%N = N

      if(serial) then
         me%ntasks = 1
         allocate(me%I_glob(0:me%ntasks-1,me%N))
         forall(ik=1:me%N) me%I_glob(0,ik) = ik
         allocate(me%N_loc(0:me%ntasks-1))
         me%N_loc(0) = me%N
         return
      end if

      nralloc = 0
      do itask=0,ntasks-1
         remainder = me%N - nralloc
         buckets = ntasks - itask
         bins(itask) = ceiling(remainder/dble(buckets))
         nralloc = nralloc + bins(itask)
      end do

      offset = 0
      do itask=1,ntasks-1
         offset(itask) = sum(bins(0:itask-1))
      end do

      if(taskid == master.and. info) then
         write(output_unit,'(50("-"))')
         write(output_unit,*) 'offset'
         write(output_unit,'(50("-"))')
         do itask=0,ntasks-1
            print*,itask,offset(itask)
         end do
         write(output_unit,'(50("-"))')
      end if

      allocate(me%N_loc(0:me%ntasks-1))
      do itask=0,ntasks-1
         me%N_loc(itask) = bins(itask)
      end do

      if(taskid == master.and. info) then
         write(output_unit,'(50("-"))')
         write(output_unit,*) 'Nk(loc)'
         write(output_unit,'(50("-"))')
         do itask=0,ntasks-1
            print*,itask,me%N_loc(itask)
         end do
         write(output_unit,'(50("-"))')
         write(output_unit,*)
      end if

      binmax = maxval(bins)

      allocate(me%I_glob(0:ntasks-1,binmax))
      do itask=0,ntasks-1
         do ik=1,me%N_loc(itask)
            me%I_glob(itask,ik) = ik + offset(itask)
         end do
      end do

   end subroutine Init
!-------------------------------------------------------------------------------------- 
   subroutine Clean(me)
      class(dist_array1d_t) :: me

      if(allocated(me%N_loc)) deallocate(me%N_loc)
      if(allocated(me%I_glob)) deallocate(me%I_glob)

   end subroutine Clean
!-------------------------------------------------------------------------------------- 
   pure integer function Indx_Loc2Glob(me,taskid,i_loc)
      !! returns the global index corresponding to the local index i_loc 
      class(dist_array1d_t),intent(in) :: me
      integer,intent(in) :: taskid !! MPI rank index
      integer,intent(in) :: i_loc !! local index

      Indx_Loc2Glob = me%I_glob(taskid,i_loc)

   end function Indx_Loc2Glob
!--------------------------------------------------------------------------------------


!======================================================================================
end module Marray1d_Dist

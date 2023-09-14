module wan_overlap
!! Provides tools for working with non-orthogonal Wannier functions.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
   use wan_fourier,only: fourier_R_to_k
   implicit none
   include '../units_inc.f90'
   include '../formats.h'
!--------------------------------------------------------------------------------------
   private
   public :: wann90_ovlp_t
!--------------------------------------------------------------------------------------
   type :: wann90_ovlp_t
      !! Class to represent the overlap matrix of Wannier functions in a similar way
      !! as the Wannier Hamiltonian.
      integer                                    :: num_wann,nrpts
      integer,allocatable,dimension(:)           :: ndegen
      integer,allocatable,dimension(:,:)         :: irvec
      real(dp),dimension(3,3)                    :: real_lattice,recip_lattice
      complex(dp),allocatable,dimension(:,:,:)   :: S_r
   contains
      procedure,public  :: ReadFromW90
      procedure,public  :: SaveToW90
      procedure,public  :: Set
      procedure,public  :: get_Smat
      procedure,public  :: Clean
#if WITHHDF5
      procedure,public  :: ReadFromHDF5
      procedure,public  :: SaveToHDF5
#endif      
   end type wann90_ovlp_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Clean(me)
      class(wann90_ovlp_t) :: me

      if(allocated(me%ndegen)) deallocate(me%ndegen)
      if(allocated(me%irvec)) deallocate(me%irvec)
      if(allocated(me%S_r)) deallocate(me%S_r)

   end subroutine Clean
!--------------------------------------------------------------------------------------
   subroutine Set(me,w90)
      class(wann90_ovlp_t) :: me
      type(wann90_ovlp_t),intent(in) :: w90

      me%num_wann = w90%num_wann
      me%nrpts = w90%nrpts
      me%real_lattice = w90%real_lattice

      if(allocated(me%ndegen)) deallocate(me%ndegen)
      if(allocated(me%irvec)) deallocate(me%irvec)
      if(allocated(me%S_r)) deallocate(me%S_r)

      allocate(me%ndegen(me%nrpts))
      allocate(me%irvec(me%nrpts,3))
      allocate(me%S_r(me%num_wann,me%num_wann,me%nrpts))

      me%ndegen = w90%ndegen
      me%irvec = w90%irvec
      me%S_r = w90%S_r

   end subroutine Set
!--------------------------------------------------------------------------------------   
   subroutine ReadFromW90(me,fname)
      class(wann90_ovlp_t) :: me
      character(len=*),intent(in) :: fname
      integer :: rst,qst,i,j,irpt,ndx1,ndx2
      real(dp) :: a,b
      integer :: file_unit

      open(newunit=file_unit, file=trim(fname), status='old', action='read')

      read(file_unit, *)  ! Date and time
      !
      ! lattice vectors
      !
      read(file_unit, *) me%real_lattice(1, :) !a_1
      read(file_unit, *) me%real_lattice(2, :) !a_2
      read(file_unit, *) me%real_lattice(3, :) !a_3
      
      ! convert to atomic units
      me%real_lattice = me%real_lattice / BohrAngstrom

      read(file_unit, *) me%num_wann
      read(file_unit, *) me%nrpts
      rst=mod(me%nrpts,15)
      qst=int(me%nrpts/15)
      allocate(me%ndegen(me%nrpts),me%irvec(me%nrpts,3))

      ! read WS degeneracies
      do i=1,qst
         read(file_unit,*) (me%ndegen(j+(i-1)*15),j=1,15)
      end do
      if(rst.ne.0) read(file_unit,*) (me%ndegen(j+qst*15),j=1,rst)

      ! read real-space overlap (no spinup-spindw hybridizations assumed)
      allocate(me%S_r(me%num_wann,me%num_wann,me%nrpts))
      do irpt = 1,me%nrpts
         read(file_unit,*) me%irvec(irpt,1),me%irvec(irpt,2),me%irvec(irpt,3)
         do i=1,me%num_wann
            do j=1,me%num_wann
               read(file_unit,*) ndx1, ndx2, a, b
               me%S_r(j,i,irpt) = cmplx(a,b,kind=dp)
            end do
         end do
      end do

   end subroutine ReadFromW90
!--------------------------------------------------------------------------------------  
 subroutine SaveToW90(me,fname)
      !! Save the Wannier Hamiltonian to file, matching the format of Wannier90
      character(len=*),parameter :: fmt_time='(a,"/",a,"/",a," at ",a,":",a,":",a)'
      class(wann90_ovlp_t)  :: me
      character(len=*),intent(in) :: fname
      character(len=8) :: date
      character(len=10) :: time
      integer :: i,j,irpt
      integer :: file_unit
      
      open(newunit=file_unit, file=trim(fname), status='replace')

      call date_and_time(date=date,time=time)
      write(file_unit, fmt_time) date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
      !
      ! lattice vectors
      !
      write(file_unit, *) BohrAngstrom * me%real_lattice(1, :) !a_1
      write(file_unit, *) BohrAngstrom * me%real_lattice(2, :) !a_2
      write(file_unit, *) BohrAngstrom * me%real_lattice(3, :) !a_3
      
      write(file_unit, *) me%num_wann
      write(file_unit, *) me%nrpts
      write(file_unit, '(15I5)') (me%ndegen(i), i=1, me%nrpts)
      !
      ! <0n|S|Rm>
      !
      do irpt = 1, me%nrpts
         write(file_unit, '(/,3I5)') me%irvec(irpt,:)
         do i = 1, me%num_wann
            do j = 1, me%num_wann
               write(file_unit, '(2I5,3x,2(E15.8,1x))') j, i, me%S_r(j, i, irpt)
            end do
         end do
      end do

      close(file_unit)

   end subroutine SaveToW90
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine SaveToHDF5(me,fname)
      use scitools_hdf5_utils
      class(wann90_ovlp_t)  :: me
      character(len=*),intent(in) :: fname
      integer(HID_t) :: file_id
      real(dp),allocatable :: rdata(:,:,:)

      call hdf_open_file(file_id, trim(fname), STATUS='NEW')

      call hdf_write_attribute(file_id,'','num_wann', me%num_wann)
      call hdf_write_attribute(file_id,'','nrpts', me%nrpts)

      call hdf_write_dataset(file_id,'ndegen',me%ndegen)
      call hdf_write_dataset(file_id,'real_lattice',me%real_lattice)
      call hdf_write_dataset(file_id,'irvec',me%irvec)

      allocate(rdata(me%num_wann,me%num_wann,me%nrpts))

      rdata = dble(me%S_r)
      call hdf_write_dataset(file_id,'S_r_real',rdata)
      rdata = aimag(me%S_r)
      call hdf_write_dataset(file_id,'S_r_imag',rdata)

      deallocate(rdata)

      call hdf_close_file(file_id)

   end subroutine SaveToHDF5
#endif
!--------------------------------------------------------------------------------------
#if WITHHDF5
   subroutine ReadFromHDF5(me,fname)
      !! Reads the Wannier Hamiltonian from HDF5 binary format
      use scitools_hdf5_utils
      class(wann90_ovlp_t)  :: me
      character(len=*),intent(in) :: fname
      integer(HID_t) :: file_id
      integer :: ir
      real(dp),allocatable :: rdata(:,:,:)
      real(dp),allocatable :: d_pos_r(:,:,:,:)

      call hdf_open_file(file_id, trim(fname), STATUS='OLD', ACTION='READ')

      call hdf_read_attribute(file_id,'','num_wann', me%num_wann)
      call hdf_read_attribute(file_id,'','nrpts', me%nrpts)

      allocate(me%ndegen(me%nrpts))
      call hdf_read_dataset(file_id,'ndegen',me%ndegen)
      call hdf_read_dataset(file_id,'real_lattice',me%real_lattice)

      allocate(me%irvec(me%nrpts,3))
      call hdf_read_dataset(file_id,'irvec',me%irvec)

      allocate(rdata(me%num_wann,me%num_wann,me%nrpts))
      allocate(me%S_r(me%num_wann,me%num_wann,me%nrpts))

      call hdf_read_dataset(file_id,'S_r_real',rdata)
      me%S_r = rdata
      call hdf_read_dataset(file_id,'S_r_imag',rdata)
      me%S_r = me%S_r + iu * rdata

      deallocate(rdata)

      call hdf_close_file(file_id)

   end subroutine ReadFromHDF5
#endif
!--------------------------------------------------------------------------------------
   function get_Smat(me,kpt) result(Sk)
      class(wann90_ovlp_t) :: me
      real(dp),intent(in) :: kpt(3)
      complex(dp) :: Sk(me%num_wann,me%num_wann)

      call fourier_R_to_k(kpt, me%irvec, me%ndegen, me%S_r, Sk)

   end function get_Smat
!--------------------------------------------------------------------------------------

!======================================================================================
end module wan_overlap

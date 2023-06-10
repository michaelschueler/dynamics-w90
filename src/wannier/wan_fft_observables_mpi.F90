module wan_fft_observables_mpi
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use mpi
   use omp_lib
   use Mdebug
   use scitools_def,only: dp,iu,one,zero
   use scitools_linalg,only: get_large_size,util_matmul
   use scitools_array1d_dist,only: dist_array1d_t
   use wan_hamiltonian,only: wann90_tb_t
   use wan_fft_ham_mpi,only: wann_fft_t
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   private
   public :: Wann_FFT_Observables_dip
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Wann_FFT_Observables_dip(ham,ham_fft,kdist,AF,EF,Rhok,Ekin,Etot,Jcurr,Jgrad,Dip,dipole_current)
      type(wann90_tb_t),intent(in) :: ham
      type(wann_fft_t),intent(in)  :: ham_fft
      type(dist_array1d_t),intent(in) :: kdist
      real(dp),intent(in)          :: AF(3)
      real(dp),intent(in)          :: EF(3)
      complex(dp),intent(in)       :: Rhok(:,:,:)
      real(dp),intent(out)         :: Ekin
      real(dp),intent(out)         :: Etot
      real(dp),intent(out)         :: Jcurr(3)
      real(dp),intent(out)         :: Jgrad(3)
      real(dp),intent(out)         :: Dip(3)
      logical,intent(in),optional  :: dipole_current
      logical :: dipole_current_
      logical :: large_size
      integer :: nbnd,Nk,Nk_loc,ik,idir
      real(dp) :: Ekin_loc,Etot_loc
      real(dp),dimension(3) :: Ared,Jpol,Jgrad_loc,Jpol_loc,Dip_loc
      complex(dp),allocatable,dimension(:,:) :: Hkt,DRhok_dt
      complex(dp),allocatable,dimension(:,:,:) :: Hk
      complex(dp),allocatable,dimension(:,:,:,:) :: grad_Hk,Dk
      integer :: ntasks,taskid,ierr

      print*, "Wann_FFT_Observables_dip"

      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

      dipole_current_ = .false.
      if(present(dipole_current)) dipole_current_ = dipole_current

      nbnd = ham_fft%nwan
      Nk = ham_fft%nkpts
      Nk_loc = kdist%N_loc(taskid)
      print*, "[Wann_FFT_Observables_dip ] Nk ", Nk
      print*, "[Wann_FFT_Observables_dip ] Nk_loc ", Nk_loc

      large_size = get_large_size(nbnd)

      call assert_shape(Rhok, [nbnd,nbnd,Nk_loc], "Wann_FFT_Observables_dip", "Rhok")

      Ared = ham%get_kreduced(AF)

      allocate(Hk(nbnd,nbnd,Nk_loc))
      call ham_fft%GetHam(kdist, Hk, Ar=Ared)

      allocate(grad_Hk(nbnd,nbnd,3,Nk_loc))
      call ham_fft%GetGradHam(kdist, grad_Hk, Ar=Ared)

      allocate(Dk(nbnd,nbnd,3,Nk_loc))
      call ham_fft%GetDipole(kdist, Dk, Ar=Ared)

      Ekin_loc = 0.0_dp
      Etot_loc = 0.0_dp
      Jgrad_loc = 0.0_dp
      Jpol_loc = 0.0_dp
      Dip_loc = 0.0_dp

      allocate(Hkt(nbnd,nbnd))
      allocate(DRhok_dt(nbnd,nbnd))

      do ik=1,Nk_loc
         Hkt = Hk(:,:,ik)
         do idir=1,3
            Hkt(:,:) = Hkt(:,:) - EF(idir) * Dk(:,:,idir,ik)
         end do

         Ekin_loc = Ekin_loc + DTRAB(nbnd, Hk(:,:,ik), Rhok(:,:,ik)) / Nk
         Etot_loc = Etot_loc + DTRAB(nbnd, Hkt, Rhok(:,:,ik)) / Nk

         do idir=1,3
            Jgrad_loc(idir) = Jgrad_loc(idir) + DTRAB(nbnd, grad_Hk(:,:,idir,ik), Rhok(:,:,ik)) / Nk
            Dip_loc(idir) = Dip_loc(idir) + DTRAB(nbnd, Dk(:,:,idir,ik), Rhok(:,:,ik)) / Nk
         end do

         if(dipole_current_) then
            DRhok_dt = zero
            call util_matmul(Hkt,Rhok(:,:,ik),DRhok_dt,alpha=-iu,large_size=large_size)
            call util_matmul(Rhok(:,:,ik),Hkt,DRhok_dt,alpha=iu,beta=one,large_size=large_size)    

            do idir=1,3
               Jpol_loc(idir) = Jpol_loc(idir) + DTRAB(nbnd, Dk(:,:,idir,ik), DRhok_dt) / Nk
            end do

         end if

      end do
      deallocate(Hkt,DRhok_dt)

      call MPI_ALLREDUCE(Ekin_loc,Ekin,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(Etot_loc,Etot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(Jgrad_loc,Jgrad,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(Jpol_loc,Jpol,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(Dip_loc,Dip,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)

      Jcurr = Jgrad + Jpol

      deallocate(Hk)
      deallocate(grad_Hk)
      deallocate(Dk)

   end subroutine Wann_FFT_Observables_dip
!--------------------------------------------------------------------------------------
   pure real(dp) function DTRAB(n,A,B)
   !! Computes \(\mathrm{Tr}[A B]\) for two square matrices \(A,B\).
      integer,intent(in) :: n !! the rank of the square matrices \(A,B\)
      complex(dp),intent(in) :: A(:,:),B(:,:) !! the square matrices \(A,B\)
      integer :: i

      DTRAB = 0.0_dp
      do i=1,n
         DTRAB = DTRAB + dble(sum(A(i,:)*B(:,i)))
      end do

   end function DTRAB

!--------------------------------------------------------------------------------------


!======================================================================================
end module wan_fft_observables_mpi
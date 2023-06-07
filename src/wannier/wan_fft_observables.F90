module wan_fft_observables
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use omp_lib
   use Mdebug
   use scitools_def,only: dp,iu,one,zero
   use scitools_linalg,only: get_large_size,util_matmul
   use wan_hamiltonian,only: wann90_tb_t
   use wan_fft_ham,only: wann_fft_t
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   private
   public :: Wann_FFT_Observables_dip
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Wann_FFT_Observables_dip(ham,ham_fft,AF,EF,Rhok,Ekin,Etot,Jcurr,Jgrad,Dip,dipole_current)
      type(wann90_tb_t),intent(in) :: ham
      type(wann_fft_t),intent(in)  :: ham_fft
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
      integer :: nbnd,Nk,ik,idir
      real(dp),dimension(3) :: Ared,Jpol(3)
      complex(dp),allocatable,dimension(:,:) :: Hkt,DRhok_dt
      complex(dp),allocatable,dimension(:,:,:) :: Hk
      complex(dp),allocatable,dimension(:,:,:,:) :: grad_Hk,Dk

      dipole_current_ = .false.
      if(present(dipole_current)) dipole_current_ = dipole_current

      nbnd = ham_fft%nwan
      Nk = ham_fft%nkpts

      large_size = get_large_size(nbnd)

      call assert_shape(Rhok, [nbnd,nbnd,Nk], "Wann_FFT_Observables_dip", "Rhok")

      print*, "Wann_FFT_Observables_dip"

      Ared = ham%get_kreduced(AF)

      allocate(Hk(nbnd,nbnd,Nk))
      call ham_fft%GetHam_Dressed(Ared, Hk)

      allocate(grad_Hk(nbnd,nbnd,3,Nk))
      call ham_fft%GetGradHam_Dressed(Ared, grad_Hk)

      allocate(Dk(nbnd,nbnd,3,Nk))
      call ham_fft%GetDipole_Dressed(Ared, Dk)

      Ekin = 0.0_dp
      Etot = 0.0_dp
      Jcurr = 0.0_dp
      Jgrad = 0.0_dp
      Jpol = 0.0_dp
      Dip = 0.0_dp

      !$OMP PARALLEL PRIVATE(ik,Hkt,DRhok_dt)
      allocate(Hkt(nbnd,nbnd))
      allocate(DRhok_dt(nbnd,nbnd))
      !$OMP DO REDUCTION(+:Ekin,Etot,Jgrad,Jpol,Dip)
      do ik=1,Nk
         Hkt = Hk(:,:,ik)
         do idir=1,3
            Hkt(:,:) = Hkt(:,:) - EF(idir) * Dk(:,:,idir,ik)
         end do

         Ekin = Ekin + DTRAB(nbnd, Hk(:,:,ik), Rhok(:,:,ik)) / Nk
         Etot = Etot + DTRAB(nbnd, Hkt, Rhok(:,:,ik)) / Nk

         do idir=1,3
            Jgrad(idir) = Jgrad(idir) + DTRAB(nbnd, grad_Hk(:,:,idir,ik), Rhok(:,:,ik)) / Nk
            Dip(idir) = Dip(idir) + DTRAB(nbnd, Dk(:,:,idir,ik), Rhok(:,:,ik)) / Nk
         end do

         if(dipole_current_) then
            DRhok_dt = zero
            call util_matmul(Hkt,Rhok(:,:,ik),DRhok_dt,alpha=-iu,large_size=large_size)
            call util_matmul(Rhok(:,:,ik),Hkt,DRhok_dt,alpha=iu,beta=one,large_size=large_size)    

            do idir=1,3
               Jpol(idir) = Jpol(idir) + DTRAB(nbnd, Dk(:,:,idir,ik), DRhok_dt) / Nk
            end do

         end if

      end do
      !$OMP END DO
      deallocate(Hkt,DRhok_dt)
      !$OMP END PARALLEL

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
end module wan_fft_observables
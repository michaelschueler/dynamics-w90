module wan_equilibrium
!! Contains utilities to construct the equilibrium density matrix.
!======================================================================================
   use Mdebug
   use scitools_def,only: dp, zero, one, iu, nfermi
   use scitools_linalg,only: Eigh,util_rotate,util_rotate_cc,get_large_size
   use scitools_root,only: brent
   use wan_hamiltonian,only: wann90_tb_t
   implicit none
!--------------------------------------------------------------------------------------
   private 
   public :: GetChemicalPotential
   public :: Wann_GenRho_eq, Wann_GenRho_eq_calc
   public :: Wann_GenRhok_eq, Wann_GenRhok_eq_calc
#ifdef MPI
   public :: GetChemicalPotential_mpi
#endif
!--------------------------------------------------------------------------------------
   ! internal wrapp variables
   type param_t
      integer :: Nk
      real(dp) :: filling,Beta
      real(dp),dimension(:,:),pointer :: Ek
   end type param_t

   type(param_t) :: p
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------   
   real(dp) function part_func(mu)
      real(dp),intent(in) :: mu
      real(dp) :: num_part

      num_part = sum(nfermi(p%Beta,p%Ek - mu))/p%Nk
      part_func = num_part - p%filling

   end function part_func
!--------------------------------------------------------------------------------------
   function GetChemicalPotential(Nk,Ek,Beta,filling) result(mu)
      real(dp),parameter :: mu_tol=1.0e-8_dp
      integer,intent(in) :: Nk
      real(dp),target,intent(in) :: Ek(:,:)
      real(dp),intent(in) :: Beta
      real(dp),intent(in) :: filling
      real(dp) :: mu
      real(dp) :: Emin,Emax

      p%Beta = Beta
      p%filling = filling
      p%Ek => Ek
      p%Nk = Nk

      Emin = minval(Ek); Emax = maxval(Ek)
      mu = brent(part_func,Emin,Emax,mu_tol)

   end function GetChemicalPotential
!--------------------------------------------------------------------------------------   
#ifdef MPI
   real(dp) function part_func_mpi(mu)
      use mpi
      real(dp),intent(in) :: mu
      real(dp) :: np_loc,num_part
      integer :: ierr

      np_loc = sum(nfermi(p%Beta,p%Ek - mu))/p%Nk
      call MPI_ALLREDUCE(np_loc,num_part,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
      part_func_mpi = num_part - p%filling

   end function part_func_mpi
#endif
!--------------------------------------------------------------------------------------
#ifdef MPI
   function GetChemicalPotential_mpi(Nk,Ek,Beta,filling) result(mu)
      use mpi
      real(dp),parameter :: mu_tol=1.0e-8_dp
      integer,intent(in) :: Nk
      real(dp),target,intent(in) :: Ek(:,:)
      real(dp),intent(in) :: Beta
      real(dp),intent(in) :: filling
      real(dp) :: mu
      real(dp) :: Emin,Emax,Emin_loc,Emax_loc
      integer :: ierr

      p%Beta = Beta
      p%filling = filling
      p%Ek => Ek
      p%Nk = Nk

      Emin_loc = minval(Ek); Emax_loc = maxval(Ek)
      call MPI_ALLREDUCE(Emin_loc,Emin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(Emax_loc,Emax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD, ierr)

      mu = brent(part_func_mpi,Emin,Emax,mu_tol)

   end function GetChemicalPotential_mpi
#endif
!--------------------------------------------------------------------------------------
   subroutine Wann_GenRhok_eq(w90,Nk,kpts,Mu,Beta,Rhok,band_basis)
   !! Generates the equilibrium density matrix \(\rho_\mathrm{eq}(\mathbf{k})\) for given k-points,
   !! either in band or in Wannier basis. The occupations are determined by the Fermi-Dirac distribution.
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian containing the recip. lattice vectors
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      real(dp),intent(in)          :: Mu !! The chemical potential
      real(dp),intent(in)          :: Beta !! The inverse temperature
      complex(dp),intent(inout)    :: Rhok(:,:,:) !! density matrix \(\rho_\mathrm{eq}(\mathbf{k})\)
      logical,intent(in),optional  :: band_basis !! if `.true.`, the band basis is assumed, otherwise
                                                 !! the density matrix is constructed in Wannier basis
      logical :: bands_=.false.
      logical :: large_size
      integer :: ik,j
      real(dp),dimension(w90%num_wann) :: En
      complex(dp),dimension(w90%num_wann,w90%num_wann) :: Hk,Rhod,Qk

      large_size = get_large_size(w90%num_wann)

      if(present(band_basis)) bands_ = band_basis
      do ik=1,Nk
         rhod = zero
         Hk = w90%get_ham([kpts(ik,1),kpts(ik,2),kpts(ik,3)])
         call eigh(Hk,En,Qk)
         do j=1,w90%num_wann
            rhod(j,j) = nfermi(Beta,En(j)-Mu)
         end do   
         if(bands_) then
            Rhok(:,:,ik) = RhoD
         else
            Rhok(:,:,ik) =  util_rotate_cc(w90%num_wann,Qk,RhoD,large_size=large_size)
         end if
      end do

   end subroutine Wann_GenRhok_eq
!--------------------------------------------------------------------------------------
   subroutine Wann_GenRhok_eq_calc(nbnd,Nk,Hk,Mu,Beta,Rhok,band_basis)
   !! Generates the equilibrium density matrix \(\rho_\mathrm{eq}(\mathbf{k})\)
   !! either in band or in Wannier basis. The occupations are determined by the Fermi-Dirac distribution.
   !! Here we assume that the Hamiltonian \(H_0(\mathbf{k})\) has been precomputed and is given
   !! as input.
      integer,intent(in)           :: nbnd !! The number of bands/orbitals
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      complex(dp),intent(in)       :: Hk(:,:,:) !! Hamiltonian \(H_0(\mathbf{k})\),dimension [nbnd,nbnd,Nk]
      real(dp),intent(in)          :: Mu !! The chemical potential
      real(dp),intent(in)          :: Beta !! The inverse temperature
      complex(dp),intent(inout)    :: Rhok(:,:,:) !! density matrix \(\rho_\mathrm{eq}(\mathbf{k})\)
      logical,intent(in),optional  :: band_basis !! if `.true.`, the band basis is assumed, otherwise
                                                 !! the density matrix is constructed in Wannier basis
      logical :: bands_=.false.
      logical :: large_size
      integer :: ik,j
      real(dp),dimension(nbnd) :: En
      complex(dp),dimension(nbnd,nbnd) :: Rhod,Qk

      large_size = get_large_size(nbnd)

      if(present(band_basis)) bands_ = band_basis
      do ik=1,Nk
         rhod = zero
         call eigh(Hk(:,:,ik),En,Qk)
         do j=1,nbnd
            rhod(j,j) = nfermi(Beta,En(j)-Mu)
         end do   
         if(bands_) then
            Rhok(:,:,ik) = RhoD
         else
            Rhok(:,:,ik) =  util_rotate_cc(nbnd,Qk,RhoD,large_size=large_size)
         end if
      end do

   end subroutine Wann_GenRhok_eq_calc
!--------------------------------------------------------------------------------------
   subroutine Wann_GenRho_eq(w90,Nk,kpt,Mu,Beta,Rho,band_basis)
   !! Generates the equilibrium density matrix \(\rho_\mathrm{eq}(\mathbf{k})\) at a single k-point,
   !! either in band or in Wannier basis. The occupations are determined by the Fermi-Dirac distribution.
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian containing the recip. lattice vectors
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpt(3) !! The k-point where the density matrix is computed
      real(dp),intent(in)          :: Mu !! The chemical potential
      real(dp),intent(in)          :: Beta !! The inverse temperature
      complex(dp),intent(inout)    :: Rho(:,:) !! density matrix \(\rho_\mathrm{eq}(\mathbf{k})\)
      logical,intent(in),optional  :: band_basis !! if `.true.`, the band basis is assumed, otherwise
                                                 !! the density matrix is constructed in Wannier basis
      logical :: bands_=.false.
      logical :: large_size
      integer :: j
      real(dp),dimension(w90%num_wann) :: En
      complex(dp),dimension(w90%num_wann,w90%num_wann) :: Hk,Rhod,Qk

      large_size = get_large_size(w90%num_wann)

      if(present(band_basis)) bands_ = band_basis

      rhod = zero
      Hk = w90%get_ham(kpt)
      call eigh(Hk,En,Qk)
      do j=1,w90%num_wann
         rhod(j,j) = nfermi(Beta,En(j)-Mu)
      end do   
      if(bands_) then
         Rho = RhoD
      else
         Rho = util_rotate_cc(w90%num_wann,Qk,RhoD,large_size=large_size)
      end if

   end subroutine Wann_GenRho_eq
!--------------------------------------------------------------------------------------
   subroutine Wann_GenRho_eq_calc(nbnd,H0,Mu,Beta,Rho,band_basis)
   !! Generates the equilibrium density matrix \(\rho_\mathrm{eq}(\mathbf{k})\) at a single k-point
   !! either in band or in Wannier basis. The occupations are determined by the Fermi-Dirac distribution.
   !! Here we assume that the Hamiltonian \(H_0(\mathbf{k})\) has been precomputed and is given
   !! as input.
      integer,intent(in)           :: nbnd !! The number of bands/orbitals
      complex(dp),intent(in)       :: H0(:,:) !! Hamiltonian \(H_0(\mathbf{k})\),dimension [nbnd,nbnd]
      real(dp),intent(in)          :: Mu !! The chemical potential
      real(dp),intent(in)          :: Beta !! The inverse temperature
      complex(dp),intent(inout)    :: Rho(:,:) !! density matrix \(\rho_\mathrm{eq}(\mathbf{k})\)
      logical,intent(in),optional  :: band_basis !! if `.true.`, the band basis is assumed, otherwise
                                                 !! the density matrix is constructed in Wannier basis
      logical :: bands_=.false.
      logical :: large_size
      integer :: j
      real(dp),dimension(nbnd) :: En
      complex(dp),dimension(nbnd,nbnd) :: Rhod,Qk

      large_size = get_large_size(nbnd)

      if(present(band_basis)) bands_ = band_basis

      rhod = zero
      call eigh(H0,En,Qk)
      do j=1,nbnd
         rhod(j,j) = nfermi(Beta,En(j)-Mu)
      end do   
      if(bands_) then
         Rho = RhoD
      else
         Rho =  util_rotate_cc(nbnd,Qk,RhoD,large_size=large_size)
      end if

   end subroutine Wann_GenRho_eq_calc
!--------------------------------------------------------------------------------------

!======================================================================================
end module wan_equilibrium
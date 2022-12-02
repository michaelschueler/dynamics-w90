module wan_dynamics
!======================================================================================
   use Mdebug
   use scitools_def,only: dp, zero, one, iu, nfermi
   use scitools_linalg,only: Eigh,util_matmul,util_rotate,util_rotate_cc,get_large_size
   use scitools_evol,only: GenU_CF2,GenU_CF4,UnitaryStepFBW
   use scitools_rungekutta,only: ODE_step_rk5
   use wan_hamiltonian,only: wann90_tb_t
   implicit none
!--------------------------------------------------------------------------------------
   ! private
   ! public :: Wann_GenHk, Wann_GenGradkHk, Wann_GenVelok, Wann_GenDipk, Wann_GetHk_dip
   ! public :: Wann_GenRhok_eq, Wann_GenRhok_eq_calc
   ! public :: Wann_Rhok_timestep_velo, Wann_Rhok_timestep_dip
   ! public :: Wann_Current_para_velo, Wann_Current_dia_velo, Wann_Current_Intra_velo 
   ! public :: Wann_Pol_velo, Wann_Pol_dip, Wann_Current_dip
   ! public :: Wann_KineticEn, Wann_TotalEn_velo, Wann_TotalEn_dip
   ! public :: Wann_Current_para_velo_calc,Wann_KineticEn_calc,Wann_TotalEn_velo_calc
   ! public :: Wann_Rhok_timestep_velo_calc, Wann_Current_Intra_velo_calc
   ! public :: Wann_Current_dip_kpt, Wann_Current_velo_kpt
   ! public :: Wann_Pol_dip_calc, Wann_Current_dip_calc
   ! public :: Wann_GetDk_dip, Wann_GetGradHk_dip
   ! public :: Wann_DTRAB_kpts
   ! public :: Wann_timestep_RelaxTime
   ! public :: Wann_SpinCurrent_para_velo, Wann_SpinCurrent_dia_velo
   ! public :: Wann_SpinCurrent_dip, Wann_SpinCurrent_dip_calc
!--------------------------------------------------------------------------------------
   ! interface Wann_Rhok_timestep_dip
   !    module procedure :: Wann_Rhok_timestep_dip_field, Wann_Rhok_timestep_dip_free
   ! end interface Wann_Rhok_timestep_dip

   ! interface Wann_timestep_RelaxTime
   !    module procedure :: Wann_timestep_RelaxTime_dip, Wann_timestep_RelaxTime_velo_calc
   ! end interface Wann_timestep_RelaxTime
!--------------------------------------------------------------------------------------
   integer,parameter :: qc = 1
   real(dp),parameter :: c1=0.5_dp-sqrt(3.0_dp)/6.0_dp
   real(dp),parameter :: c2=0.5_dp+sqrt(3.0_dp)/6.0_dp

   abstract interface
      subroutine vecpot_efield_func(t,A,E)
         import :: dp
         implicit none
         real(dp),intent(in) :: t
         real(dp),intent(out) :: A(3),E(3)
      end subroutine vecpot_efield_func
   end interface
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine Wann_GenHk(w90,Nk,kpts,Hk)
   !! Generates the Wannier Hamiltonian \(H(\mathbf{k})\) for given k-points
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian containing the recip. lattice vectors
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:)  !! List of k-points, dimension [Nk,3]
      complex(dp),intent(inout)    :: Hk(:,:,:) !! Hamiltonian \(H(\mathbf{k})\) in Wannier basis,
                                                !! dimension [nbnd,nbnd,Nk]
      integer :: ik  

      do ik=1,Nk
         Hk(:,:,ik) = w90%get_ham([kpts(ik,1),kpts(ik,2),kpts(ik,3)])   
      end do

   end subroutine Wann_GenHk
!--------------------------------------------------------------------------------------
   subroutine Wann_GenGradkHk(w90,Nk,kpts,grad_Hk)
   !! Generates the gradient Wannier Hamiltonian \(\nabla_{\mathbf{k}} H(\mathbf{k})\) for given k-points
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian containing the recip. lattice vectors
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      complex(dp),intent(inout)    :: grad_Hk(:,:,:,:) !! Gradient of the Hamiltonian 
                                                       !! \(\nabla_{\mathbf{k}} H(\mathbf{k})\) in Wannier basis,
                                                       !! dimension [nbnd,nbnd,Nk,3]
      integer :: ik  

      do ik=1,Nk
         grad_Hk(:,:,:,ik) = w90%get_gradk_ham([kpts(ik,1),kpts(ik,2),kpts(ik,3)])   
      end do

   end subroutine Wann_GenGradkHk
!--------------------------------------------------------------------------------------
   subroutine Wann_GenVelok(w90,Nk,kpts,velok)
   !! Generates the velocity matrix elements \(\mathbf{v}_{\alpha\alpha^\prime}(\mathbf{k})\) 
   !! (band basis) for given k-points
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian containing the recip. lattice vectors
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      complex(dp),intent(inout)    :: velok(:,:,:,:) !! velocity matrix elements \(\mathbf{v}_{\alpha\alpha^\prime}(\mathbf{k})\),
                                                     !! dimension [nbnd,nbnd,Nk,3]
      complex(dp) :: vk(w90%num_wann,w90%num_wann,3)
      integer :: ik  

      do ik=1,Nk
         vk = w90%get_velocity([kpts(ik,1),kpts(ik,2),kpts(ik,3)])
         velok(:,:,ik,1:3) = vk(:,:,1:3)
      end do

   end subroutine Wann_GenVelok
!--------------------------------------------------------------------------------------
   subroutine Wann_GenDipk(w90,Nk,kpts,Dipk)
   !! Generates the dipole matrix elements \(\mathbf{D}_{jj^\prime}(\mathbf{k})\) 
   !! (Wannier basis) for given k-points
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian containing the recip. lattice vectors
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      complex(dp),intent(inout)    :: Dipk(:,:,:,:) !! velocity matrix elements \(\mathbf{D}_{jj^\prime}(\mathbf{k})\),
                                                    !! dimension [nbnd,nbnd,Nk,3]
      integer :: ik  
      complex(dp) :: Dk(w90%num_wann,w90%num_wann,3)

      do ik=1,Nk
         Dk = w90%get_dipole([kpts(ik,1),kpts(ik,2),kpts(ik,3)])
         Dipk(:,:,1:3,ik) = Dk(:,:,1:3)
      end do

   end subroutine Wann_GenDipk
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
   subroutine Wann_timestep_RelaxTime_dip(w90,Nk,kpts,tstp,dt,field,T1,T2,Beta,Mu,Rhok,empirical)
   !! Performs the time step \(\rho(\mathbf{k},t) \rightarrow \rho(\mathbf{k},t + \Delta t) \)
   !! assuming phenomenological dissipative dynamics. The Hamiltonian and density matrix 
   !! are treated in the dipole gauge.
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian containing the recip. lattice vectors
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      integer,intent(in)           :: tstp !! the time step where the density matrix is known
      real(dp),intent(in)          :: dt !! time step size 
      procedure(vecpot_efield_func),pointer :: field !! External field: vector potential + electric field
      real(dp),intent(in)          :: T1 !! The time constant for relaxation to equilibrium
      real(dp),intent(in)          :: T2 !! Time time constant for dephasing of off-diagonal components
      real(dp),intent(in)          :: Beta !! effective temperature for equilibrium state
      real(dp),intent(in)          :: Mu !! effective chemical potential for equilibrium state
      complex(dp),intent(inout)    :: Rhok(:,:,:) !! On input: density matrix at `tstp`; on output:
                                                  !! density matrix at `tstp+1`
      logical,intent(in),optional  :: empirical !! if `.true.` the dipole term is ignore (equivalent to
                                                !! the Peieerls substitution)
      logical :: empirical_
      logical :: large_size
      integer :: nbnd,ik,k_index
      real(dp) :: gm1,gm2,gm12
      real(dp) :: kpt(3)
      complex(dp),allocatable :: Rhok_old(:,:)

      empirical_ = .false.
      if(present(empirical)) empirical_ = empirical

      gm1 = 1.0_dp / T1
      gm2 = 1.0_dp / T2
      gm12 = gm1 - gm2

      nbnd = size(Rhok,dim=1)
      call assert_shape(Rhok,[nbnd,nbnd,Nk],"Wann_timestep_RelaxTime_dip","Rhok")

      large_size = get_large_size(nbnd)

      allocate(Rhok_old(nbnd,nbnd))

      do ik=1,Nk
         k_index = ik
         kpt = kpts(ik,:)
         Rhok_old = Rhok(:,:,ik)
         Rhok(:,:,ik) = ODE_step_RK5(nbnd,tstp,dt,deriv_dipole_gauge,Rhok_old)
      end do

      deallocate(Rhok_old)
      !......................................................
      contains
      !......................................................
      function deriv_dipole_gauge(nst,t,yt) result(dydt)
         integer,intent(in) :: nst
         real(dp),intent(in) :: t 
         complex(dp),intent(in) :: yt(:,:)
         complex(dp) :: dydt(nst,nst)
         integer :: i
         real(dp),dimension(nst) :: epsk,occk
         complex(dp),dimension(nst,nst) :: ht,ht_yt,yt_ht,Uk
         complex(dp),dimension(nst,nst) :: rho_off,rhod,rho_eq,rho12
         real(dp) :: AF(3),EF(3)

         AF = 0.0_dp; EF = 0.0_dp
         if(associated(field)) call field(t,AF,EF)

         ht = Wann_GetHk_dip(w90,AF,EF,kpt,reducedA=.false.,&
            Peierls_only=empirical_)
         call EigH(ht,epsk,Uk)

         occk = nfermi(Beta,epsk - Mu)

         RhoD = zero
         do i=1,nst
            RhoD(i,i) = occk(i)
         end do
         rho_eq = util_rotate_cc(nst,Uk,RhoD,large_size=large_size)

         rho_off = util_rotate(nst,Uk,yt,large_size=large_size)
         do i=1,nst
            rho_off(i,i) = zero
         end do
         rho12 = util_rotate_cc(nst,Uk,rho_off,large_size=large_size)

         call util_matmul(ht,yt,ht_yt,large_size=large_size)
         call util_matmul(yt,ht,yt_ht,large_size=large_size)
         dydt = -iu * (ht_yt - yt_ht)
         dydt = dydt - gm1 * (yt - Rho_eq)
         dydt = dydt + gm12 * rho12

      end function deriv_dipole_gauge
      !....................................................
   end subroutine Wann_timestep_RelaxTime_dip
!--------------------------------------------------------------------------------------
   subroutine Wann_timestep_RelaxTime_velo_calc(nbnd,Nk,Hk,vk,tstp,dt,field,T1,T2,Beta,Mu,Rhok)
   !! Performs the time step \(\rho(\mathbf{k},t) \rightarrow \rho(\mathbf{k},t + \Delta t) \)
   !! assuming phenomenological dissipative dynamics. The Hamiltonian and density matrix 
   !! are treated in the velocity gauge with precomputed Hamiltonian \(H_0(\mathbf{k})\) and velocity matrix
   !! elements `\(v^\mu(\mathbf{k}) \ , \ \mu=x,y,z\)`.
      integer,intent(in)           :: nbnd !! number of bands/orbitals
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      complex(dp),intent(in)       :: Hk(:,:,:) !! Hamiltonian \(H_0(\mathbf{k})\),dimension [nbnd,nbnd,Nk]
      complex(dp),intent(in)       :: vk(:,:,:,:) !! velocity matrix elements \(v^\mu(\mathbf{k})\), 
                                                  !! dimension [nbnd,nbnd,Nk,3]
      integer,intent(in)           :: tstp !! the time step where the density matrix is known
      real(dp),intent(in)          :: dt !! time step size 
      procedure(vecpot_efield_func),pointer :: field !! External field: vector potential + electric field
      real(dp),intent(in)          :: T1 !! The time constant for relaxation to equilibrium
      real(dp),intent(in)          :: T2 !! Time time constant for dephasing of off-diagonal components
      real(dp),intent(in)          :: Beta !! effective temperature for equilibrium state
      real(dp),intent(in)          :: Mu !! effective chemical potential for equilibrium state
      complex(dp),intent(inout)    :: Rhok(:,:,:) !! On input: density matrix at `tstp`; on output:
                                                  !! density matrix at `tstp+1`
      logical :: large_size
      integer :: ik,k_index
      real(dp) :: gm1,gm2,gm12
      complex(dp),allocatable :: Rhok_old(:,:)

      gm1 = 1.0_dp / T1
      gm2 = 1.0_dp / T2
      gm12 = gm1 - gm2

      call assert_shape(Hk,[nbnd,nbnd,Nk],"Wann_timestep_RelaxTime_velo_calc","Hk")
      call assert_shape(vk,[nbnd,nbnd,Nk,3],"Wann_timestep_RelaxTime_velo_calc","vk")
      call assert_shape(Rhok,[nbnd,nbnd,Nk],"Wann_timestep_RelaxTime_velo_calc","Rhok")

      large_size = get_large_size(nbnd)

      allocate(Rhok_old(nbnd,nbnd))

      do ik=1,Nk
         k_index = ik
         Rhok_old = Rhok(:,:,ik)
         Rhok(:,:,ik) = ODE_step_RK5(nbnd,tstp,dt,deriv_velocity_gauge,Rhok_old)
      end do

      deallocate(Rhok_old)
      !......................................................
      contains
      !......................................................
      function deriv_velocity_gauge(nst,t,yt) result(dydt)
         integer,intent(in) :: nst
         real(dp),intent(in) :: t 
         complex(dp),intent(in) :: yt(:,:)
         complex(dp) :: dydt(nst,nst)
         integer :: i,j
         real(dp),dimension(nst) :: epsk,occk
         complex(dp),dimension(nst,nst) :: ht,ht_yt,yt_ht,Uk
         complex(dp),dimension(nst,nst) :: rho_off,rhod,rho_eq,rho12
         real(dp) :: AF(3),EF(3)
         integer :: idir

         AF = 0.0_dp; EF = 0.0_dp
         if(associated(field)) call field(t,AF,EF)

         ht = Hk(:,:,k_index)
         do idir=1,3
            ht = ht - qc * AF(idir) * vk(:,:,k_index,idir)
         end do
         call EigH(ht,epsk,Uk)

         occk = nfermi(Beta,epsk - Mu)

         RhoD = zero
         do i=1,nst
            RhoD(i,i) = occk(i)
         end do
         rho_eq = util_rotate_cc(nst,Uk,RhoD,large_size=large_size)

         rho_off = util_rotate(nst,Uk,yt,large_size=large_size)
         do i=1,nst
            rho_off(i,i) = zero
         end do
         rho12 = util_rotate_cc(nst,Uk,rho_off,large_size=large_size)

         call util_matmul(ht,yt,ht_yt,large_size=large_size)
         call util_matmul(yt,ht,yt_ht,large_size=large_size)
         dydt = -iu * (ht_yt - yt_ht)
         dydt = dydt - gm1 * (yt - Rho_eq)
         dydt = dydt + gm12 * rho12

      end function deriv_velocity_gauge
      !....................................................

   end subroutine Wann_timestep_RelaxTime_velo_calc
!--------------------------------------------------------------------------------------
   subroutine Wann_Rhok_timestep_velo(w90,Nk,kpts,tstp,dt,field,Rhok)
      type(wann90_tb_t),intent(in) :: w90
      integer,intent(in)           :: Nk
      real(dp),intent(in)          :: kpts(:,:)
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      procedure(vecpot_efield_func),pointer :: field
      complex(dp),intent(inout)    :: Rhok(:,:,:)
      logical :: large_size
      integer :: ik
      real(dp) :: EF_1(3),AF_1(3),EF_2(3),AF_2(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann) :: H1,H2,Rho_old,Udt
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: vk
      
      large_size = get_large_size(w90%num_wann)

      AF_1 = 0.0_dp; AF_2 = 0.0_dp
      EF_1 = 0.0_dp; EF_2 = 0.0_dp
      if(associated(field)) then
         call field((tstp + c1)*dt,AF_1,EF_1)
         call field((tstp + c2)*dt,AF_2,EF_2)
      end if

      do ik=1,Nk
         Rho_old = Rhok(:,:,ik)
         H1 = w90%get_ham_diag([kpts(ik,1),kpts(ik,2),kpts(ik,3)])
         H2 = H1
         vk = w90%get_velocity([kpts(ik,1),kpts(ik,2),kpts(ik,3)])
         H1 = H1 - qc*(AF_1(1)*vk(:,:,1) + AF_1(2)*vk(:,:,2) + AF_1(3)*vk(:,:,3))
         H2 = H2 - qc*(AF_2(1)*vk(:,:,1) + AF_2(2)*vk(:,:,2) + AF_2(3)*vk(:,:,3))
         call GenU_CF4(dt,H1,H2,Udt)
         call UnitaryStepFBW(w90%num_wann,Udt,Rho_old,Rhok(:,:,ik),large_size=large_size)
      end do

   end subroutine Wann_Rhok_timestep_velo
!--------------------------------------------------------------------------------------
   subroutine Wann_Rhok_timestep_velo_calc(nbnd,Nk,Hk,vk,tstp,dt,field,Rhok)
      integer,intent(in)           :: nbnd,Nk
      complex(dp),intent(in)       :: Hk(:,:,:)
      complex(dp),intent(in)       :: vk(:,:,:,:)
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      procedure(vecpot_efield_func),pointer :: field
      complex(dp),intent(inout)    :: Rhok(:,:,:)
      logical :: large_size
      integer :: ik
      real(dp) :: EF_1(3),AF_1(3),EF_2(3),AF_2(3)
      complex(dp),dimension(nbnd,nbnd) :: H1,H2,Rho_old,Udt
 
      large_size = get_large_size(nbnd)

      AF_1 = 0.0_dp; AF_2 = 0.0_dp
      EF_1 = 0.0_dp; EF_2 = 0.0_dp
      if(associated(field)) then
         call field((tstp + c1)*dt,AF_1,EF_1)
         call field((tstp + c2)*dt,AF_2,EF_2)
      end if

      do ik=1,Nk
         Rho_old = Rhok(:,:,ik)
         H1 = Hk(:,:,ik) - qc*(AF_1(1)*vk(:,:,ik,1) + AF_1(2)*vk(:,:,ik,2) + AF_1(3)*vk(:,:,ik,3))
         H2 = Hk(:,:,ik) - qc*(AF_2(1)*vk(:,:,ik,1) + AF_2(2)*vk(:,:,ik,2) + AF_1(3)*vk(:,:,ik,3))
         call GenU_CF4(dt,H1,H2,Udt)
         call UnitaryStepFBW(nbnd,Udt,Rho_old,Rhok(:,:,ik),large_size=large_size)
      end do

   end subroutine Wann_Rhok_timestep_velo_calc
!--------------------------------------------------------------------------------------
   function Wann_GetHk_dip(w90,Avec,Efield,kpt,reducedA,Peierls_only) result(Hk)
      type(wann90_tb_t),intent(in) :: w90
      real(dp),intent(in)          :: Avec(3),Efield(3)
      real(dp),intent(in)          :: kpt(3)
      logical,intent(in),optional  :: reducedA
      logical,intent(in),optional  :: Peierls_only
      logical :: reducedA_,peierls_
      real(dp) :: Ared(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann) :: Hk
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: Dk

      reducedA_ = .false.
      if(present(reducedA)) reducedA_ = reducedA

      peierls_ = .false.
      if(present(Peierls_only)) peierls_ = Peierls_only

      if(.not.reducedA) then
         Ared = Cart_to_red(w90,Avec)
      else
         Ared = Avec
      end if

      Hk = w90%get_ham([kpt(1)-Ared(1),kpt(2)-Ared(2),kpt(3)-Ared(3)]) 
      if(peierls_) return
      Dk = w90%get_dipole([kpt(1)-Ared(1),kpt(2)-Ared(2),kpt(3)-Ared(3)])
      Hk = Hk - qc*(Dk(:,:,1) * EField(1) + Dk(:,:,2) * EField(2) + Dk(:,:,3) * EField(3))

   end function Wann_GetHk_dip
!--------------------------------------------------------------------------------------
   function Wann_GetDk_dip(w90,Avec,Efield,kpt,reducedA) result(Dk)
      type(wann90_tb_t),intent(in) :: w90
      real(dp),intent(in)          :: Avec(3),Efield(3)
      real(dp),intent(in)          :: kpt(3)
      logical,intent(in),optional  :: reducedA
      logical :: reducedA_
      real(dp) :: Ared(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: Dk

      reducedA_ = .false.
      if(present(reducedA)) reducedA_ = reducedA

      if(.not.reducedA) then
         Ared = Cart_to_red(w90,Avec)
      else
         Ared = Avec
      end if

      Dk = w90%get_dipole([kpt(1)-Ared(1),kpt(2)-Ared(2),kpt(3)-Ared(3)])

   end function Wann_GetDk_dip
!--------------------------------------------------------------------------------------
   function Wann_GetGradHk_dip(w90,Avec,Efield,kpt,reducedA) result(grad_hk)
      type(wann90_tb_t),intent(in) :: w90
      real(dp),intent(in)          :: Avec(3),Efield(3)
      real(dp),intent(in)          :: kpt(3)
      logical,intent(in),optional  :: reducedA
      logical :: reducedA_
      real(dp) :: Ared(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: grad_hk

      reducedA_ = .false.
      if(present(reducedA)) reducedA_ = reducedA

      if(.not.reducedA) then
         Ared = Cart_to_red(w90,Avec)
      else
         Ared = Avec
      end if

      grad_hk = w90%get_gradk_ham([kpt(1)-Ared(1),kpt(2)-Ared(2),kpt(3)-Ared(3)])

   end function Wann_GetGradHk_dip
!--------------------------------------------------------------------------------------
   subroutine Wann_Rhok_timestep_dip_field(w90,Nk,kpts,tstp,dt,field,Rhok,Peierls_only)
      type(wann90_tb_t),intent(in) :: w90
      integer,intent(in)           :: Nk
      real(dp),intent(in)          :: kpts(:,:)
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      procedure(vecpot_efield_func),pointer :: field
      complex(dp),intent(inout)    :: Rhok(:,:,:)
      logical,intent(in),optional  :: Peierls_only
      logical :: peierls_
      logical :: large_size
      integer :: ik
      real(dp) :: EF_1(3),AF_1(3),EF_2(3),AF_2(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann) :: H1,H2,Rho_old,Udt
      
      large_size = get_large_size(w90%num_wann)

      peierls_ = .false.
      if(present(Peierls_only)) peierls_ = Peierls_only

      AF_1 = 0.0_dp; AF_2 = 0.0_dp
      EF_1 = 0.0_dp; EF_2 = 0.0_dp
      call field((tstp + c1)*dt,AF_1,EF_1)
      call field((tstp + c2)*dt,AF_2,EF_2)

      do ik=1,Nk
         Rho_old = Rhok(:,:,ik)
         H1 = Wann_GetHk_dip(w90,AF_1,EF_1,kpts(ik,:),reducedA=.false.,Peierls_only=peierls_)
         H2 = Wann_GetHk_dip(w90,AF_2,EF_2,kpts(ik,:),reducedA=.false.,Peierls_only=peierls_)
         call GenU_CF4(dt,H1,H2,Udt)
         call UnitaryStepFBW(w90%num_wann,Udt,Rho_old,Rhok(:,:,ik),large_size=large_size)
      end do

   end subroutine Wann_Rhok_timestep_dip_field
!--------------------------------------------------------------------------------------
   subroutine Wann_Rhok_timestep_dip_free(w90,Nk,kpts,tstp,dt,AF,Rhok)
      type(wann90_tb_t),intent(in) :: w90
      integer,intent(in)           :: Nk
      real(dp),intent(in)          :: kpts(:,:)
      integer,intent(in)           :: tstp
      real(dp),intent(in)          :: dt
      real(dp),intent(in)          :: AF(3)
      complex(dp),intent(inout)    :: Rhok(:,:,:)
      logical :: large_size
      integer :: ik
      complex(dp),dimension(w90%num_wann,w90%num_wann) :: H1,Rho_old,Udt
      
      large_size = get_large_size(w90%num_wann)

      do ik=1,Nk
         Rho_old = Rhok(:,:,ik)
         H1 = Wann_GetHk_dip(w90,AF,[0.0_dp,0.0_dp,0.0_dp],kpts(ik,:),reducedA=.false.)
         call GenU_CF2(dt,H1,Udt)
         call UnitaryStepFBW(w90%num_wann,Udt,Rho_old,Rhok(:,:,ik),large_size=large_size)
      end do

   end subroutine Wann_Rhok_timestep_dip_free
!--------------------------------------------------------------------------------------
   function Wann_Current_para_velo(w90,Nk,kpts,Rhok) result(Jpara)
      type(wann90_tb_t),intent(in) :: w90
      integer,intent(in)           :: Nk
      real(dp),intent(in)          :: kpts(:,:)
      complex(dp),intent(in)       :: Rhok(:,:,:)
      real(dp)                     :: Jpara(3)
      integer :: ik
      real(dp) :: Jk(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: vk

      Jpara = 0.0_dp
      do ik=1,Nk
         vk = w90%get_velocity([kpts(ik,1),kpts(ik,2),kpts(ik,3)])
         Jk(1) = qc * DTRAB(w90%num_wann,vk(:,:,1),Rhok(:,:,ik))/Nk
         Jk(2) = qc * DTRAB(w90%num_wann,vk(:,:,2),Rhok(:,:,ik))/Nk
         Jk(3) = qc * DTRAB(w90%num_wann,vk(:,:,3),Rhok(:,:,ik))/Nk
         Jpara = Jpara + Jk
      end do

   end function Wann_Current_para_velo
!--------------------------------------------------------------------------------------
   function Wann_Current_para_velo_calc(nbnd,Nk,vk,Rhok) result(Jpara)
      integer,intent(in)           :: nbnd,Nk
      complex(dp),intent(in)       :: vk(:,:,:,:)
      complex(dp),intent(in)       :: Rhok(:,:,:)
      real(dp)                     :: Jpara(3)
      integer :: ik
      real(dp) :: Jk(3)

      Jpara = 0.0_dp

      do ik=1,Nk
         Jk(1) = qc * DTRAB(nbnd,vk(:,:,ik,1),Rhok(:,:,ik))/Nk
         Jk(2) = qc * DTRAB(nbnd,vk(:,:,ik,2),Rhok(:,:,ik))/Nk
         Jk(3) = qc * DTRAB(nbnd,vk(:,:,ik,3),Rhok(:,:,ik))/Nk
         Jpara = Jpara + Jk
      end do

   end function Wann_Current_para_velo_calc
!--------------------------------------------------------------------------------------
   function Wann_SpinCurrent_para_velo(nbnd,Nk,vk,vectk,Rhok) result(Jpara)
      character(len=*),parameter :: proc_name="Wann_SpinCurrent_para_velo"
      integer,intent(in)           :: nbnd
      integer,intent(in)           :: Nk
      complex(dp),intent(in)       :: vk(:,:,:,:)
      complex(dp),intent(in)       :: vectk(:,:,:)
      complex(dp),intent(in)       :: Rhok(:,:,:)
      real(dp)                     :: Jpara(3,3)
      integer :: idir

      call assert_shape(vk, [nbnd, nbnd, Nk, 3], proc_name, "vk")
      call assert_shape(vectk, [nbnd, nbnd, Nk], proc_name, "vectk")
      call assert_shape(Rhok, [nbnd, nbnd, Nk], proc_name, "Rhok")

      do idir=1,3
         Jpara(:,idir) = Wann_DTRAB_kpts_spin(nbnd,Nk,Rhok,vk(:,:,:,idir),rot_mat=vectk)
      end do

   end function Wann_SpinCurrent_para_velo
!--------------------------------------------------------------------------------------
   function Wann_Current_dia_velo(Nk,Avec,Rhok) result(current)
      use scitools_linalg,only: trace
      character(len=*),parameter :: proc_name="Wann_Current_dia_velo"
      integer,intent(in)           :: Nk
      real(dp),intent(in)          :: Avec(3)
      complex(dp),intent(in)       :: Rhok(:,:,:)
      real(dp)                     :: current(3)
      integer :: nbnd,ik
      real(dp) :: Jk(3),npart
 
      nbnd = size(Rhok, dim=1)
      call assert_shape(Rhok, [nbnd, nbnd, Nk], proc_name, "Rhok")

      npart = 0.0_dp
      do ik=1,Nk
         npart = npart + dble(trace(Rhok(:,:,ik)))/Nk
      end do

      current = -qc**2 * npart * Avec

   end function Wann_Current_dia_velo
!--------------------------------------------------------------------------------------
   function Wann_SpinCurrent_dia_velo(nbnd,Nk,Avec,vectk,Rhok) result(Jpara)
      character(len=*),parameter :: proc_name="Wann_SpinCurrent_dia_velo"
      integer,intent(in)           :: nbnd
      integer,intent(in)           :: Nk
      real(dp),intent(in)          :: Avec(3)
      complex(dp),intent(in)       :: vectk(:,:,:)
      complex(dp),intent(in)       :: Rhok(:,:,:)
      real(dp)                     :: Jpara(3,3)
      integer :: idir,i,ik
      complex(dp) :: Aeye(nbnd,nbnd,3)

      call assert_shape(vectk, [nbnd, nbnd, Nk], proc_name, "vectk")
      call assert_shape(Rhok, [nbnd, nbnd, Nk], proc_name, "Rhok")

      Aeye = zero
      do idir=1,3
         do i=1,nbnd
            Aeye(i,i,idir) = - Avec(idir)
         end do
      end do

      Jpara = 0.0_dp
      do ik=1,Nk
         do idir=1,3
            Jpara(:,idir) = Jpara(:,idir) + Wann_DTRAB_spin(nbnd,Rhok(:,:,ik),Aeye(:,:,idir),&
               rot_mat=vectk(:,:,ik))
         end do
      end do

      Jpara = Jpara / Nk

   end function Wann_SpinCurrent_dia_velo
!--------------------------------------------------------------------------------------
   subroutine Wann_Current_velo_kpt(nbnd,Nk,vk,Avec,Rhok,Jk)
      use scitools_linalg,only: trace
      integer,intent(in)     :: nbnd,Nk
      complex(dp),intent(in) :: vk(:,:,:,:)
      real(dp),intent(in)    :: Avec(3)
      complex(dp),intent(in) :: Rhok(:,:,:)
      real(dp),intent(inout) :: Jk(:,:)
      integer  :: ik
      real(dp) :: npart

      do ik=1,Nk
         Jk(1,ik) = qc * DTRAB(nbnd,vk(:,:,ik,1),Rhok(:,:,ik))/Nk
         Jk(2,ik) = qc * DTRAB(nbnd,vk(:,:,ik,2),Rhok(:,:,ik))/Nk
         Jk(3,ik) = qc * DTRAB(nbnd,vk(:,:,ik,3),Rhok(:,:,ik))/Nk
         npart = dble(trace(Rhok(:,:,ik)))/Nk
         Jk(:,ik) = Jk(:,ik) - qc**2 * npart * Avec
      end do

   end subroutine Wann_Current_velo_kpt
!--------------------------------------------------------------------------------------
   function Wann_Current_Intra_velo(w90,Nk,kpts,Rhok) result(Jcurr)
      type(wann90_tb_t),intent(in) :: w90
      integer,intent(in)           :: Nk
      real(dp),intent(in)          :: kpts(:,:)
      complex(dp),intent(in)       :: Rhok(:,:,:)
      real(dp)                     :: Jcurr(3)
      integer :: ik,i
      real(dp) :: Jk(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: vk

      Jcurr = 0.0_dp
      do ik=1,Nk
         vk = w90%get_velocity([kpts(ik,1),kpts(ik,2),kpts(ik,3)])
         Jk = 0.0_dp
         do i=1,w90%num_wann
            Jk(1) = Jk(1) + qc * dble(vk(i,i,1)*Rhok(i,i,ik))/Nk
            Jk(2) = Jk(2) + qc * dble(vk(i,i,2)*Rhok(i,i,ik))/Nk
            Jk(3) = Jk(3) + qc * dble(vk(i,i,3)*Rhok(i,i,ik))/Nk
         end do
         Jcurr = Jcurr + Jk
      end do

   end function Wann_Current_Intra_velo
!--------------------------------------------------------------------------------------
   function Wann_Current_Intra_velo_calc(nbnd,Nk,vk,Rhok) result(Jcurr)
      integer,intent(in)           :: nbnd,Nk
      complex(dp),intent(in)       :: vk(:,:,:,:)
      complex(dp),intent(in)       :: Rhok(:,:,:)
      real(dp)                     :: Jcurr(3)
      integer :: ik,i
      real(dp) :: Jk(3)

      Jcurr = 0.0_dp
      do ik=1,Nk
         Jk = 0.0_dp
         do i=1,nbnd
            Jk(1) = Jk(1) + qc * dble(vk(i,i,ik,1)*Rhok(i,i,ik))/Nk
            Jk(2) = Jk(2) + qc * dble(vk(i,i,ik,2)*Rhok(i,i,ik))/Nk
            Jk(3) = Jk(3) + qc * dble(vk(i,i,ik,3)*Rhok(i,i,ik))/Nk
         end do
         Jcurr = Jcurr + Jk
      end do

   end function Wann_Current_Intra_velo_calc
!--------------------------------------------------------------------------------------
   function Wann_Pol_velo(w90,Nk,kpts,Rhok,band_basis) result(dipole)
   !! Computes the expectation value of the dipole operator with respect to the 
   !! time-dependent density matrix in velocity gauge \(\rho^\mathrm{VG}(\mathbf{k}),t)\):
   !! $$\mathbf{P}(t) = \frac{q}{N} \sum_{\mathbf{k}} \mathrm{Tr}[ 
   !! \mathbf{D}(\mathbf{k}) \rho^\mathrm{VG}(\mathbf{k}),t) ] . $$
   !! The dipole matrix elements \(\mathbf{D}(\mathbf{k})\) are
   !! generate on-the-fly.
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian including the lattice geometry
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      complex(dp),intent(in)       :: Rhok(:,:,:)  !! density matrix at time \(t\), dimension [nbnd,nbnd,Nk]
      logical,intent(in),optional  :: band_basis !! if `.true.` we assume the density matrix is
                                                 !! expressed in the band basis
      real(dp)                     :: dipole(3)
      logical :: band_basis_=.false.
      integer :: ik
      real(dp) :: Dipk(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: Dk

      if(present(band_basis)) band_basis_ = band_basis

      dipole = 0.0_dp
      do ik=1,Nk
         Dk = w90%get_dipole([kpts(ik,1),kpts(ik,2),kpts(ik,3)],band_basis=band_basis_)
         Dipk(1) = qc * DTRAB(w90%num_wann,Dk(:,:,1),Rhok(:,:,ik))/Nk
         Dipk(2) = qc * DTRAB(w90%num_wann,Dk(:,:,2),Rhok(:,:,ik))/Nk
         Dipk(3) = qc * DTRAB(w90%num_wann,Dk(:,:,3),Rhok(:,:,ik))/Nk
         dipole = dipole + Dipk
      end do
      

   end function Wann_Pol_velo
!--------------------------------------------------------------------------------------
   function Wann_Pol_dip(w90,Nk,kpts,Avec,Rhok,rot_mat) result(dipole)
   !! Computes the expectation value of the dipole operator with respect to the 
   !! time-dependent density matrix in dipole gauge \(\rho^\mathrm{DG}(\mathbf{k}),t)\):
   !! $$\mathbf{P}(t) = \frac{q}{N} \sum_{\mathbf{k}} \mathrm{Tr}[ 
   !! \mathbf{D}(\mathbf{k}- q\mathbf{A}(t)) \rho^\mathrm{DG}(\mathbf{k}),t) ] . $$
   !! The dipole matrix elements \(\mathbf{D}(\mathbf{k}- q\mathbf{A}(t))\) are
   !! generate on-the-fly.
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian including the lattice geometry
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      real(dp),intent(in)          :: Avec(3) !! Vector potential \(\mathbf{A}(t)\) at time \(t\)
      complex(dp),intent(in)       :: Rhok(:,:,:) !! density matrix at time \(t\), dimension [nbnd,nbnd,Nk]
      complex(dp),intent(in),optional :: rot_mat(:,:,:) !! rotation matrix to allow a change of basis
      real(dp)                     :: dipole(3)
      integer :: ik,idir
      real(dp) :: Ared(3),kAred(3),Dipk(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: Dk

      Ared = Cart_to_red(w90,Avec)

      dipole = 0.0_dp

      if(present(rot_mat)) then
         do ik=1,Nk
            kAred = kpts(ik,1:3)-Ared(1)
            Dk = w90%get_dipole(kAred,band_basis=.false.)
            do idir=1,3
               Dk(:,:,idir) = util_rotate(w90%num_wann,rot_mat(:,:,ik),Dk(:,:,idir))
               Dipk(idir) = qc * DTRAB(w90%num_wann,Dk(:,:,idir),Rhok(:,:,ik))/Nk
            end do
            dipole = dipole + Dipk
         end do
      else
         do ik=1,Nk
            kAred = kpts(ik,1:3)-Ared(1)
            Dk = w90%get_dipole(kAred,band_basis=.false.)
            Dipk(1) = qc * DTRAB(w90%num_wann,Dk(:,:,1),Rhok(:,:,ik))/Nk
            Dipk(2) = qc * DTRAB(w90%num_wann,Dk(:,:,2),Rhok(:,:,ik))/Nk
            Dipk(3) = qc * DTRAB(w90%num_wann,Dk(:,:,3),Rhok(:,:,ik))/Nk
            dipole = dipole + Dipk
         end do
      end if
      
   end function Wann_Pol_dip
!--------------------------------------------------------------------------------------
   function Wann_Pol_dip_calc(nbnd,Nk,Dk,Rhok) result(dipole)
   !! Computes the expectation value of the dipole operator with respect to the 
   !! time-dependent density matrix \(\rho(\mathbf{k}),t)\):
   !! $$\mathbf{P}(t) = \frac{q}{N} \sum_{\mathbf{k}} \mathrm{Tr}[ 
   !! \mathbf{D}(\mathbf{k}) \rho(\mathbf{k}),t) ] . $$
   !! The dipole matrix elements \(\mathbf{D}(\mathbf{k})\) are input here.
      integer,intent(in)           :: nbnd !! number of bands/orbitals
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      complex(dp),intent(in)       :: Dk(:,:,:,:) !! dipole matrix elements, dimension [nbnd,nbnd,Nk,3]
      complex(dp),intent(in)       :: Rhok(:,:,:) !! density matrix, dimension [nbnd,nbnd,Nk
      real(dp)                     :: dipole(3)
      integer :: ik
      real(dp) :: Dipk(3)
 
      dipole = 0.0_dp

      do ik=1,Nk
         Dipk(1) = qc * DTRAB(nbnd,Dk(:,:,ik,1),Rhok(:,:,ik))/Nk
         Dipk(2) = qc * DTRAB(nbnd,Dk(:,:,ik,2),Rhok(:,:,ik))/Nk
         Dipk(3) = qc * DTRAB(nbnd,Dk(:,:,ik,3),Rhok(:,:,ik))/Nk
         dipole = dipole + Dipk
      end do
    

   end function Wann_Pol_dip_calc
!--------------------------------------------------------------------------------------
   subroutine Wann_Current_dip_kpt(w90,Nk,kpts,Avec,EF,Rhok,Jk)
   !! Computes the time-dependent and momentum-resolved current with respect to the 
   !! time-dependent density matrix in dipole gauge \(\rho^\mathrm{DG}(\mathbf{k}),t\).
   !! The current is given by the sum of
   !! two contributions \(\mathbf{J}(\mathbf{k},t) = \mathbf{J}^{(1)}(\mathbf{k},t) 
   !! + \mathbf{J}^{(2)}(\mathbf{k},t)\).
   !! We split the two contributions as follows:
   !! $$\mathbf{J}^{(1)}(\mathbf{k},t) = \frac{q}{N}
   !! \mathrm{Tr}[\nabla_{\mathbf{k}} H_0(\mathbf{k}-q\mathbf{A}(t)) 
   !! \rho^\mathrm{DG}(\mathbf{k}),t)], $$
   !! $$ \mathbf{J}^{(2)}(\mathbf{k},t) = \frac{q}{N} \mathrm{Tr}\left[ 
   !! \mathbf{D}(\mathbf{k}-q\mathbf{A}(t)) \frac{d}{dt} \rho^\mathrm{DG}(\mathbf{k}),t) \right] .$$
   !! Here the field-field Hamiltonian \(H_0(\mathbf{k})\) and the dipole matrix elements
   !! \(\mathbf{D}(\mathbf{k})\) are computed on-the-fly. 
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian including the lattice geometry
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      real(dp),intent(in)          :: Avec(3) !! Vector potential \(\mathbf{A}(t)\) at time \(t\)
      real(dp),intent(in)          :: EF(3)  !! Electric field \(\mathbf{A}(t)\) at time \(t\)
      complex(dp),intent(in)       :: Rhok(:,:,:) !! density matrix at time \(t\), dimension [nbnd,nbnd,Nk]
      real(dp),intent(inout)       :: Jk(:,:) !! k-resolved current at time \(t\), dimension [3,Nk]
      integer :: ik,idir
      real(dp) :: Ared(3),kAred(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann)   :: DRhok_dt
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: grad_Hk,Dk

      Ared = Cart_to_red(w90,Avec)

      do ik=1,Nk
         kAred = kpts(ik,:) - Ared
         grad_Hk = w90%get_gradk_ham(kAred) 
         Jk(1,ik) = qc * DTRAB(w90%num_wann,grad_Hk(:,:,1),Rhok(:,:,ik))/Nk
         Jk(2,ik) = qc * DTRAB(w90%num_wann,grad_Hk(:,:,2),Rhok(:,:,ik))/Nk
         Jk(3,ik) = qc * DTRAB(w90%num_wann,grad_Hk(:,:,3),Rhok(:,:,ik))/Nk
         DRhok_dt = Get_Drhok_Dt_dip(w90,Avec,EF,Rhok(:,:,ik),[kpts(ik,1),kpts(ik,2),kpts(ik,3)])
         Dk = w90%get_dipole(kAred)
         Jk(1,ik) = Jk(1,ik) + qc * DTRAB(w90%num_wann,Dk(:,:,1),DRhok_dt)/Nk
         Jk(2,ik) = Jk(2,ik) + qc * DTRAB(w90%num_wann,Dk(:,:,2),DRhok_dt)/Nk
         Jk(3,ik) = Jk(3,ik) + qc * DTRAB(w90%num_wann,Dk(:,:,3),DRhok_dt)/Nk
      end do

   end subroutine Wann_Current_dip_kpt
!--------------------------------------------------------------------------------------
   function Wann_Current_dip(w90,Nk,kpts,Avec,EF,Rhok,dipole_current,rot_mat) result(jcurr)
   !! Computes the time-dependent current with respect to the time-dependent density matrix
   !! in dipole gauge \(\rho^\mathrm{DG}(\mathbf{k}),t\). The current is given by the sum of
   !! two contributions \(\mathbf{J}(t) = \mathbf{J}^{(1)}(t) + \mathbf{J}^{(2)}(t)\).
   !! We split the two contributions as follows:
   !! $$\mathbf{J}^{(1)}(t) = \frac{q}{N} \sum_{\mathbf{k}} 
   !! \left(\mathrm{Tr}[\nabla_{\mathbf{k}} H_0(\mathbf{k}-q\mathbf{A}(t)) 
   !! \rho^\mathrm{DG}(\mathbf{k}),t)], $$
   !! $$ \mathbf{J}^{(2)}(t) = \frac{q}{N} \sum_{\mathbf{k}} \left(\mathrm{Tr}[ 
   !! \mathbf{D}(\mathbf{k}-q\mathbf{A}(t)) \frac{d}{dt} \rho^\mathrm{DG}(\mathbf{k}),t) \right] .$$
   !! Here the field-field Hamiltonian \(H_0(\mathbf{k})\) and the dipole matrix elements
   !! \(\mathbf{D}(\mathbf{k})\) are computed on-the-fly. 
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian including the lattice geometry
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      real(dp),intent(in)          :: Avec(3) !! Vector potential \(\mathbf{A}(t)\) at time \(t\)
      real(dp),intent(in)          :: EF(3)  !! Electric field \(\mathbf{A}(t)\) at time \(t\)
      complex(dp),intent(in)       :: Rhok(:,:,:) !! density matrix at time \(t\), dimension [nbnd,nbnd,Nk]
      logical,intent(in),optional  :: dipole_current !! if `.false.`, the contribution \(\mathbf{J}^{(2)}(t)\)
                                                     !! is neglected 
      complex(dp),intent(in),optional :: rot_mat(:,:,:) !! rotation matrix to allow a change of basis
      real(dp)                     :: jcurr(3)
      logical :: dip_curr
      integer :: ik,idir
      real(dp) :: Ared(3),kAred(3),Jk(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann)   :: DRhok_dt
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: grad_Hk,Dk

      dip_curr = .true.
      if(present(dipole_current)) dip_curr = dipole_current

      Ared = Cart_to_red(w90,Avec)

      Jcurr = 0.0_dp

      if(present(rot_mat)) then
         do ik=1,Nk
            kAred = kpts(ik,:) - Ared
            grad_Hk = w90%get_gradk_ham(kAred) 
            do idir=1,3
               grad_Hk(:,:,idir) = util_rotate(w90%num_wann,rot_mat(:,:,ik),grad_Hk(:,:,idir))
               Jk(idir) = qc * DTRAB(w90%num_wann,grad_Hk(:,:,idir),Rhok(:,:,ik))/Nk
            end do
            Jcurr = Jcurr + Jk
            if(dip_curr) then
               Dk = w90%get_dipole(kAred)
               DRhok_dt = Get_Drhok_Dt_dip(w90,Avec,EF,Rhok(:,:,ik),[kpts(ik,1),kpts(ik,2),kpts(ik,3)],&
                  rot_mat=rot_mat(:,:,ik))
               do idir=1,3
                  Dk(:,:,idir) = util_rotate(w90%num_wann,rot_mat(:,:,ik),Dk(:,:,idir))
                  Jk(idir) = qc * DTRAB(w90%num_wann,Dk(:,:,idir),DRhok_dt)/Nk
               end do
               Jcurr = Jcurr + Jk
            end if
         end do
      else
         do ik=1,Nk
            kAred = kpts(ik,:) - Ared
            grad_Hk = w90%get_gradk_ham(kAred) 
            do idir=1,3
               Jk(idir) = qc * DTRAB(w90%num_wann,grad_Hk(:,:,1),Rhok(:,:,ik))/Nk
            end do
            Jcurr = Jcurr + Jk
            if(dip_curr) then
               Dk = w90%get_dipole(kAred)
               DRhok_dt = Get_Drhok_Dt_dip(w90,Avec,EF,Rhok(:,:,ik),[kpts(ik,1),kpts(ik,2),kpts(ik,3)])
               do idir=1,3
                  Jk(idir) = qc * DTRAB(w90%num_wann,Dk(:,:,idir),DRhok_dt)/Nk
               end do
               Jcurr = Jcurr + Jk
            end if
         end do
      end if

   end function Wann_Current_dip
!--------------------------------------------------------------------------------------
   function Wann_SpinCurrent_dip(w90,Nk,kpts,Avec,EF,Rhok,dipole_current) result(Jspin)
   !! Computes the time-dependent current with respect to the time-dependent density matrix
   !! in dipole gauge \(\rho^\mathrm{DG}(\mathbf{k}),t\). The current is given by the sum of
   !! two contributions \(\mathbf{J}(t) = \mathbf{J}^{(1)}(t) + \mathbf{J}^{(2)}(t)\).
   !! We split the two contributions as follows:
   !! $$\mathbf{J}^{(1)}(t) = \frac{q}{N} \sum_{\mathbf{k}} 
   !! \left(\mathrm{Tr}[\nabla_{\mathbf{k}} H_0(\mathbf{k}-q\mathbf{A}(t)) 
   !! \rho^\mathrm{DG}(\mathbf{k}),t)], $$
   !! $$ \mathbf{J}^{(2)}(t) = \frac{q}{N} \sum_{\mathbf{k}} \left(\mathrm{Tr}[ 
   !! \mathbf{D}(\mathbf{k}-q\mathbf{A}(t)) \frac{d}{dt} \rho^\mathrm{DG}(\mathbf{k}),t) \right] .$$
   !! Here the field-field Hamiltonian \(H_0(\mathbf{k})\) and the dipole matrix elements
   !! \(\mathbf{D}(\mathbf{k})\) are computed on-the-fly. 
      character(len=*),parameter :: proc_name="Wann_SpinCurrent_dip"
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian including the lattice geometry
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      real(dp),intent(in)          :: Avec(3) !! Vector potential \(\mathbf{A}(t)\) at time \(t\)
      real(dp),intent(in)          :: EF(3)  !! Electric field \(\mathbf{A}(t)\) at time \(t\)
      complex(dp),intent(in)       :: Rhok(:,:,:) !! density matrix at time \(t\), dimension [nbnd,nbnd,Nk]
      logical,intent(in),optional  :: dipole_current !! if `.false.`, the contribution \(\mathbf{J}^{(2)}(t)\)
                                                     !! is neglected 
      real(dp)                     :: Jspin(3,3)
      logical :: dip_curr
      integer :: ik,idir
      real(dp) :: Ared(3),kAred(3)
      complex(dp),dimension(w90%num_wann,w90%num_wann)   :: DRhok_dt
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: grad_Hk,Dk

      call assert_shape(kpts, [Nk, 3], proc_name, "kpts")
      call assert_shape(Rhok, [w90%num_wann,w90%num_wann,Nk], proc_name, "Rhok")

      dip_curr = .true.
      if(present(dipole_current)) dip_curr = dipole_current

      Ared = Cart_to_red(w90,Avec)

      Jspin = 0.0_dp
      do ik=1,Nk
         kAred = kpts(ik,:) - Ared
         grad_Hk = w90%get_gradk_ham(kAred) 
         do idir=1,3
            Jspin(:,idir) = Jspin(:,idir) + Wann_DTRAB_spin(w90%num_wann,Rhok(:,:,ik),grad_Hk(:,:,idir)) 
         end do
         if(dip_curr) then
            Dk = w90%get_dipole(kAred)
            DRhok_dt = Get_Drhok_Dt_dip(w90,Avec,EF,Rhok(:,:,ik),[kpts(ik,1),kpts(ik,2),kpts(ik,3)])
            do idir=1,3
               Jspin(:,idir) = Jspin(:,idir) + Wann_DTRAB_spin(w90%num_wann,DRhok_dt,Dk(:,:,idir))
            end do
         end if
      end do

      Jspin = Jspin / Nk

   end function Wann_SpinCurrent_dip
!--------------------------------------------------------------------------------------
   function Wann_Current_dip_calc(nbnd,Nk,Hk,grad_Hk,Dk,Rhok,dipole_current) result(Jcurr)
   !! Computes the time-dependent current with respect to the time-dependent density matrix
   !! in dipole gauge \(\rho^\mathrm{DG}(\mathbf{k}),t\). The current is given by
   !! $$\mathbf{J}(t) = \frac{q}{N} \sum_{\mathbf{k}} \left(\mathrm{Tr}[\nabla_{\mathbf{k}} H_0(\mathbf{k}) 
   !! rho^\mathrm{DG}(\mathbf{k}),t)] + mathrm{Tr}[\mathbf{D}(\mathbf{k}) 
   !! rho^\mathrm{DG}(\mathbf{k}),t)] \right).  $$
   !! Here we assume that the external fields have been switched off, and the gradient of the field-free
   !! Hamiltonian \(\nabla_{\mathbf{k}} H_0(\mathbf{k}) \) and the dipole matrix elements
   !! \(\mathbf{D}(\mathbf{k}) \) have been pre-computed and are given as input to this function.
      integer,intent(in)           :: nbnd !! number of bands/orbitals
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      complex(dp),intent(in)       :: Hk(:,:,:) !! Hamiltonian \(H_0(\mathbf{k})\),dimension [nbnd,nbnd,Nk]
      complex(dp),intent(in)       :: grad_Hk(:,:,:,:) !! gradiant of the Hamiltonian
      complex(dp),intent(in)       :: Dk(:,:,:,:)
      complex(dp),intent(in)       :: Rhok(:,:,:)
      logical,intent(in),optional  :: dipole_current
      real(dp)                     :: Jcurr(3)
      logical :: dip_curr
      logical :: large_size
      integer :: ik
      real(dp) :: Jk(3)
      complex(dp),dimension(nbnd,nbnd)   :: DRhok_dt,H_Rho,Rho_H

      dip_curr = .true.
      if(present(dipole_current)) dip_curr = dipole_current

      large_size = get_large_size(nbnd)

      Jcurr = 0.0_dp

      do ik=1,Nk
         Jk(1) = qc * DTRAB(nbnd,grad_Hk(:,:,ik,1),Rhok(:,:,ik))/Nk
         Jk(2) = qc * DTRAB(nbnd,grad_Hk(:,:,ik,2),Rhok(:,:,ik))/Nk
         Jk(3) = qc * DTRAB(nbnd,grad_Hk(:,:,ik,3),Rhok(:,:,ik))/Nk
         Jcurr = Jcurr + Jk
         if(dip_curr) then
            ! DRhok_dt = -iu*(matmul(Hk(:,:,ik),Rhok(:,:,ik)) - matmul(Rhok(:,:,ik),Hk(:,:,ik)))
            call util_matmul(Hk(:,:,ik),Rhok(:,:,ik),H_Rho,large_size=large_size)
            call util_matmul(Rhok(:,:,ik),Hk(:,:,ik),Rho_H,large_size=large_size)
            DRhok_dt = -iu * (H_Rho - Rho_H)
            Jk(1) = qc * DTRAB(nbnd,Dk(:,:,ik,1),DRhok_dt)/Nk
            Jk(2) = qc * DTRAB(nbnd,Dk(:,:,ik,2),DRhok_dt)/Nk
            Jk(3) = qc * DTRAB(nbnd,Dk(:,:,ik,3),DRhok_dt)/Nk
            Jcurr = Jcurr + Jk
         end if
      end do

   end function Wann_Current_dip_calc
!--------------------------------------------------------------------------------------
   function Wann_SpinCurrent_dip_calc(nbnd,Nk,Hk,grad_Hk,Dk,Rhok,dipole_current) result(Jspin)
   !! Computes the time-dependent current with respect to the time-dependent density matrix
   !! in dipole gauge \(\rho^\mathrm{DG}(\mathbf{k}),t\). The current is given by
   !! $$\mathbf{J}(t) = \frac{q}{N} \sum_{\mathbf{k}} \left(\mathrm{Tr}[\nabla_{\mathbf{k}} H_0(\mathbf{k}) 
   !! rho^\mathrm{DG}(\mathbf{k}),t)] + mathrm{Tr}[\mathbf{D}(\mathbf{k}) 
   !! rho^\mathrm{DG}(\mathbf{k}),t)] \right).  $$
   !! Here we assume that the external fields have been switched off, and the gradient of the field-free
   !! Hamiltonian \(\nabla_{\mathbf{k}} H_0(\mathbf{k}) \) and the dipole matrix elements
   !! \(\mathbf{D}(\mathbf{k}) \) have been pre-computed and are given as input to this function.
      character(len=*),parameter :: proc_name="Wann_SpinCurrent_dip_calc"
      integer,intent(in)           :: nbnd !! number of bands/orbitals
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      complex(dp),intent(in)       :: Hk(:,:,:) !! Hamiltonian \(H_0(\mathbf{k})\),dimension [nbnd,nbnd,Nk]
      complex(dp),intent(in)       :: grad_Hk(:,:,:,:) !! gradiant of the Hamiltonian
      complex(dp),intent(in)       :: Dk(:,:,:,:)
      complex(dp),intent(in)       :: Rhok(:,:,:)
      logical,intent(in),optional  :: dipole_current
      real(dp)                     :: Jspin(3,3)
      logical :: dip_curr
      logical :: large_size
      integer :: ik,idir
      complex(dp),dimension(nbnd,nbnd)   :: DRhok_dt,H_Rho,Rho_H

      call assert_shape(Hk, [nbnd,nbnd,Nk], proc_name, "Hk")
      call assert_shape(grad_Hk, [nbnd,nbnd,Nk,3], proc_name, "grad_Hk")
      call assert_shape(Dk, [nbnd,nbnd,Nk,3], proc_name, "Dk")
      call assert_shape(Rhok, [nbnd,nbnd,Nk], proc_name, "Rhok")

      dip_curr = .true.
      if(present(dipole_current)) dip_curr = dipole_current

      large_size = get_large_size(nbnd)

      Jspin = 0.0_dp

      do ik=1,Nk
         do idir=1,3
            Jspin(:,idir) = Jspin(:,idir) + Wann_DTRAB_spin(nbnd,Rhok(:,:,ik),grad_Hk(:,:,ik,idir))
         end do
         if(dip_curr) then
            ! DRhok_dt = -iu*(matmul(Hk(:,:,ik),Rhok(:,:,ik)) - matmul(Rhok(:,:,ik),Hk(:,:,ik)))
            call util_matmul(Hk(:,:,ik),Rhok(:,:,ik),H_Rho,large_size=large_size)
            call util_matmul(Rhok(:,:,ik),Hk(:,:,ik),Rho_H,large_size=large_size)
            DRhok_dt = -iu * (H_Rho - Rho_H)
            do idir=1,3
               Jspin(:,idir) = Jspin(:,idir) + Wann_DTRAB_spin(nbnd,DRhok_dt,Dk(:,:,ik,idir))
            end do
         end if
      end do

      Jspin = Jspin / Nk

   end function Wann_SpinCurrent_dip_calc
!--------------------------------------------------------------------------------------
   function Get_Drhok_Dt_dip(w90,Avec,Efield,Rhok,kpt,rot_mat) result(DRhok_dt)
   !! Computes the time derivative \(\frac{d}{dt}\rho(\mathbf{k},t) = -i [H(\mathbf{k},t), \rho(\mathbf{k},t) ]\)
   !! needed for the calculation of the current in dipole gauge.
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian including the lattice geometry
      real(dp),intent(in)          :: Avec(3) !! vector potential \(\mathbf{A}(t)\) at time \(t\)
      real(dp),intent(in)          :: Efield(3) !! electric field \(\mathbf{E}(t)\) at time \(t\)
      complex(dp),intent(in)       :: Rhok(:,:) !! density matrix at specific k-point `kpt`, dimension [nbnd,nbnd]
      real(dp),intent(in)          :: kpt(3) !! k-point at which the \(\frac{d}{dt}\rho(\mathbf{k},t) \) is
                                             !! evaluated
      complex(dp),intent(in),optional :: rot_mat(:,:) !! rotation matrix for the option to change basis
      complex(dp),dimension(w90%num_wann,w90%num_wann) :: DRhok_dt
      logical :: large_size
      complex(dp),dimension(w90%num_wann,w90%num_wann) :: Hk,H_Rho, Rho_H

      large_size = get_large_size(w90%num_wann)

      Hk = Wann_GetHk_dip(w90,Avec,Efield,kpt,reducedA=.false.)
      if(present(rot_mat)) then
         Hk = util_rotate(w90%num_wann,rot_mat,Hk,large_size=large_size)
      end if

      call util_matmul(Hk,Rhok,H_Rho,large_size=large_size)
      call util_matmul(Rhok,Hk,Rho_H,large_size=large_size)
      DRhok_dt = -iu * (H_Rho - Rho_H)

      ! DRhok_dt = -iu*(matmul(Hk,Rhok) - matmul(Rhok,Hk))
      
   end function Get_Drhok_Dt_dip
!--------------------------------------------------------------------------------------
   function Wann_TotalEn_velo(w90,Nk,kpts,Avec,Rhok) result(Etot)
   !! Computes the total energy with respect to the time-dependent Hamiltonian
   !! in velocity gauge (VG):
   !! $$E = \frac{1}{N}\sum_{\mathbf{k}} \mathrm{Tr}[H^\mathrm{VG}(\mathbf{k},t) \rho^\mathrm{VG}(\mathbf{k},t)] ,$$
   !! where \(H^\mathrm{VG}(\mathbf{k},t) = H_0(\mathbf{k}) - q \mathbf{A}(t)\cdot\mathbf{v}(\mathbf{k}) +
   !!  \frac{q^2 \mathbf{A}(t)^2}{2}\). 
   !! Here the field-field Hamiltonian \(H_0(\mathbf{k})\) and the velocity matrix elements
   !! \(\mathbf{v}(\mathbf{k})\) are computed on-the-fly. 
      use scitools_linalg,only: trace
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian including the lattice geometry
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      real(dp),intent(in)          :: Avec(3) !! vector potential \(\mathbf{A}(t)\) at time \(t\)
      complex(dp),intent(in)       :: Rhok(:,:,:) !! density matrix, dimension [nbnd,nbnd,Nk]
      real(dp)                     :: Etot
      integer :: ik
      real(dp) :: Ek
      complex(dp),dimension(w90%num_wann,w90%num_wann)   :: Hk
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: vk

      Etot = 0.0_dp
      do ik=1,Nk
         Hk = w90%get_ham_diag([kpts(ik,1),kpts(ik,2),kpts(ik,3)])
         vk = w90%get_velocity([kpts(ik,1),kpts(ik,2),kpts(ik,3)])
         Hk = Hk - qc*(Avec(1)*vk(:,:,1) + Avec(2)*vk(:,:,2) + Avec(3)*vk(:,:,3))
         Ek = DTRAB(w90%num_wann,Hk,Rhok(:,:,ik))
         Etot = Etot + Ek/Nk
         Etot = Etot + 0.5_dp*qc**2*norm2(Avec)**2*trace(Rhok(:,:,ik))/Nk
      end do

   end function Wann_TotalEn_velo
!--------------------------------------------------------------------------------------
   function Wann_TotalEn_velo_calc(nbnd,Nk,Hk,vk,Avec,Rhok) result(Etot)
   !! Computes the total energy with respect to the time-dependent Hamiltonian
   !! in velocity gauge (VG):
   !! $$E = \frac{1}{N}\sum_{\mathbf{k}} \mathrm{Tr}[H^\mathrm{VG}(\mathbf{k},t) \rho^\mathrm{VG}(\mathbf{k},t)] ,$$
   !! where \(H^\mathrm{VG}(\mathbf{k},t) = H_0(\mathbf{k}) - q \mathbf{A}(t)\cdot\mathbf{v}(\mathbf{k}) +
   !!  \frac{q^2 \mathbf{A}(t)^2}{2}\).
   !! Here we assume that the field-field Hamiltonian \(H_0(\mathbf{k})\) and the velocity matrix elements
   !! \(\mathbf{v}(\mathbf{k})\) have been precomputed are given as input to this function. 
      use scitools_linalg,only: trace
      integer,intent(in)           :: nbnd !! number of bands/orbitals
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      complex(dp),intent(in)       :: Hk(:,:,:) !! Hamiltonian \(H_0(\mathbf{k})\),dimension [nbnd,nbnd,Nk]
      complex(dp),intent(in)       :: vk(:,:,:,:) !! velocity matrix elements \(\mathbf{v}(\mathbf{k})\),
                                                  !! dimension [nbnd,nbnd,Nk,3]
      real(dp),intent(in)          :: Avec(3) !! vector potential \(\mathbf{A}(t)\) at time \(t\)
      complex(dp),intent(in)       :: Rhok(:,:,:) !! density matrix, dimension [nbnd,nbnd,Nk]
      real(dp)                     :: Etot
      integer :: ik
      real(dp) :: Ek

      Etot = 0.0_dp
      do ik=1,Nk
         Ek = DTRAB(nbnd,Hk(:,:,ik),Rhok(:,:,ik))
         Etot = Etot + Ek/Nk
         Etot = Etot + 0.5_dp*qc**2*norm2(Avec)**2*trace(Rhok(:,:,ik))/Nk
      end do

   end function Wann_TotalEn_velo_calc
!--------------------------------------------------------------------------------------
   function Wann_KineticEn(w90,Nk,kpts,Rhok,band_basis) result(Ekin)
   !! Computes the kinetic (field-free) energy 
   !! $$E = \frac{1}{N}\sum_{\mathbf{k}} \mathrm{Tr}[H_0(\mathbf{k}) \rho(\mathbf{k})]. $$
   !! Here the Hamiltonian \(H_0(\mathbf{k})\) is constructed on-the-fly, either in Wannier or 
   !! band basis.
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian including the lattice geometry
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      complex(dp),intent(in)       :: Rhok(:,:,:) !! density matrix, dimension [nbnd,nbnd,Nk]
      logical,intent(in),optional  :: band_basis !! if `.true.`, the band basis (diagonal Hamiltonian)
                                                 !! is assumed
      real(dp)                     :: Ekin
      logical :: bands_ = .false.
      integer :: ik
      real(dp) :: Ek
      complex(dp),dimension(w90%num_wann,w90%num_wann)   :: Hk

      if(present(band_basis)) bands_ = band_basis

      Ekin = 0.0_dp

      if(bands_) then
         do ik=1,Nk
            Hk = w90%get_ham_diag([kpts(ik,1),kpts(ik,2),kpts(ik,3)]) 
            Ek = DTRAB(w90%num_wann,Hk,Rhok(:,:,ik))
            Ekin = Ekin + Ek/Nk
         end do
      else
         do ik=1,Nk
            Hk = w90%get_ham([kpts(ik,1),kpts(ik,2),kpts(ik,3)]) 
            Ek = DTRAB(w90%num_wann,Hk,Rhok(:,:,ik))
            Ekin = Ekin + Ek/Nk
         end do
      end if

   end function Wann_KineticEn
!--------------------------------------------------------------------------------------
   function Wann_KineticEn_calc(nbnd,Nk,Hk,Rhok) result(Ekin)
   !! Computes the kinetic (field-free) energy 
   !! $$E = \frac{1}{N}\sum_{\mathbf{k}} \mathrm{Tr}[H_0(\mathbf{k}) \rho(\mathbf{k})]. $$
   !! Here we use the pre-computed Hamiltonian \(H_0(\mathbf{k})\) for fast evaluation.
      integer,intent(in)           :: nbnd !! number of bands/orbitals
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      complex(dp),intent(in)       :: Hk(:,:,:) !! Hamiltonian \(H_0(\mathbf{k})\),dimension [nbnd,nbnd,Nk]
      complex(dp),intent(in)       :: Rhok(:,:,:) !! density matrix, dimension [nbnd,nbnd,Nk]
      real(dp)                     :: Ekin
      integer :: ik
      real(dp) :: Ek

      Ekin = 0.0_dp
      do ik=1,Nk
         Ek = DTRAB(nbnd,Hk(:,:,ik),Rhok(:,:,ik))
         Ekin = Ekin + Ek/Nk
      end do

   end function Wann_KineticEn_calc
!--------------------------------------------------------------------------------------
   function Wann_TotalEn_dip(w90,Nk,kpts,Avec,Efield,Rhok) result(Etot)
   !! Computes the total energy with respect to the time-dependent Hamiltonian
   !! in dipole gauge (DG):
   !! $$E = \frac{1}{N}\sum_{\mathbf{k}} \mathrm{Tr}[H^\mathrm{DG}(\mathbf{k},t) \rho^\mathrm{DG}(\mathbf{k},t)] $$
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian including the lattice geometry
      integer,intent(in)           :: Nk !! The number of k-points / supercells
      real(dp),intent(in)          :: kpts(:,:) !! List of k-points, dimension [Nk,3]
      real(dp),intent(in)          :: Avec(3) !! Vector potential \(\mathbf{A}(t)\) at time \(t\)
      real(dp),intent(in)          :: Efield(3)  !! Electric field \(\mathbf{A}(t)\) at time \(t\)
      complex(dp),intent(in)       :: Rhok(:,:,:) !! density matrix at time \(t\), dimension [nbnd,nbnd,Nk]
      real(dp)                     :: Etot
      integer :: ik
      real(dp) :: Ared(3),kAred(3),Ek
      complex(dp),dimension(w90%num_wann,w90%num_wann)   :: Hk
      complex(dp),dimension(w90%num_wann,w90%num_wann,3) :: Dk

      Ared = Cart_to_red(w90,Avec)

      Etot = 0.0_dp
      do ik=1,Nk
         kAred = kpts(ik,:) - Ared
         Hk = w90%get_ham(kAred) 
         Dk = w90%get_dipole(kAred)
         Hk = Hk - qc*(Dk(:,:,1) * Efield(1) + Dk(:,:,2) * Efield(2) + Dk(:,:,3) * EField(3))
         Ek = DTRAB(w90%num_wann,Hk,Rhok(:,:,ik))
         Etot = Etot + Ek/Nk
      end do

   end function Wann_TotalEn_dip
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
   real(dp) function Wann_DTRAB_kpts(nbnd,Nk,A,B)
   !! Computes the expectation value 
   !! $$ \langle A B \rangle = \frac{1}{N}\sum_{\mathrm{k}} \mathrm{Tr}[A_\mathbf{k} B_\mathbf{k}] $$
   !! with two operators \(A_\mathbf{k}, B_\mathbf{k}\).
      integer,intent(in) :: nbnd !! number of bands/orbitals = dimension of matrices \(A_\mathbf{k}, B_\mathbf{k}\)
      integer,intent(in) :: Nk !! number of k-points \(N\)
      complex(dp),intent(in) :: A(:,:,:),B(:,:,:) !! operators \(A_\mathbf{k}, B_\mathbf{k}\), 
                                                  !! dimension [nbnd,nbnd,Nk]
      integer :: i,ik
      real(dp) :: val

      !$OMP PARALLEL DO REDUCTION(+:val) COLLAPSE(2)
      do ik=1,Nk
         do i=1,nbnd
            val = val + dble(sum(A(i,:,ik)*B(:,i,ik))) / Nk
         end do
      end do

      Wann_DTRAB_kpts = val

   end function Wann_DTRAB_kpts
!--------------------------------------------------------------------------------------
   function Wann_DTRAB_spin(n,A,B,rot_mat) result(val)
   !! Computes \(\mathrm{Tr}[A \{B,\sigma_\nu\}]\) for two square matrices \(A,B\) 
   !! that are understood in spin-orbital space; \(\sigma_\nu\) are the Pauli matrices.
      integer,intent(in) :: n !! the rank of the square matrices \(A,B\)
      complex(dp),intent(in) :: A(:,:),B(:,:) !! the square matrices \(A,B\)
      complex(dp),intent(in),optional :: rot_mat(:,:) !! rotation matrix to change basis
                                                      !! of the spin operators
      real(dp) :: val(3)
      logical :: large_size
      integer :: norb,i,j,nu
      complex(dp),allocatable :: sgm(:,:,:),sgmk(:,:),Bsp(:,:)

      call assert_shape(A, [n,n], "Wann_DTRAB_spin", "A")
      call assert_shape(B, [n,n], "Wann_DTRAB_spin", "B")
      if(present(rot_mat)) then
         call assert_shape(rot_mat, [n,n], "Wann_DTRAB_spin", "rot_mat")
      end if

      large_size = get_large_size(n)

      norb = nint(n / 2.0_dp)
      allocate(sgm(n,n,3))

      call GetSigma(norb,sgm)

      allocate(Bsp(n,n))
      allocate(sgmk(n,n))

      val = 0.0_dp
      do nu=1,3

         if(present(rot_mat)) then
            sgmk = util_rotate(n,rot_mat,sgm(:,:,nu),large_size=large_size)
         else
            sgmk = sgm(:,:,nu)
         end if

         Bsp = zero
         call util_matmul(B,sgmk,Bsp,alpha=one,beta=zero,large_size=large_size)
         call util_matmul(sgmk,B,Bsp,alpha=one,beta=one,large_size=large_size)
         do i=1,n
            val(nu) = val(nu) + dble(sum(A(i,:) * Bsp(:,i)))
         end do
      end do

      deallocate(sgm)
      deallocate(sgmk)
      deallocate(Bsp)

   end function Wann_DTRAB_spin
!--------------------------------------------------------------------------------------
   function Wann_DTRAB_kpts_spin(nbnd,Nk,A,B,rot_mat) result(val)
   !! Computes the expectation value 
   !! $$ \langle A \{B,\sigma_\nu\} \rangle = \frac{1}{N}\sum_{\mathrm{k}} 
   !! \mathrm{Tr}[A_\mathbf{k} \{B_\mathbf{k},\sigma_\nu\}] $$
   !! with two operators \(A_\mathbf{k}, B_\mathbf{k}\) in spin-orbital (or band) space.
   !! \(\sigma_\nu\) are the Pauli matrices.
      integer,intent(in) :: nbnd !! number of bands/orbitals = dimension of matrices \(A_\mathbf{k}, B_\mathbf{k}\)
      integer,intent(in) :: Nk !! number of k-points \(N\)
      complex(dp),intent(in) :: A(:,:,:),B(:,:,:) !! operators \(A_\mathbf{k}, B_\mathbf{k}\), 
                                                  !! dimension [nbnd,nbnd,Nk]
      complex(dp),intent(in),optional :: rot_mat(:,:,:) !! rotation matrix to change basis
                                                      !! of the spin operators
      real(dp) :: val(3)
      logical :: large_size
      integer :: norb,i,j,ik,nu
      complex(dp),allocatable :: sgm(:,:,:),sgmk(:,:),Bsp(:,:)

      call assert_shape(A, [nbnd, nbnd, Nk], "Wann_DTRAB_kpts_spin", "A")
      call assert_shape(B, [nbnd, nbnd, Nk], "Wann_DTRAB_kpts_spin", "B")
      if(present(rot_mat)) then
         call assert_shape(rot_mat, [nbnd, nbnd, Nk], "Wann_DTRAB_kpts_spin", "rot_mat")
      end if

      large_size = get_large_size(nbnd)

      norb = nint(nbnd / 2.0_dp)
      allocate(sgm(nbnd,nbnd,3))

      call GetSigma(norb,sgm)

      val = 0.0_dp


      allocate(Bsp(nbnd,nbnd))
      allocate(sgmk(nbnd,nbnd))

      do ik=1,Nk
         do nu=1,3

            if(present(rot_mat)) then
               sgmk = util_rotate(nbnd,rot_mat(:,:,ik),sgmk,large_size=large_size)
            else
               sgmk = sgm(:,:,nu)
            end if

            Bsp = zero
            call util_matmul(B(:,:,ik),sgmk,Bsp,alpha=one,beta=zero,large_size=large_size)
            call util_matmul(sgmk,B(:,:,ik),Bsp,alpha=one,beta=one,large_size=large_size)
            do i=1,nbnd
               val(nu) = val(nu) + dble(sum(A(i,:,ik) * Bsp(:,i))) / Nk
            end do
         end do
      end do

      deallocate(sgmk)
      deallocate(Bsp)
      deallocate(sgm)

   end function Wann_DTRAB_kpts_spin
!--------------------------------------------------------------------------------------
   function Cart_to_red(w90,kvec) result(kred)
   !! converts a k-vector given in cartisian coordinates to a k-point in reduced coordinates
      type(wann90_tb_t),intent(in) :: w90 !! Wannier Hamiltonian containing the recip. lattice vectors
      real(dp),intent(in)          :: kvec(3) !! k-vector in cartesian coordinates
      real(dp)                     :: kred(3) !! k-point in reducec coordinates

      kred = matmul(w90%recip_reduced(1:3,1:3), kvec)

   end function Cart_to_red
!--------------------------------------------------------------------------------------
   subroutine GetSigma(norb,sgm)
      integer,intent(in) :: norb
      complex(dp),intent(inout) :: sgm(:,:,:)
      integer :: i,j

      sgm = zero

      do i=1,norb
         do j=1,norb
            sgm(2*i-1, 2*j, 1) = one
            sgm(2*i, 2*j-1, 1) = one
            sgm(2*i-1, 2*j, 2) = -iu
            sgm(2*i, 2*j-1, 2) = iu
            sgm(2*i-1, 2*j-1, 3) = one
            sgm(2*i, 2*j, 3) = -one
         end do
      end do


   end subroutine GetSigma
!--------------------------------------------------------------------------------------

!======================================================================================
end module wan_dynamics

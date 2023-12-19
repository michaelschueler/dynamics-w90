module Marpes_calc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero
   use scitools_utils,only: str, linspace, savetxt, stop_error
   use scitools_vector_bsplines,only: cplx_matrix_spline_t
   use wan_latt_kpts,only: kpoints_t
   use wan_hamiltonian,only: wann90_tb_t
   use wan_overlap,only: wann90_ovlp_t
   use wan_compress,only: PruneHoppings
   use wan_slab,only: Wannier_BulkToSlab, Overlap_BulkToSlab
   use wan_orbitals,only: wannier_orbs_t
   use wan_utils,only: Batch_Diagonalize_t
   use pes_scatt_input, only: scatt_input_t
   use pes_radialwf,only: radialwf_t
   use pes_scattwf,only: scattwf_t
   use pes_radialintegral,only: radialinteg_t
   use pes_main, only: PES_Intensity, PES_Slab_Intensity, &
                       PES_Bulk_Intensity, PES_AtomicIntegrals_lambda
   use io_params,only: HamiltonianParams_t, PESParams_t
   use io_hamiltonian,only: ReadHamiltonian, ReadOverlap
   use io_orbitals,only: ReadWannierOrbitals
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: arpes_calc_t
!--------------------------------------------------------------------------------------
   type :: arpes_calc_t
      logical     :: lambda_mode=.false.,slab_mode=.false.,dipole_approx=.true.
      logical     :: bulk_mode = .false.
      logical     :: orthogonal_basis=.true.
      logical     :: scatt_from_input = .false.
      integer     :: nbnd,norb
      integer     :: nlayer=0
      integer     :: bulk_numpoints_kz
      integer     :: gauge
      integer     :: Nepe
      integer     :: lmax,radint_nk,radint_nr
      real(dp)    :: wphot,MuChem,Eshift,lambda_esc,eta_smear, Vinner
      real(dp)    :: qphot(3)=[0.0_dp,0.0_dp,0.0_dp]
      complex(dp) :: polvec(3)
      integer     :: Nk
      integer, allocatable, dimension(:)    :: excluded_layers
      real(dp),allocatable,dimension(:)   :: Epe
      real(dp),allocatable,dimension(:,:) :: kpts,spect
      type(scatt_input_t)  :: scatt_input
      type(wann90_tb_t)    :: ham
      type(wann90_ovlp_t)  :: ovlp
      type(wannier_orbs_t) :: orbs
      type(scattwf_t),allocatable,dimension(:)            :: chis
      type(radialinteg_t),allocatable,dimension(:)        :: radints
      type(cplx_matrix_spline_t),allocatable,dimension(:) :: bessel_integ
   contains
      procedure,public  :: Init
      procedure,public  :: CalcIntegrals
      procedure,private :: CalcIntegrals_radial
      procedure,private :: CalcIntegrals_lambda
      procedure,public  :: CalcPES
      procedure, private :: CalcPES_2D
      procedure, private :: CalcPES_Bulk
      procedure, private :: CalcPES_Slab
      procedure,public  :: WriteSpectrum
      procedure,private :: WriteSpectrum_txt
#ifdef WITHHDF5
      procedure,private :: WriteSpectrum_hdf5
#endif      
   end type arpes_calc_t
!--------------------------------------------------------------------------------------
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   integer,parameter :: gauge_len=0, gauge_mom=1
   integer, parameter :: scattwf_pw = 0, scattwf_coul = 1, scattwf_inp = 2
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine WannOrb_to_RadialWF(orbs,iorb,rwf)
      !! This is a wrapper that converts a the radial part of the Wannier orbitals 
      !! (type wannier_orbs_t) to a radial wave-function (type radialwf_t)
      integer,parameter :: wf_slater=0, wf_grid=1
      type(wannier_orbs_t),intent(in) :: orbs
      integer,intent(in)              :: iorb
      type(radialwf_t),intent(out)    :: rwf

      if(orbs%wf_type == wf_slater) then
         call rwf%InitSlater(orbs%Zorb(iorb),orbs%N_indx(iorb),orbs%L_indx(iorb))
      else
         call rwf%InitGrid(orbs%rs,orbs%Rrad(:,iorb))
      end if

   end subroutine WannOrb_to_RadialWF
!--------------------------------------------------------------------------------------
   subroutine Init(me,par_ham,par_pes,kp)
      class(arpes_calc_t) :: me
      type(HamiltonianParams_t),intent(in) :: par_ham
      type(PESParams_t),intent(in)         :: par_pes
      type(kpoints_t),intent(in)           :: kp
      integer :: ik,iorb,ilay
      real(dp) :: comp_rate
      real(dp) :: kvec(3),kred(3)
      type(wann90_tb_t) :: ham_tmp
      type(wann90_ovlp_t) :: ovlp_tmp

      call ReadHamiltonian(par_ham%file_ham,ham_tmp,file_xyz=par_ham%file_xyz)
      if(par_ham%energy_thresh > 0.0_dp) then
         call PruneHoppings(par_ham%energy_thresh,ham_tmp,me%ham,comp_rate)
         write(output_unit,fmt_info) "compression rate: "//str(nint(100 * comp_rate)) // "%"
      else
         call me%ham%Set(ham_tmp)
      end if
      call ham_tmp%Clean()

      if(len_trim(par_ham%file_ovlp) > 0) then
         call ReadOverlap(par_ham%file_ovlp,ovlp_tmp)
         if(par_ham%ovlp_thresh > 0.0_dp) then
            call PruneHoppings(par_ham%ovlp_thresh,ovlp_tmp,me%ovlp,comp_rate)
            write(output_unit,fmt_info) "overlap compression rate: "//str(nint(100 * comp_rate)) // "%"
         else
            call me%ovlp%Set(ovlp_tmp)
         end if
         call ovlp_tmp%Clean()   
         me%orthogonal_basis = .false.      
      end if

      if(par_ham%slab_mode .and. par_ham%slab_nlayer > 0 .and. .not. par_pes%bulk_mode) then
         write(output_unit,fmt_info) "building slab with "//str(par_ham%slab_nlayer)//" layers"
         call ham_tmp%Set(me%ham)
         call Wannier_BulkToSlab(ham_tmp,par_ham%slab_nlayer,me%ham)
         me%slab_mode = .true.
         me%nlayer = par_ham%slab_nlayer
         call ham_tmp%Clean()

         if(.not.me%orthogonal_basis) then
            call ovlp_tmp%Set(me%ovlp)
            call Overlap_BulkToSlab(ovlp_tmp,par_ham%slab_nlayer,me%ovlp)
            call ovlp_tmp%Clean()
         end if

      end if

      if (par_pes%bulk_mode) then
         write (output_unit, fmt_info) "bulk mode: building infinity periodic slab"
         me%bulk_mode = .true.
         me%bulk_numpoints_kz = par_pes%bulk_numpoints_kz
      end if

      me%nbnd = me%ham%num_wann
      me%MuChem = par_ham%MuChem

      call ReadWannierOrbitals(par_pes%file_orbs,me%orbs)
      me%gauge = par_pes%gauge
      me%Nepe = par_pes%Nepe
      me%wphot = par_pes%wphot
      me%Eshift = par_pes%Eshift
      me%Vinner = par_pes%Vinner
      me%lambda_esc = par_pes%lambda_esc
      me%eta_smear = par_pes%eta_smear    
      me%polvec = par_pes%polvec
      me%radint_nk = par_pes%radint_numpoints_k
      me%radint_nr = par_pes%radint_numpoints_r
      me%lmax = par_pes%expansion_lmax
      me%lambda_mode = par_pes%lambda_orbital_term

      me%dipole_approx = par_pes%dipole_approximation
      if(.not. me%dipole_approx) me%qphot=par_pes%qmom_phot

      me%scatt_from_input = len_trim(par_pes%file_scatt) > 0 &
                            .and. par_pes%scatt_type == scattwf_inp

      if (me%scatt_from_input .and. me%lambda_mode) then
         call stop_error("scattering input not implemented with imaginary part of wave-vector")
      end if

      if (me%scatt_from_input) then
         call me%scatt_input%ReadFromFile(par_pes%file_scatt)
      end if

      if(par_ham%exclude_orbitals) then
         do iorb=1,size(par_ham%orbs_excl, dim=1)
            me%orbs%weight(par_ham%orbs_excl(iorb)) = 0.0_dp
         end do
      end if

      me%norb = me%orbs%norb

      allocate(me%Epe(me%Nepe))
      if(me%Nepe == 1) then
         allocate(me%Epe(1))
         me%Epe(1) = par_pes%Epe_min
      else
         me%Epe = linspace(par_pes%Epe_min, par_pes%Epe_max, me%Nepe)
      end if

      allocate (me%chis(me%norb))
      if (me%scatt_from_input) then
         do iorb = 1, me%norb
            call me%chis(iorb)%Init(scattwf_inp, scatt_input=me%scatt_input, &
                                    scatt_iorb=iorb)
         end do
      else
         do iorb = 1, me%norb
            call me%chis(iorb)%Init(par_pes%scatt_type, me%orbs%Zscatt(iorb))
         end do
      end if


      me%Nk = kp%Nk
      allocate(me%kpts(me%Nk,2))
      if(par_pes%kpts_reduced) then
         do ik=1,me%Nk
            me%kpts(ik,1:2) = me%ham%recip_lattice(1,1:2) * kp%kpts(ik,1) + &
               me%ham%recip_lattice(2,1:2) * kp%kpts(ik,2) 
         end do
      else
         me%kpts(1:me%Nk,1:2) = kp%kpts(1:me%Nk,1:2)
      end if

      select case(par_pes%gauge)
      case(gauge_len)
         write(output_unit,fmt_info) "light-matter coupling: length gauge"
      case(gauge_mom)
         write(output_unit,fmt_info) "light-matter coupling: velocity gauge"         
      end select

      select case (par_pes%scatt_type)
      case (scattwf_pw)
         write (output_unit, fmt_info) "final states: plane waves"
      case (scattwf_coul)
         write (output_unit, fmt_info) "final states: Coulomb waves"
      case (scattwf_inp)
         write (output_unit, fmt_info) "radial integrals & phase shifts from input"
      end select   

   end subroutine Init
!--------------------------------------------------------------------------------------
   subroutine CalcIntegrals(me)
      class(arpes_calc_t) :: me

      if(me%lambda_mode) then
         call me%CalcIntegrals_lambda()
      else
         call me%CalcIntegrals_radial()
      end if

   end subroutine CalcIntegrals
!--------------------------------------------------------------------------------------
   subroutine CalcIntegrals_radial(me)
      class(arpes_calc_t) :: me
      integer :: iorb
      type(radialwf_t) :: rwf
      real(dp) :: kmin,kmax

      if (minval(me%Epe) < 1.0e-8_dp) then
         call stop_error("photoelectron energy minimum is negative")
      end if

      kmin = 0.95_dp * sqrt(2.0_dp * minval(me%Epe))
      kmax = 1.05_dp * sqrt(2.0_dp * maxval(me%Epe))  

      allocate(me%radints(me%nbnd))

      if (me%scatt_from_input) then
         do iorb = 1, me%norb
            call me%radints(iorb)%SetFromInput(me%orbs%L_indx(iorb), me%scatt_input, iorb, &
                                               gauge=me%gauge)
         end do
      else
         do iorb = 1, me%norb
            call WannOrb_to_RadialWF(me%orbs, iorb, rwf)
            call me%radints(iorb)%Init(me%orbs%L_indx(iorb), kmin, kmax, me%chis(iorb), rwf, &
                                       nk=me%radint_nk, gauge=me%gauge)
            call rwf%Clean()
         end do
      end if


   end subroutine CalcIntegrals_radial
!--------------------------------------------------------------------------------------
   subroutine CalcIntegrals_lambda(me)
      class(arpes_calc_t) :: me
      integer :: iorb
      type(radialwf_t) :: rwf
      real(dp) :: kmin,kmax

      if (minval(me%Epe) < 1.0e-8_dp) then
         call stop_error("photoelectron energy minimum is negative")
      end if

      kmin = 0.95_dp * sqrt(2.0_dp * minval(me%Epe))
      kmax = 1.05_dp * sqrt(2.0_dp * maxval(me%Epe))  

      call PES_AtomicIntegrals_lambda(me%orbs,me%chis,me%lambda_esc,me%lmax,kmin,kmax,me%bessel_integ,&
         gauge=me%gauge,Nr=me%radint_nr,Nk=me%radint_nk)

   end subroutine CalcIntegrals_lambda
!--------------------------------------------------------------------------------------
   subroutine CalcPES(me)
      class(arpes_calc_t) :: me

      if (.not. allocated(me%spect)) allocate (me%spect(me%Nepe, me%Nk))

      if (me%bulk_mode) then
         call me%CalcPES_Bulk()
      elseif (me%slab_mode) then
         call me%CalcPES_Slab()
      else
         call me%CalcPES_2D()
      end if

   end subroutine CalcPES
!--------------------------------------------------------------------------------------
   subroutine CalcPES_Bulk(me)
      class(arpes_calc_t) :: me
      integer :: ik, ikz, iepe
      real(dp) :: kpar(2), kqpar(2), kpt(3)
      real(dp), allocatable, dimension(:)        :: kz_red
      real(dp), allocatable, dimension(:, :)      :: epsk
      complex(dp), allocatable, dimension(:, :)   :: Hk, Sk
      complex(dp), allocatable, dimension(:, :, :) :: vectk
      type(Batch_Diagonalize_t) :: batch_diag

      allocate (epsk(me%nbnd, me%bulk_numpoints_kz))
      allocate (vectk(me%nbnd, me%nbnd, me%bulk_numpoints_kz))
      allocate (Hk(me%nbnd, me%nbnd))

      if (.not. me%orthogonal_basis) allocate (Sk(me%nbnd, me%nbnd))

      call batch_diag%Init(me%nbnd)

      allocate (kz_red(me%bulk_numpoints_kz))
      do ikz = 1, me%bulk_numpoints_kz
         kz_red(ikz) = -0.5_dp + (ikz - 1)/dble(me%bulk_numpoints_kz)
      end do

      kpt = 0.0_dp
      do ik = 1, me%Nk
         kpar(1:2) = me%kpts(ik, 1:2)

         if(.not. me%dipole_approx) then
            kqpar(1:2) = kpar(1:2) - me%qphot(1:2)
            kpt(1:2) = me%ham%get_kreduced(kqpar(1:2))
         else
            kpt(1:2) = me%ham%get_kreduced(kpar(1:2))
         end if

         do ikz = 1, me%bulk_numpoints_kz
            kpt(3) = kz_red(ikz)
            Hk = me%ham%get_ham(kpt)
            if (.not. me%orthogonal_basis) Sk = me%ovlp%get_Smat(kpt)

            if (me%orthogonal_basis) then
               call batch_diag%Diagonalize(Hk, epsk(:, ikz), vectk(:, :, ikz))
            else
               call batch_diag%DiagonalizeGen(Hk, Sk, epsk(:, ikz), vectk(:, :, ikz))
            end if
            epsk(:,ikz) = epsk(:,ikz) + me%Eshift
         end do

         if (me%lambda_mode) then
            !$OMP PARALLEL
            !$OMP DO
            do iepe = 1, me%Nepe
               me%spect(iepe, ik) = PES_Bulk_Intensity( &
                                    me%ham, &
                                    me%chis, &
                                    me%lmax, &
                                    me%bessel_integ, &
                                    kpar, &
                                    me%wphot, &
                                    me%polvec, &
                                    me%Epe(iepe), &
                                    epsk, &
                                    vectk, &
                                    me%MuChem, &
                                    me%Vinner, &
                                    me%lambda_esc, &
                                    kz_red, &
                                    me%eta_smear, &
                                    qphot=me%qphot)
            end do
            !$OMP END DO
            !$OMP END PARALLEL
         else
            !$OMP PARALLEL
            !$OMP DO
            do iepe = 1, me%Nepe
               me%spect(iepe, ik) = PES_Bulk_Intensity( &
                                    me%orbs, &
                                    me%ham, &
                                    me%chis, &
                                    me%radints, &
                                    kpar, &
                                    me%wphot, &
                                    me%polvec, &
                                    me%Epe(iepe), &
                                    epsk, &
                                    vectk, &
                                    me%MuChem, &
                                    me%Vinner, &
                                    me%lambda_esc, &
                                    kz_red, &
                                    me%eta_smear, &
                                    me%gauge, &
                                    qphot=me%qphot)
            end do
            !$OMP END DO
            !$OMP END PARALLEL
         end if

      end do

      deallocate (kz_red)
      deallocate (epsk, Hk, vectk)
      if (allocated(Sk)) deallocate (Sk)

      call batch_diag%Clean()

   end subroutine CalcPES_Bulk
!--------------------------------------------------------------------------------------
   subroutine CalcPES_2D(me)
      class(arpes_calc_t) :: me
      integer :: ik, iepe
      real(dp) :: kpar(2), kqpar(2), kpt(3), kqpt(3)
      real(dp), allocatable, dimension(:)        :: epsk
      complex(dp), allocatable, dimension(:, :)   :: Hk, Sk, vectk
      type(Batch_Diagonalize_t) :: batch_diag

      allocate (epsk(me%nbnd), Hk(me%nbnd, me%nbnd), vectk(me%nbnd, me%nbnd))

      if (.not. me%orthogonal_basis) allocate (Sk(me%nbnd, me%nbnd))

      call batch_diag%Init(me%nbnd)

      kpt = 0.0_dp; kqpt = 0.0_dp
      do ik = 1, me%Nk
         kpar(1:2) = me%kpts(ik, 1:2)

         if(.not. me%dipole_approx) then
            kqpar(1:2) = kpar(1:2) - me%qphot(1:2)
            kqpt(1:2) = me%ham%get_kreduced(kqpar(1:2))
            Hk = me%ham%get_ham(kqpt)
            if (.not. me%orthogonal_basis) Sk = me%ovlp%get_Smat(kqpt)
         else
            kpt(1:2) = me%ham%get_kreduced(kpar(1:2))
            Hk = me%ham%get_ham(kpt)
            if (.not. me%orthogonal_basis) Sk = me%ovlp%get_Smat(kpt)
         end if

         if (me%orthogonal_basis) then
            call batch_diag%Diagonalize(Hk, epsk, vectk)
         else
            call batch_diag%DiagonalizeGen(Hk, Sk, epsk, vectk)
         end if
         epsk = epsk + me%Eshift

         if (me%lambda_mode) then
            !$OMP PARALLEL
            !$OMP DO
            do iepe = 1, me%Nepe
               me%spect(iepe, ik) = PES_Intensity( &
                                    me%ham, &
                                    me%chis, &
                                    me%lmax, &
                                    me%bessel_integ, &
                                    kpar, &
                                    me%wphot, &
                                    me%polvec, &
                                    me%Epe(iepe), &
                                    epsk, &
                                    vectk, &
                                    me%MuChem, &
                                    me%Vinner, &
                                    me%lambda_esc, &
                                    me%eta_smear, &
                                    qphot=me%qphot)
            end do
            !$OMP END DO
            !$OMP END PARALLEL
         else
            !$OMP PARALLEL
            !$OMP DO
            do iepe = 1, me%Nepe
               me%spect(iepe, ik) = PES_Intensity( &
                                    me%orbs, &
                                    me%ham, &
                                    me%chis, &
                                    me%radints, &
                                    kpar, &
                                    me%wphot, &
                                    me%polvec, &
                                    me%Epe(iepe), &
                                    epsk, &
                                    vectk, &
                                    me%MuChem, &
                                    me%Vinner, &
                                    me%lambda_esc, &
                                    me%eta_smear, &
                                    me%gauge, &
                                    qphot=me%qphot)
            end do
            !$OMP END DO
            !$OMP END PARALLEL
         end if

      end do

      deallocate (epsk, Hk, vectk)
      if (allocated(Sk)) deallocate (Sk)

      call batch_diag%Clean()

   end subroutine CalcPES_2D
!--------------------------------------------------------------------------------------
   subroutine CalcPES_Slab(me)
      class(arpes_calc_t) :: me
      integer :: ik, iepe
      real(dp) :: kpar(2), kqpar(2), kpt(3), kqpt(3)
      real(dp), allocatable, dimension(:)        :: epsk
      complex(dp), allocatable, dimension(:, :)   :: Hk, Sk, vectk
      type(Batch_Diagonalize_t) :: batch_diag

      allocate (epsk(me%nbnd), Hk(me%nbnd, me%nbnd), vectk(me%nbnd, me%nbnd))

      if (.not. me%orthogonal_basis) allocate (Sk(me%nbnd, me%nbnd))

      call batch_diag%Init(me%nbnd)

      kpt = 0.0_dp; kqpt = 0.0_dp
      do ik = 1, me%Nk
         kpar(1:2) = me%kpts(ik, 1:2)

         if(.not. me%dipole_approx) then
            kqpar(1:2) = kpar(1:2) - me%qphot(1:2)
            kqpt(1:2) = me%ham%get_kreduced(kqpar(1:2))
            Hk = me%ham%get_ham(kqpt)
            if (.not. me%orthogonal_basis) Sk = me%ovlp%get_Smat(kqpt)
         else
            kpt(1:2) = me%ham%get_kreduced(kpar(1:2))
            Hk = me%ham%get_ham(kpt)
            if (.not. me%orthogonal_basis) Sk = me%ovlp%get_Smat(kpt)
         end if

         if (me%orthogonal_basis) then
            call batch_diag%Diagonalize(Hk, epsk, vectk)
         else
            call batch_diag%DiagonalizeGen(Hk, Sk, epsk, vectk)
         end if
         epsk = epsk + me%Eshift

         if (me%lambda_mode) then
            if (allocated(me%excluded_layers)) then
               !$OMP PARALLEL
               !$OMP DO
               do iepe = 1, me%Nepe
                  me%spect(iepe, ik) = PES_Slab_Intensity( &
                                       me%ham, &
                                       me%nlayer, &
                                       me%chis, &
                                       me%lmax, &
                                       me%bessel_integ, &
                                       kpar, &
                                       me%wphot, &
                                       me%polvec, &
                                       me%Epe(iepe), &
                                       epsk, &
                                       vectk, &
                                       me%MuChem, &
                                       me%Vinner, &
                                       me%lambda_esc, &
                                       me%eta_smear, &
                                       qphot=me%qphot, &
                                       excluded_layers=me%excluded_layers)
               end do
               !$OMP END DO
               !$OMP END PARALLEL
            else
               !$OMP PARALLEL
               !$OMP DO
               do iepe = 1, me%Nepe
                  me%spect(iepe, ik) = PES_Slab_Intensity( &
                                       me%ham, &
                                       me%nlayer, &
                                       me%chis, &
                                       me%lmax, &
                                       me%bessel_integ, &
                                       kpar, &
                                       me%wphot, &
                                       me%polvec, &
                                       me%Epe(iepe), &
                                       epsk, &
                                       vectk, &
                                       me%MuChem, &
                                       me%Vinner, &
                                       me%lambda_esc, &
                                       me%eta_smear, &
                                       qphot=me%qphot)
               end do
               !$OMP END DO
               !$OMP END PARALLEL
            end if
         else
            if (allocated(me%excluded_layers)) then
               !$OMP PARALLEL
               !$OMP DO
               do iepe = 1, me%Nepe
                  me%spect(iepe, ik) = PES_Slab_Intensity( &
                                       me%orbs, &
                                       me%ham, &
                                       me%nlayer, &
                                       me%chis, &
                                       me%radints, &
                                       kpar, &
                                       me%wphot, &
                                       me%polvec, &
                                       me%Epe(iepe), &
                                       epsk, &
                                       vectk, &
                                       me%MuChem, &
                                       me%vinner, &
                                       me%lambda_esc, &
                                       me%eta_smear, &
                                       me%gauge, &
                                       qphot=me%qphot, &
                                       excluded_layers=me%excluded_layers)
               end do
               !$OMP END DO
               !$OMP END PARALLEL
            else
               !$OMP PARALLEL
               !$OMP DO
               do iepe = 1, me%Nepe
                  me%spect(iepe, ik) = PES_Slab_Intensity( &
                                       me%orbs, &
                                       me%ham, &
                                       me%nlayer, &
                                       me%chis, &
                                       me%radints, &
                                       kpar, &
                                       me%wphot, &
                                       me%polvec, &
                                       me%Epe(iepe), &
                                       epsk, &
                                       vectk, &
                                       me%MuChem, &
                                       me%vinner, &
                                       me%lambda_esc, &
                                       me%eta_smear, &
                                       me%gauge, &
                                       qphot=me%qphot)
               end do
               !$OMP END DO
               !$OMP END PARALLEL
            end if
         end if

      end do

      deallocate (epsk, Hk, vectk)
      if (allocated(Sk)) deallocate (Sk)

      call batch_diag%Clean()

   end subroutine CalcPES_Slab
!--------------------------------------------------------------------------------------
   subroutine WriteSpectrum(me,prefix)
      class(arpes_calc_t) :: me
      character(len=*),intent(in) :: prefix

#ifdef WITHHDF5
      call me%WriteSpectrum_hdf5(prefix)
#else
      call me%WriteSpectrum_txt(prefix)
#endif

   end subroutine WriteSpectrum 
!--------------------------------------------------------------------------------------
   subroutine WriteSpectrum_txt(me,prefix)
      class(arpes_calc_t) :: me
      character(len=*),intent(in) :: prefix

      call savetxt(trim(prefix)//'_pes.txt', me%spect, transp=.true.)

   end subroutine WriteSpectrum_txt
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine WriteSpectrum_hdf5(me,prefix)
      use scitools_hdf5_utils
      class(arpes_calc_t) :: me
      character(len=*),intent(in) :: prefix
      integer(HID_t) :: file_id

      call hdf_open_file(file_id, trim(prefix)//'_pes.h5', STATUS='NEW')      
      call hdf_write_dataset(file_id,'epe',me%Epe)
      call hdf_write_dataset(file_id,'spect',me%spect)     
      call hdf_close_file(file_id)

   end subroutine WriteSpectrum_hdf5
#endif
!--------------------------------------------------------------------------------------

!======================================================================================    
end module Marpes_calc

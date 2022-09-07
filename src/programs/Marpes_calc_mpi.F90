module Marpes_calc_mpi
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use mpi
   use Mdebug
   use Mdef,only: dp,iu,zero
   use Mutils,only: str
   use Mvector_bsplines,only: cplx_matrix_spline_t
   use Mlatt_utils,only: utility_Cart2Red_2D                   
   use Mham_w90,only: wann90_tb_t
   use Mwann_compress,only: PruneHoppings
   use Mwann_slab,only: Wannier_BulkToSlab
   use Mwannier_orbitals,only: wannier_orbs_t
   use Mradialwf,only: radialwf_t
   use Mscattwf,only: scattwf_t
   use Mradialintegral,only: radialinteg_t
   use Mpes_intensity,only: PES_Intensity, PES_Slab_Intensity, &
      PES_AtomicIntegrals_lambda_mpi
   use Mio_params,only: HamiltonianParams_t, PESParams_t
   use Mio_hamiltonian,only: ReadHamiltonian
   use Mio_orbitals,only: ReadWannierOrbitals
   use Marray1d_dist,only: dist_array1d_t,GetDisplSize1D
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: arpes_calc_t
!--------------------------------------------------------------------------------------
   type arpes_calc_t
      logical     :: lambda_mode=.false.,slab_mode=.false.
      integer     :: nbnd,norb
      integer     :: nlayer=0
      integer     :: gauge
      integer     :: Nepe      
      integer     :: lmax,radint_nk,radint_nr
      real(dp)    :: wphot,MuChem,Eshift,lambda_esc,eta_smear
      complex(dp) :: polvec(3)       
      integer     :: Nk,Nk_loc
      real(dp),allocatable,dimension(:)   :: Epe
      real(dp),allocatable,dimension(:,:) :: kpts,kpts_loc,spect
      type(wann90_tb_t)    :: ham
      type(wannier_orbs_t) :: orbs
      type(scattwf_t),allocatable,dimension(:)     :: chis
      type(radialinteg_t),allocatable,dimension(:) :: radints
      type(cplx_matrix_spline_t),allocatable,dimension(:) :: bessel_integ
   contains
      procedure,public  :: Init
      procedure,public  :: CalcIntegrals
      procedure,private :: CalcIntegrals_radial
      procedure,private :: CalcIntegrals_lambda
      procedure,public  :: CalcPES
      procedure,public  :: WriteSpectrum     
   end type arpes_calc_t
!--------------------------------------------------------------------------------------
   character(len=*),parameter :: fmt_info='(" Info: ",a)'                             
   integer,parameter :: gauge_len=0, gauge_mom=1
   integer,parameter :: scatt_pw=0, scatt_coul=1
   ! .. parallelization ..
   integer,parameter  :: master=0,from_master=1,from_worker=2
   integer :: ntasks,taskid,ierr
   logical :: on_root
   integer :: status(MPI_STATUS_SIZE)
   type(dist_array1d_t),private :: kdist
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
   subroutine Init(me,par_ham,par_pes,kpts)
      use Mutils,only: linspace
      class(arpes_calc_t) :: me
      type(HamiltonianParams_t),intent(in) :: par_ham
      type(PESParams_t),intent(in)         :: par_pes
      real(dp),intent(in)                  :: kpts(:,:)
      integer :: ik,ik_glob,iorb,ilay
      real(dp) :: comp_rate
      real(dp) :: kvec(3),kred(3)
      type(radialwf_t) :: rwf
      type(wann90_tb_t) :: ham_tmp

      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
      on_root = taskid == master

      call ReadHamiltonian(par_ham%file_ham,ham_tmp,file_xyz=par_ham%file_xyz)
      if(par_ham%energy_thresh > 0.0_dp) then
         call PruneHoppings(par_ham%energy_thresh,ham_tmp,me%ham,comp_rate)
         if(on_root) then
            write(output_unit,fmt_info) "compression rate: "//str(nint(100 * comp_rate)) // "%"
         end if
      else
         call me%ham%Set(ham_tmp)
      end if
      call ham_tmp%Clean()

      if(par_ham%slab_mode .and. par_ham%slab_nlayer > 0) then
         if(on_root) then 
            write(output_unit,fmt_info) "building slab with "//str(par_ham%slab_nlayer)//" layers"
         end if
         call ham_tmp%Set(me%ham)
         call Wannier_BulkToSlab(ham_tmp,par_ham%slab_nlayer,me%ham,ijmax=par_ham%slab_max_zhop)
         me%slab_mode = .true.
         me%nlayer = par_ham%slab_nlayer
      end if
      call ham_tmp%Clean()

      me%nbnd = me%ham%num_wann
      me%MuChem = par_ham%MuChem


      call ReadWannierOrbitals(par_pes%file_orbs,me%orbs)
      me%gauge = par_pes%gauge
      me%Nepe = par_pes%Nepe
      me%wphot = par_pes%wphot
      me%Eshift = par_pes%Eshift
      me%lambda_esc = par_pes%lambda_esc
      me%eta_smear = par_pes%eta_smear    
      me%polvec = par_pes%polvec
      me%radint_nk = par_pes%radint_numpoints_k
      me%radint_nr = par_pes%radint_numpoints_r
      me%lmax = par_pes%expansion_lmax
      me%lambda_mode = par_pes%lambda_orbital_term

      me%norb = me%orbs%norb

      allocate(me%Epe(me%Nepe))
      if(me%Nepe == 1) then
         allocate(me%Epe(1))
         me%Epe(1) = par_pes%Epe_min
      else
         me%Epe = linspace(par_pes%Epe_min, par_pes%Epe_max, me%Nepe)
      end if

      allocate(me%chis(me%norb))
      do iorb=1,me%norb
         call me%chis(iorb)%Init(par_pes%scatt_type,me%orbs%Zscatt(iorb))
      end do

      me%Nk = size(kpts, dim=1)
      allocate(me%kpts(me%Nk,2))
      if(par_pes%kpts_reduced) then
         do ik=1,me%Nk
            me%kpts(ik,1:2) = me%ham%recip_lattice(1,1:2) * kpts(ik,1) + &
               me%ham%recip_lattice(2,1:2) * kpts(ik,2) 
         end do
      else
         me%kpts(1:me%Nk,1:2) = kpts(1:me%Nk,1:2)
      end if

      call kdist%Init(ntasks,taskid,me%Nk)
      me%Nk_loc = kdist%N_loc(taskid)

      allocate(me%kpts_loc(me%Nk_loc,2))
      do ik=1,me%Nk_loc
         ik_glob = kdist%Indx_Loc2Glob(taskid,ik)
         me%kpts_loc(ik,1:2) = me%kpts(ik_glob,1:2)
      end do

      if(on_root) then
         select case(par_pes%gauge)
         case(gauge_len)
            write(output_unit,fmt_info) "light-matter coupling: length gauge"
         case(gauge_mom)
            write(output_unit,fmt_info) "light-matter coupling: velocity gauge"         
         end select

         select case(par_pes%scatt_type)
         case(scatt_pw)
            write(output_unit,fmt_info) "final states: plane waves"
         case(scatt_coul)
            write(output_unit,fmt_info) "final states: Coulomb waves"         
         end select      
      end if

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
      integer :: iorb,ilay
      type(radialwf_t) :: rwf
      real(dp) :: kmin,kmax

      kmin = 0.999999_dp * sqrt(2.0_dp * minval(me%Epe))
      kmax = 1.000001_dp * sqrt(2.0_dp * maxval(me%Epe))  

      allocate(me%radints(me%norb))

      do iorb=1,me%norb
         call WannOrb_to_RadialWF(me%orbs,iorb,rwf)
         call me%radints(iorb)%Init(me%orbs%L_indx(iorb),kmin,kmax,me%chis(iorb),rwf,&
            nk=me%radint_nk,gauge=me%gauge)
         call rwf%Clean()
      end do

   end subroutine CalcIntegrals_radial
!--------------------------------------------------------------------------------------
   subroutine CalcIntegrals_lambda(me)
      class(arpes_calc_t) :: me
      integer :: iorb
      type(radialwf_t) :: rwf
      real(dp) :: kmin,kmax

      kmin = 0.999999_dp * sqrt(2.0_dp * minval(me%Epe))
      kmax = 1.000001_dp * sqrt(2.0_dp * maxval(me%Epe))  

      call PES_AtomicIntegrals_lambda_mpi(me%orbs,me%chis,me%lambda_esc,me%lmax,kmin,kmax,me%bessel_integ,&
         gauge=me%gauge,Nr=me%radint_nr,Nk=me%radint_nk)

   end subroutine CalcIntegrals_lambda
!--------------------------------------------------------------------------------------
   subroutine CalcPES(me)
      use Mlinalg,only: EigH
      class(arpes_calc_t) :: me
      integer :: ik,iepe
      real(dp) :: kpar(2),kpt(3)
      real(dp),allocatable :: epsk(:)
      complex(dp),allocatable :: Hk(:,:),vectk(:,:)

      allocate(epsk(me%nbnd),Hk(me%nbnd,me%nbnd),vectk(me%nbnd,me%nbnd))
      allocate(me%spect(me%Nepe,me%Nk_loc))

      kpt = 0.0_dp
      do ik=1,me%Nk_loc
         kpar(1:2) = me%kpts_loc(ik,1:2)
         kpt(1:2) = utility_Cart2Red_2D(me%ham%recip_reduced,kpar)

         Hk = me%ham%get_ham(kpt)
         call EigH(Hk,epsk,vectk)
         epsk = epsk + me%Eshift

         print*, taskid, epsk

         if(me%slab_mode) then
            if(me%lambda_mode) then
               !$OMP PARALLEL
               !$OMP DO
               do iepe=1,me%Nepe
                  me%spect(iepe,ik) = PES_Slab_Intensity(me%ham,me%nlayer,me%chis,me%lmax,me%bessel_integ,kpar,me%wphot,&
                     me%polvec,me%Epe(iepe),epsk,vectk,me%MuChem,me%lambda_esc,me%eta_smear)
               end do
               !$OMP END DO
               !$OMP END PARALLEL
            else
               !$OMP PARALLEL
               !$OMP DO
               do iepe=1,me%Nepe
                  me%spect(iepe,ik) = PES_Slab_Intensity(me%orbs,me%ham,me%nlayer,me%chis,me%radints,kpar,me%wphot,&
                     me%polvec,me%Epe(iepe),epsk,vectk,me%MuChem,me%lambda_esc,me%eta_smear,me%gauge)
               end do
               !$OMP END DO
               !$OMP END PARALLEL
            end if
         else
            if(me%lambda_mode) then
               !$OMP PARALLEL
               !$OMP DO
               do iepe=1,me%Nepe
                  me%spect(iepe,ik) = PES_Intensity(me%ham,me%chis,me%lmax,me%bessel_integ,kpar,me%wphot,&
                     me%polvec,me%Epe(iepe),epsk,vectk,me%MuChem,me%lambda_esc,me%eta_smear)
               end do
               !$OMP END DO
               !$OMP END PARALLEL
            else
               !$OMP PARALLEL
               !$OMP DO
               do iepe=1,me%Nepe
                  me%spect(iepe,ik) = PES_Intensity(me%orbs,me%ham,me%chis,me%radints,kpar,me%wphot,&
                     me%polvec,me%Epe(iepe),epsk,vectk,me%MuChem,me%lambda_esc,me%eta_smear,me%gauge)
               end do
               !$OMP END DO
               !$OMP END PARALLEL
            end if
         end if
      end do

      deallocate(epsk,Hk,vectk)

   end subroutine CalcPES
!--------------------------------------------------------------------------------------
   subroutine WriteSpectrum(me,prefix)
      class(arpes_calc_t) :: me
      character(len=*),intent(in) :: prefix
      integer :: nsize
      integer,allocatable :: displ(:),size_loc(:)
      real(dp),allocatable :: spect(:,:)

      allocate(displ(0:ntasks-1),size_loc(0:ntasks-1))
      call GetDisplSize1D(kdist%N_loc,me%Nepe,displ,nsize,size_loc)

      if(on_root) allocate(spect(me%Nepe,me%Nk))
      call MPI_Gatherv(me%spect,nsize,MPI_DOUBLE_PRECISION,spect,size_loc,displ,&
         MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD, ierr)

      if(on_root) then
#ifdef WITHHDF5
         call WriteSpectrum_hdf5(prefix,me%Epe,spect)
#else
         call WriteSpectrum_txt(prefix,spect)
#endif
      end if

      deallocate(displ,size_loc)
      if(on_root) deallocate(spect)

   end subroutine WriteSpectrum 
!--------------------------------------------------------------------------------------
   subroutine WriteSpectrum_txt(prefix,spect)
      use Mutils,only: savetxt
      character(len=*),intent(in) :: prefix
      real(dp),intent(in) :: spect(:,:)

      call savetxt(trim(prefix)//'_pes.txt', spect, transp=.true.)

   end subroutine WriteSpectrum_txt
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine WriteSpectrum_hdf5(prefix,Epe,spect)
      use Mhdf5_utils
      character(len=*),intent(in) :: prefix
      real(dp),intent(in) :: Epe(:)
      real(dp),intent(in) :: spect(:,:)
      integer(HID_t) :: file_id

      call hdf_open_file(file_id, trim(prefix)//'_pes.h5', STATUS='NEW')      
      call hdf_write_dataset(file_id,'epe',Epe)
      call hdf_write_dataset(file_id,'spect',spect)     
      call hdf_close_file(file_id)

   end subroutine WriteSpectrum_hdf5
#endif
!--------------------------------------------------------------------------------------

!======================================================================================    
end module Marpes_calc_mpi

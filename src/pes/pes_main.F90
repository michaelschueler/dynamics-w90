module pes_main
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use scitools_def,only: dp, iu, zero, one, gauss, save_exp
   use scitools_utils,only: linspace, stop_error
   use scitools_linalg,only: get_large_size, util_matmul, util_zgemm, EigvalsHE
   use scitools_bsplines,only: spline1d_t
   use scitools_vector_bsplines,only: cplx_vector_spline_t, cplx_matrix_spline_t
   use pes_angcoeff,only: Transform_Y2X
   use pes_radialwf,only: radialwf_t
   use pes_scattwf,only: scattwf_t
   use pes_radialintegral,only: radialinteg_t
   use pes_matel,only: ScattMatrixElement_Momentum, ScattMatrixElement_Length, &
      ApplyDipoleOp, Dipole_lambda_projection, ScattMatrixElement_Lambda, &
      dipole_lambda_BesselTransform
   use wan_orbitals,only: wannier_orbs_t
   use wan_hamiltonian,only: wann90_tb_t
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: PES_MatrixElements, PES_Intensity
   public :: PES_AtomicIntegrals_lambda, PES_Intensity_besselinteg
   public :: PES_Slab_MatrixElements, PES_Slab_Intensity
#ifdef MPI
   public :: PES_AtomicIntegrals_lambda_mpi
#endif

   interface PES_MatrixElements
      module procedure PES_MatrixElements_precomp, PES_MatrixElements_besselinteg
   end interface PES_MatrixElements

   interface PES_Slab_MatrixElements
      module procedure PES_Slab_MatrixElements_precomp, PES_Slab_MatrixElements_besselinteg
   end interface PES_Slab_MatrixElements

   interface PES_Intensity
      module procedure PES_Intensity_precomp, PES_Intensity_besselinteg
   end interface PES_Intensity   

   interface PES_Slab_Intensity
      module procedure PES_Slab_Intensity_precomp, PES_Slab_Intensity_besselinteg
   end interface PES_Slab_Intensity 
!--------------------------------------------------------------------------------------
   integer,parameter :: gauge_len=0, gauge_mom=1
   integer,parameter :: wf_slater=0, wf_grid=1
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine WannOrb_to_RadialWF(orbs,iorb,rwf)
      !! This is a wrapper that converts a the radial part of the Wannier orbitals 
      !! (type wannier_orbs_t) to a radial wave-function (type radialwf_t)
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
   subroutine PES_AtomicIntegrals_lambda(orbs,scwfs,lam,lmax,kmin,kmax,bessel_integ,gauge,Nr,Nk,phi)
      type(wannier_orbs_t),intent(in) :: orbs
      type(scattwf_t),intent(in)      :: scwfs(:)
      real(dp),intent(in)             :: lam
      integer,intent(in)              :: lmax
      real(dp),intent(in)             :: kmin,kmax
      type(cplx_matrix_spline_t),allocatable,intent(out) :: bessel_integ(:)
      integer,intent(in),optional     :: gauge
      integer,intent(in),optional     :: Nr
      integer,intent(in),optional     :: Nk
      real(dp),intent(in),optional    :: phi
      integer :: gauge_,Nr_,nk_
      real(dp) :: phi_
      integer :: norb,nbnd,iorb,ibnd,mabs,idir,ncmp,sig,ir
      real(dp),allocatable :: rs(:)
      complex(dp),allocatable,dimension(:,:) :: integlm
      complex(dp),allocatable,dimension(:,:,:) :: lam_dip_Rfun_p1,lam_dip_Rfun_m1,lam_dip_Rfun
      type(radialwf_t) :: rwf
      type(cplx_matrix_spline_t) :: dipwf_m1,dipwf_p1
      type(cplx_matrix_spline_t) :: lam_dip

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      Nr_ = 256
      if(present(Nr)) Nr_ = Nr

      Nk_ = 40
      if(present(Nk)) Nk_ = Nk

      phi_ = 0.0_dp
      if(present(phi)) phi_ = phi

      norb = orbs%norb
      nbnd = norb

      allocate(bessel_integ(norb))

      ncmp = (lmax + 1)**2
      allocate(lam_dip_Rfun_m1(Nr_,ncmp,3)); lam_dip_Rfun_m1 = zero
      allocate(lam_dip_Rfun_p1(Nr_,ncmp,3)); lam_dip_Rfun_p1 = zero
      allocate(lam_dip_Rfun(Nr_,ncmp,3)); lam_dip_Rfun = zero

      do iorb=1,norb
         ! if(orbs%weight(iorb) < 1.0e-5_dp) cycle

         call WannOrb_to_RadialWF(orbs,iorb,rwf)

         rs = linspace(0.0_dp, rwf%Rmax, Nr_)

         lam_dip_Rfun_m1 = zero
         lam_dip_Rfun_p1 = zero
         lam_dip_Rfun = zero

         if(orbs%weight(iorb) > 1.0e-5_dp) then

            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               sig = 1
               if(mod(orbs%M_indx(iorb),2) /= 0) sig = -1

               ! m = |m0|
               call ApplyDipoleOp(rwf,orbs%L_indx(iorb),mabs,dipwf_m1,dipwf_p1,gauge_,phi=phi_)

               !$OMP PARALLEL PRIVATE(integlm)
               allocate(integlm(ncmp,3))
               !$OMP DO
               do ir=1,Nr_
                  call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)+1,dipwf_p1,lam,rs(ir),integlm)
                  lam_dip_Rfun_p1(ir,:,:) = integlm
               end do
               !$OMP END DO

               if(orbs%L_indx(iorb) > 0) then
                  !$OMP DO
                  do ir=1,Nr_
                     call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)-1,dipwf_m1,lam,rs(ir),integlm) 
                     lam_dip_Rfun_m1(ir,:,:) = integlm
                  end do
                  !$OMP END DO
               end if
               deallocate(integlm)
               !$OMP END PARALLEL

               lam_dip_Rfun = lam_dip_Rfun_p1 + lam_dip_Rfun_m1

               ! m = -|m0|
               if(mabs /= 0) then
                  lam_dip_Rfun_m1 = zero
                  lam_dip_Rfun_p1 = zero

                  call ApplyDipoleOp(rwf,orbs%L_indx(iorb),-mabs,dipwf_m1,dipwf_p1,gauge_,phi=phi_)

                  !$OMP PARALLEL PRIVATE(integlm)
                  allocate(integlm(ncmp,3))
                  !$OMP DO
                  do ir=1,Nr_
                     call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)+1,dipwf_p1,lam,rs(ir),integlm)
                     lam_dip_Rfun_p1(ir,:,:) = integlm
                  end do
                  !$OMP END DO

                  if(orbs%L_indx(iorb) > 0) then
                     !$OMP DO
                     do ir=1,Nr_
                        call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)-1,dipwf_m1,lam,rs(ir),integlm) 
                        lam_dip_Rfun_m1(ir,:,:) = integlm
                     end do
                      !$OMP END DO
                  end if
                  deallocate(integlm)
                  !$OMP END PARALLEL

                  ! transform to real basis
                  if(orbs%M_indx(iorb) < 0) then
                     lam_dip_Rfun = iu/sqrt(2.0_dp) * (lam_dip_Rfun_p1 + lam_dip_Rfun_m1 - sig * lam_dip_Rfun)
                  else
                     lam_dip_Rfun = 1.0_dp/sqrt(2.0_dp) * (lam_dip_Rfun_p1 + lam_dip_Rfun_m1 + sig * lam_dip_Rfun)
                  end if

               end if

            else

               call ApplyDipoleOp(rwf,orbs%L_indx(iorb),orbs%M_indx(iorb),dipwf_m1,dipwf_p1,gauge_,phi=phi_)
               !$OMP PARALLEL PRIVATE(integlm)
               allocate(integlm(ncmp,3))
               !$OMP DO
               do ir=1,Nr_
                  call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)+1,dipwf_p1,lam,rs(ir),integlm)
                  lam_dip_Rfun_p1(ir,:,:) = integlm
               end do
               !$OMP END DO

               if(orbs%L_indx(iorb) > 0) then
                  !$OMP DO
                  do ir=1,Nr_
                     call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)-1,dipwf_m1,lam,rs(ir),integlm) 
                     lam_dip_Rfun_m1(ir,:,:) = integlm
                  end do
                  !$OMP END DO
               end if
               deallocate(integlm)
               !$OMP END PARALLEL

               lam_dip_Rfun = lam_dip_Rfun_p1 + lam_dip_Rfun_m1

            end if

         end if

         ! multiply orbital weight
         lam_dip_Rfun = lam_dip_Rfun * orbs%weight(iorb)

         ! construct spline object
         call lam_dip%Init(rs, lam_dip_Rfun, ncmp, 3)

         ! radial integrals
         call dipole_lambda_BesselTransform(lmax,lam_dip,scwfs(iorb),kmin,kmax,nk_,bessel_integ(iorb))

         ! clean up for next orbital
         call rwf%Clean()
         call lam_dip%Clean()
         deallocate(rs)

      end do

      deallocate(lam_dip_Rfun_p1)
      deallocate(lam_dip_Rfun_m1)
      deallocate(lam_dip_Rfun)

   end subroutine PES_AtomicIntegrals_lambda
!--------------------------------------------------------------------------------------



!--------------------------------------------------------------------------------------
#ifdef MPI
   subroutine PES_AtomicIntegrals_lambda_mpi(orbs,scwfs,lam,lmax,kmin,kmax,bessel_integ,gauge,Nr,Nk,phi)
      use mpi
      use scitools_array1d_dist,only: dist_array1d_t,GetDisplSize1D 
      use pes_matel,only: dipole_lambda_BesselTransform_mpi
      type(wannier_orbs_t),intent(in) :: orbs
      type(scattwf_t),intent(in)      :: scwfs(:)
      real(dp),intent(in)             :: lam
      integer,intent(in)              :: lmax
      real(dp),intent(in)             :: kmin,kmax
      type(cplx_matrix_spline_t),allocatable,intent(out) :: bessel_integ(:)
      integer,intent(in),optional     :: gauge
      integer,intent(in),optional     :: Nr
      integer,intent(in),optional     :: Nk
      real(dp),intent(in),optional    :: phi
      integer :: gauge_,Nr_,nk_
      real(dp) :: phi_
      integer :: Nr_loc
      integer :: norb,nbnd,iorb,ibnd,mabs,idir,ncmp,sig,ir,ir_glob
      real(dp),allocatable :: rs(:),rs_loc(:)
      complex(dp),allocatable,dimension(:,:,:) :: lam_dip_Rfun_p1,lam_dip_Rfun_m1
      complex(dp),allocatable,dimension(:,:,:) :: lam_dip_Rfun_loc, lam_dip_Rfun_trs, lam_dip_Rfun
      type(radialwf_t) :: rwf
      type(cplx_matrix_spline_t) :: dipwf_m1,dipwf_p1
      type(cplx_matrix_spline_t) :: lam_dip
      ! .. parallelization ..
      integer,parameter  :: master=0,from_master=1,from_worker=2
      integer :: ntasks,taskid,ierr
      integer :: nsize
      integer,allocatable :: size_loc(:),displ(:)
      type(dist_array1d_t) :: rdist

      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      Nr_ = 256
      if(present(Nr)) Nr_ = Nr

      call rdist%Init(ntasks,taskid,Nr_)
      Nr_loc = rdist%N_loc(taskid)

      Nk_ = 40
      if(present(Nk)) Nk_ = Nk

      phi_ = 0.0_dp
      if(present(phi)) phi_ = phi

      norb = orbs%norb
      nbnd = norb

      allocate(bessel_integ(norb))

      ncmp = (lmax + 1)**2
      allocate(lam_dip_Rfun_m1(ncmp,3,Nr_loc)); lam_dip_Rfun_m1 = zero
      allocate(lam_dip_Rfun_p1(ncmp,3,Nr_loc)); lam_dip_Rfun_p1 = zero
      allocate(lam_dip_Rfun_loc(ncmp,3,Nr_loc)); lam_dip_Rfun_loc = zero
      allocate(lam_dip_Rfun_trs(ncmp,3,Nr_)); lam_dip_Rfun_trs = zero
      allocate(lam_dip_Rfun(Nr_,ncmp,3)); lam_dip_Rfun = zero

      allocate(displ(0:ntasks-1),size_loc(0:ntasks-1))
      call GetDisplSize1D(rdist%N_loc,ncmp*3,displ,nsize,size_loc)

      do iorb=1,norb
         ! if(orbs%weight(iorb) < 1.0e-5_dp) cycle

         call WannOrb_to_RadialWF(orbs,iorb,rwf)

         rs = linspace(0.0_dp, rwf%Rmax, Nr_)
         allocate(rs_loc(Nr_loc))
         do ir=1,Nr_loc
            ir_glob = rdist%Indx_Loc2Glob(taskid,ir)
            rs_loc(ir) = rs(ir_glob)
         end do

         lam_dip_Rfun_m1 = zero
         lam_dip_Rfun_p1 = zero
         lam_dip_Rfun_loc = zero

         if(orbs%weight(iorb) > 1.0e-5_dp) then

            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               sig = 1
               if(mod(orbs%M_indx(iorb),2) /= 0) sig = -1

               ! m = |m0|
               call ApplyDipoleOp(rwf,orbs%L_indx(iorb),mabs,dipwf_m1,dipwf_p1,gauge_,phi=phi_)

               do ir=1,Nr_loc
                  call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)+1,dipwf_p1,lam,rs_loc(ir),lam_dip_Rfun_p1(:,:,ir))
               end do

               if(orbs%L_indx(iorb) > 0) then
                  do ir=1,Nr_loc
                     call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)-1,dipwf_m1,lam,rs_loc(ir),lam_dip_Rfun_m1(:,:,ir)) 
                  end do
               end if

               lam_dip_Rfun_loc = lam_dip_Rfun_p1 + lam_dip_Rfun_m1

               ! m = -|m0|
               if(mabs /= 0) then
                  lam_dip_Rfun_m1 = zero
                  lam_dip_Rfun_p1 = zero

                  call ApplyDipoleOp(rwf,orbs%L_indx(iorb),-mabs,dipwf_m1,dipwf_p1,gauge_,phi=phi_)

                  do ir=1,Nr_loc
                     call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)+1,dipwf_p1,lam,rs_loc(ir),lam_dip_Rfun_p1(:,:,ir) )
                  end do

                  if(orbs%L_indx(iorb) > 0) then
                     do ir=1,Nr_loc
                        call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)-1,dipwf_m1,lam,rs_loc(ir),lam_dip_Rfun_m1(:,:,ir)) 
                     end do
                  end if

                  ! transform to real basis
                  if(orbs%M_indx(iorb) < 0) then
                     lam_dip_Rfun_loc = iu/sqrt(2.0_dp) * (lam_dip_Rfun_p1 + lam_dip_Rfun_m1 - sig * lam_dip_Rfun_loc)
                  else
                     lam_dip_Rfun_loc = 1.0_dp/sqrt(2.0_dp) * (lam_dip_Rfun_p1 + lam_dip_Rfun_m1 + sig * lam_dip_Rfun_loc)
                  end if

               end if

            else

               call ApplyDipoleOp(rwf,orbs%L_indx(iorb),orbs%M_indx(iorb),dipwf_m1,dipwf_p1,gauge_,phi=phi_)
               do ir=1,Nr_loc
                  call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)+1,dipwf_p1,lam,rs_loc(ir),lam_dip_Rfun_p1(:,:,ir))
               end do

               if(orbs%L_indx(iorb) > 0) then
                  do ir=1,Nr_loc
                     call Dipole_lambda_projection(lmax,orbs%L_indx(iorb)-1,dipwf_m1,lam,rs_loc(ir),lam_dip_Rfun_m1(:,:,ir)) 
                  end do
               end if

               lam_dip_Rfun_loc = lam_dip_Rfun_p1 + lam_dip_Rfun_m1

            end if

         end if

         ! multiply orbital weight
         lam_dip_Rfun_loc = lam_dip_Rfun_loc * orbs%weight(iorb)

         ! construct spline object

         call MPI_Allgatherv(lam_dip_Rfun_loc,nsize,MPI_DOUBLE_COMPLEX,lam_dip_Rfun_trs,size_loc,displ,&
            MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)

         do ir=1,Nr_
            lam_dip_Rfun(ir,:,:) = lam_dip_Rfun_trs(:,:,ir)
         end do

         call lam_dip%Init(rs, lam_dip_Rfun, ncmp, 3)

         ! radial integrals
         call dipole_lambda_BesselTransform_mpi(lmax,lam_dip,scwfs(iorb),kmin,kmax,nk_,bessel_integ(iorb))

         ! clean up for next orbital
         call rwf%Clean()
         call lam_dip%Clean()
         deallocate(rs)
         deallocate(rs_loc)

      end do

      deallocate(lam_dip_Rfun_p1)
      deallocate(lam_dip_Rfun_m1)
      deallocate(lam_dip_Rfun_loc)
      deallocate(lam_dip_Rfun_trs)
      deallocate(lam_dip_Rfun)
      deallocate(displ,size_loc)

      call rdist%Clean()

   end subroutine PES_AtomicIntegrals_lambda_mpi
#endif
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
   subroutine PES_MatrixElements_precomp(orbs,wann,scwfs,radints,kvec,vectk,lam,Matel,gauge,phi)
      type(wannier_orbs_t),intent(in)   :: orbs
      type(wann90_tb_t),intent(in)      :: wann
      type(scattwf_t),intent(in)        :: scwfs(:)
      type(radialinteg_t),intent(in)    :: radints(:)
      real(dp),intent(in)               :: kvec(3)
      complex(dp),intent(in)            :: vectk(:,:)
      real(dp),intent(in)               :: lam
      complex(dp),intent(inout)         :: matel(:,:)
      integer,intent(in),optional       :: gauge
      real(dp),intent(in),optional      :: phi
      integer :: gauge_
      real(dp) :: phi_
      logical :: large_size
      integer :: norb,nbnd,ibnd,iorb,mabs,idir
      complex(dp) :: mat_m(3),mat_mm(3)
      complex(dp),allocatable :: matomic(:,:) 
      complex(dp),allocatable :: vectk_phase(:,:)    

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      phi_ = 0.0_dp
      if(present(phi)) phi_ = phi

      norb = wann%num_wann
      large_size = get_large_size(norb)

      call assert(orbs%norb == norb, "PES_MatrixElements: orbs%norb == norb")
      call assert_shape(matel, [norb,3], "PES_MatrixElements", "matel")
      call assert_shape(vectk, [norb,norb], "PES_MatrixElements", "vectk")      

      allocate(matomic(norb,3)); matomic = zero
      do iorb=1,norb
         if(orbs%weight(iorb) < 1.0e-5_dp) cycle
         select case(gauge_)
         case(gauge_len)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               mat_m(1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),mabs,kvec,phi=phi_)
               mat_mm(1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),-mabs,kvec,phi=phi_)
               do idir=1,3
                  Matomic(iorb,idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
               end do
            else
               matomic(iorb,1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),&
                  orbs%M_indx(iorb),kvec,phi=phi_)
            end if
         case(gauge_mom)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               mat_m(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),mabs,kvec,phi=phi_)
               mat_mm(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),-mabs,kvec,phi=phi_)
               do idir=1,3
                  Matomic(iorb,idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
               end do
            else
               matomic(iorb,1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),&
                  orbs%M_indx(iorb),kvec,phi=phi_)
            end if         
         end select

         matomic(iorb,1:3) = matomic(iorb,1:3) * orbs%weight(iorb)

      end do

      call VectorPhase(norb,wann%coords,kvec,lam,vectk,vectk_phase)

      matel = zero
      call util_zgemm(vectk_phase, matomic, matel, transa_opt='T')

      deallocate(matomic)
      deallocate(vectk_phase)

   end subroutine PES_MatrixElements_precomp
!--------------------------------------------------------------------------------------
   subroutine PES_Slab_MatrixElements_precomp(orbs,wann,nlayer,scwfs,radints,kvec,vectk,lam,Matel,gauge,phi)
      real(dp) :: rthresh=-20.0_dp
      type(wannier_orbs_t),intent(in)   :: orbs
      type(wann90_tb_t),intent(in)      :: wann
      integer,intent(in)                :: nlayer
      type(scattwf_t),intent(in)        :: scwfs(:)
      type(radialinteg_t),intent(in) :: radints(:)
      real(dp),intent(in)               :: kvec(3)
      complex(dp),intent(in)            :: vectk(:,:)
      real(dp),intent(in)               :: lam
      complex(dp),intent(inout)         :: matel(:,:)
      integer,intent(in),optional       :: gauge
      real(dp),intent(in),optional      :: phi
      integer :: gauge_
      real(dp) :: phi_
      logical :: large_size
      integer :: norb,nbnd,ibnd,iorb,mabs,idir,ilay,j
      complex(dp) :: mat_m(3),mat_mm(3)
      complex(dp),allocatable :: matomic(:,:),matomic_layer(:,:) 
      complex(dp),allocatable :: vectk_phase(:,:) 

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      phi_ = 0.0_dp
      if(present(phi)) phi_ = phi

      norb = orbs%norb
      large_size = get_large_size(norb)

      nbnd = wann%num_wann

      if(norb * nlayer /= nbnd) then
         call stop_error("PES_Slab_MatrixElements_besselinteg: norb * nlayer /= nbnd")
      end if

      call assert_shape(matel, [nbnd,3], "PES_MatrixElements", "matel")
      call assert_shape(vectk, [nbnd,nbnd], "PES_MatrixElements", "vectk")      

      allocate(matomic_layer(norb,3)); matomic_layer = zero
      do iorb=1,norb
         if(orbs%weight(iorb) < 1.0e-5_dp) cycle
         select case(gauge_)
         case(gauge_len)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               mat_m(1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),mabs,kvec,phi=phi_)
               mat_mm(1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),-mabs,kvec,phi=phi_)
               do idir=1,3
                  matomic_layer(iorb,idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
               end do
            else
               matomic_layer(iorb,1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),&
                  orbs%M_indx(iorb),kvec,phi=phi_)
            end if
         case(gauge_mom)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               mat_m(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),mabs,kvec,phi=phi_)
               mat_mm(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),-mabs,kvec,phi=phi_)
               do idir=1,3
                  matomic_layer(iorb,idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
               end do
            else
               matomic_layer(iorb,1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),&
                  orbs%M_indx(iorb),kvec,phi=phi_)
            end if         
         end select

         matomic_layer(iorb,1:3) = matomic_layer(iorb,1:3) * orbs%weight(iorb)
      end do

      allocate(matomic(nbnd,3))
      do ilay=1,nlayer
         matomic((ilay-1)*norb+1:ilay*norb,1:3) = matomic_layer(1:norb,1:3)
      end do
      deallocate(matomic_layer)

      call VectorPhase(nbnd,wann%coords,kvec,lam,vectk,vectk_phase)

      matel = zero
      call util_zgemm(vectk_phase, matomic, matel, transa_opt='T')

      deallocate(matomic)
      deallocate(vectk_phase)

   end subroutine PES_Slab_MatrixElements_precomp
!--------------------------------------------------------------------------------------
   subroutine PES_MatrixElements_besselinteg(wann,scwfs,lmax,bessel_integ,kvec,vectk,lam,Matel,phi)
      type(wann90_tb_t),intent(in)           :: wann
      type(scattwf_t),intent(in)             :: scwfs(:)
      integer,intent(in)                     :: lmax
      type(cplx_matrix_spline_t),intent(in)  :: bessel_integ(:)
      real(dp),intent(in)                    :: kvec(3)
      complex(dp),intent(in)                 :: vectk(:,:)
      real(dp),intent(in)                    :: lam
      complex(dp),intent(inout)              :: matel(:,:)
      real(dp),intent(in),optional           :: phi
      real(dp) :: phi_
      logical :: large_size
      integer :: norb,nbnd,ibnd,iorb,mabs,idir
      complex(dp),allocatable :: matomic(:,:)  
      complex(dp),allocatable :: vectk_phase(:,:)  

      phi_ = 0.0_dp
      if(present(phi)) phi_ = phi

      norb = wann%num_wann
      large_size = get_large_size(norb)

      call assert_shape(matel, [norb,3], "PES_MatrixElements_besselinteg", "matel")
      call assert_shape(vectk, [norb,norb], "PES_MatrixElements_besselinteg", "vectk")      

      allocate(matomic(norb,3)); matomic = zero
      do iorb=1,norb
         matomic(iorb,1:3) = ScattMatrixElement_Lambda(lmax,scwfs(iorb),bessel_integ(iorb),kvec)
      end do

      call VectorPhase(norb,wann%coords,kvec,lam,vectk,vectk_phase)

      matel = zero
      call util_zgemm(vectk_phase, matomic, matel, transa_opt='T')

      deallocate(matomic)
      deallocate(vectk_phase)

   end subroutine PES_MatrixElements_besselinteg
!--------------------------------------------------------------------------------------
   subroutine PES_Slab_MatrixElements_besselinteg(wann,nlayer,scwfs,lmax,bessel_integ,kvec,vectk,lam,Matel)
      type(wann90_tb_t),intent(in)           :: wann
      integer,intent(in)                     :: nlayer
      type(scattwf_t),intent(in)             :: scwfs(:)
      integer,intent(in)                     :: lmax
      type(cplx_matrix_spline_t),intent(in)  :: bessel_integ(:)
      real(dp),intent(in)                    :: kvec(3)
      complex(dp),intent(in)                 :: vectk(:,:)
      real(dp),intent(in)                    :: lam
      complex(dp),intent(inout)              :: matel(:,:)
      logical :: large_size
      integer :: norb,nbnd,ibnd,iorb,mabs,idir,ilay,j
      complex(dp),allocatable :: matomic(:,:),matomic_layer(:,:)   
      complex(dp),allocatable :: vectk_phase(:,:)

      norb = size(bessel_integ,dim=1)
      large_size = get_large_size(norb)

      nbnd = wann%num_wann
      if(norb * nlayer /= nbnd) then
         call stop_error("PES_Slab_MatrixElements_besselinteg: norb * nlayer /= nbnd")
      end if
      call assert_shape(matel, [nbnd,3], "PES_Slab_MatrixElements_besselinteg", "matel")
      call assert_shape(vectk, [nbnd,nbnd], "PES_Slab_MatrixElements_besselinteg", "vectk")      

      allocate(matomic_layer(norb,3)); matomic_layer = zero
      do iorb=1,norb
         matomic_layer(iorb,1:3) = ScattMatrixElement_Lambda(lmax,scwfs(iorb),bessel_integ(iorb),kvec)
      end do

      allocate(matomic(nbnd,3))
      do ilay=1,nlayer
         matomic((ilay-1)*norb+1:ilay*norb,1:3) = matomic_layer(1:norb,1:3)
      end do

      call VectorPhase(nbnd,wann%coords,kvec,lam,vectk,vectk_phase)

      matel = zero
      call util_zgemm(vectk_phase, matomic, matel, transa_opt='T')

      deallocate(matomic_layer)
      deallocate(matomic)
      deallocate(vectk_phase)

   end subroutine PES_Slab_MatrixElements_besselinteg
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
   function PES_Intensity_precomp(orbs,wann,scwfs,radints,kpar,wphot,pol,Epe,epsk,vectk,&
      mu,lam,eta,gauge,qphot,phi) result(inten)
      type(wannier_orbs_t),intent(in)   :: orbs
      type(wann90_tb_t),intent(in)      :: wann
      type(scattwf_t),intent(in)        :: scwfs(:)
      type(radialinteg_t),intent(in)    :: radints(:)      
      real(dp),intent(in)               :: kpar(2)
      real(dp),intent(in)               :: wphot    
      complex(dp),intent(in)            :: pol(3)  
      real(dp),intent(in)               :: Epe
      real(dp),intent(in)               :: epsk(:)   
      complex(dp),intent(in)            :: vectk(:,:)        
      real(dp),intent(in)               :: mu      
      real(dp),intent(in)               :: lam       
      real(dp),intent(in)               :: eta           
      integer,intent(in),optional       :: gauge
      real(dp),intent(in),optional      :: qphot(3)
      real(dp),intent(in),optional      :: phi
      real(dp)                          :: inten
      integer :: gauge_    
      real(dp) :: phi_
      integer :: idir,nbnd,ibnd
      real(dp) :: Ez,kvec(3)
      complex(dp),allocatable :: matel(:,:),matel_pol(:)

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      phi_ = 0.0_dp
      if(present(phi)) phi_ = phi

      nbnd = wann%num_wann
      call assert(orbs%norb == nbnd, "PES_Intensity: orbs%norb == nbnd")
      call assert_shape(vectk,[nbnd,nbnd],"PES_Intensity","vectk")

      Ez = Epe - 0.5_dp * (kpar(1)**2 + kpar(2)**2)
      if(Ez < 1.0e-5_dp) then
         inten = 0.0_dp
         return
      end if

      kvec(1:2) = kpar
      kvec(3) = sqrt(2.0_dp * Ez)
      if(present(qphot)) then
         kvec = kvec - qphot
      end if
      if(kvec(3) < 1.0e-5_dp) then
         inten = 0.0_dp
         return
      end if     

      if(all(abs(epsk + wphot - Epe) > 6 * eta)) then
         inten = 0.0_dp
         return
      end if

      allocate(matel(nbnd,3),matel_pol(nbnd))

      call PES_MatrixElements(orbs,wann,scwfs,radints,kvec,vectk,lam,matel,gauge=gauge_,phi=phi_)

      matel_pol = zero
      do idir=1,3
         matel_pol(:) = matel_pol(:) + pol(idir) * matel(:,idir)
      end do

      inten = 0.0_dp
      do ibnd=1,nbnd
         if(epsk(ibnd) > mu) cycle
         inten = inten + abs(matel_pol(ibnd))**2 * gauss(eta, epsk(ibnd) + wphot - Epe)
      end do

      deallocate(matel,matel_pol)


   end function PES_Intensity_precomp
!--------------------------------------------------------------------------------------
   function PES_Slab_Intensity_precomp(orbs,wann,nlayer,scwfs,radints,kpar,wphot,pol,Epe,&
      epsk,vectk,mu,lam,eta,gauge,qphot,phi) result(inten)
      type(wannier_orbs_t),intent(in)   :: orbs
      type(wann90_tb_t),intent(in)      :: wann
      integer,intent(in)                :: nlayer
      type(scattwf_t),intent(in)        :: scwfs(:)
      type(radialinteg_t),intent(in)    :: radints(:)      
      real(dp),intent(in)               :: kpar(2)
      real(dp),intent(in)               :: wphot    
      complex(dp),intent(in)            :: pol(3)  
      real(dp),intent(in)               :: Epe
      real(dp),intent(in)               :: epsk(:)   
      complex(dp),intent(in)            :: vectk(:,:)        
      real(dp),intent(in)               :: mu      
      real(dp),intent(in)               :: lam       
      real(dp),intent(in)               :: eta           
      integer,intent(in),optional       :: gauge
      real(dp),intent(in),optional      :: qphot(3)
      real(dp),intent(in),optional      :: phi
      real(dp)                          :: inten
      integer :: gauge_    
      real(dp) :: phi_
      integer :: idir,norb,nbnd,ibnd
      real(dp) :: Ez,kvec(3)
      complex(dp),allocatable :: matel(:,:),matel_pol(:)

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      phi_ = 0.0_dp
      if(present(phi)) phi_ = phi

      norb = orbs%norb
      nbnd = wann%num_wann

      if(norb * nlayer /= nbnd) then
         call stop_error("PES_Slab_MatrixElements_besselinteg: norb * nlayer /= nbnd")
      end if
      call assert_shape(vectk,[nbnd,nbnd],"PES_Slab_Intensity_precomp","vectk")

      Ez = Epe - 0.5_dp * (kpar(1)**2 + kpar(2)**2)
      if(Ez < 1.0e-5_dp) then
         inten = 0.0_dp
         return
      end if

      kvec(1:2) = kpar
      kvec(3) = sqrt(2.0_dp * Ez)
      if(present(qphot)) kvec = kvec - qphot
      if(kvec(3) < 1.0e-5_dp) then
         inten = 0.0_dp
         return
      end if     

      if(all(abs(epsk + wphot - Epe) > 6 * eta)) then
         inten = 0.0_dp
         return
      end if

      allocate(matel(nbnd,3),matel_pol(nbnd))

      call PES_Slab_MatrixElements(orbs,wann,nlayer,scwfs,radints,kvec,vectk,lam,matel,gauge=gauge_,phi=phi_)

      matel_pol = zero
      do idir=1,3
         matel_pol(:) = matel_pol(:) + pol(idir) * matel(:,idir)
      end do

      inten = 0.0_dp
      do ibnd=1,nbnd
         if(epsk(ibnd) > mu) cycle
         inten = inten + abs(matel_pol(ibnd))**2 * gauss(eta, epsk(ibnd) + wphot - Epe)   
      end do

      deallocate(matel,matel_pol)


   end function PES_Slab_Intensity_precomp
!--------------------------------------------------------------------------------------
   function PES_Intensity_besselinteg(wann,scwfs,lmax,bessel_integ,kpar,wphot,pol,Epe,epsk,&
      vectk,mu,lam,eta,qphot,phi) result(inten)
      type(wann90_tb_t),intent(in)           :: wann
      type(scattwf_t),intent(in)             :: scwfs(:)
      integer,intent(in)                     :: lmax
      type(cplx_matrix_spline_t),intent(in)  :: bessel_integ(:)   
      real(dp),intent(in)                    :: kpar(2)
      real(dp),intent(in)                    :: wphot    
      complex(dp),intent(in)                 :: pol(3)  
      real(dp),intent(in)                    :: Epe
      real(dp),intent(in)                    :: epsk(:)   
      complex(dp),intent(in)                 :: vectk(:,:)        
      real(dp),intent(in)                    :: mu      
      real(dp),intent(in)                    :: lam       
      real(dp),intent(in)                    :: eta  
      real(dp),intent(in),optional           :: qphot(3)  
      real(dp),intent(in),optional           :: phi      
      real(dp)                               :: inten
      real(dp) :: phi_
      integer :: idir,nbnd,ibnd
      real(dp) :: Ez,kvec(3)
      complex(dp),allocatable :: matel(:,:),matel_pol(:)

      phi_ = 0.0_dp
      if(present(phi)) phi_ = phi

      nbnd = wann%num_wann
      call assert_shape(vectk,[nbnd,nbnd],"PES_Intensity_besselinteg","vectk")

      Ez = Epe - 0.5_dp * (kpar(1)**2 + kpar(2)**2)
      if(Ez < 1.0e-5_dp) then
         inten = 0.0_dp
         return
      end if

      kvec(1:2) = kpar
      kvec(3) = sqrt(2.0_dp * Ez)
      if(present(qphot)) kvec = kvec - qphot
      if(kvec(3) < 1.0e-5_dp) then
         inten = 0.0_dp
         return
      end if     

      if(all(abs(epsk + wphot - Epe) > 6 * eta)) then
         inten = 0.0_dp
         return
      end if

      allocate(matel(nbnd,3),matel_pol(nbnd))

      call PES_MatrixElements(wann,scwfs,lmax,bessel_integ,kvec,vectk,lam,Matel,phi=phi_)

      matel_pol = zero
      do idir=1,3
         matel_pol(:) = matel_pol(:) + pol(idir) * matel(:,idir)
      end do

      inten = 0.0_dp
      do ibnd=1,nbnd
         if(epsk(ibnd) > mu) cycle
         inten = inten + abs(matel_pol(ibnd))**2 * gauss(eta, epsk(ibnd) + wphot - Epe)   
      end do

      deallocate(matel,matel_pol)


   end function PES_Intensity_besselinteg
!--------------------------------------------------------------------------------------
   function PES_Slab_Intensity_besselinteg(wann,nlayer,scwfs,lmax,bessel_integ,kpar,wphot,pol,Epe,epsk,&
      vectk,mu,lam,eta,qphot,phi) result(inten)
      type(wann90_tb_t),intent(in)           :: wann
      integer,intent(in)                     :: nlayer
      type(scattwf_t),intent(in)             :: scwfs(:)
      integer,intent(in)                     :: lmax
      type(cplx_matrix_spline_t),intent(in)  :: bessel_integ(:)   
      real(dp),intent(in)                    :: kpar(2)
      real(dp),intent(in)                    :: wphot    
      complex(dp),intent(in)                 :: pol(3)  
      real(dp),intent(in)                    :: Epe
      real(dp),intent(in)                    :: epsk(:)   
      complex(dp),intent(in)                 :: vectk(:,:)        
      real(dp),intent(in)                    :: mu      
      real(dp),intent(in)                    :: lam       
      real(dp),intent(in)                    :: eta      
      real(dp),intent(in),optional           :: qphot(3) 
      real(dp),intent(in),optional           :: phi    
      real(dp)                               :: inten
      real(dp) :: phi_
      integer :: idir,nbnd,ibnd,norb
      real(dp) :: Ez,kvec(3)
      complex(dp),allocatable :: matel(:,:),matel_pol(:)

      phi_ = 0.0_dp
      if(present(phi)) phi_ = phi

      norb = size(bessel_integ,dim=1)
      nbnd = wann%num_wann
      call assert_shape(vectk,[nbnd,nbnd],"PES_Slab_Intensity_besselinteg","vectk")

      Ez = Epe - 0.5_dp * (kpar(1)**2 + kpar(2)**2)
      if(Ez < 1.0e-5_dp) then
         inten = 0.0_dp
         return
      end if

      kvec(1:2) = kpar
      kvec(3) = sqrt(2.0_dp * Ez)
      if(present(qphot)) kvec = kvec - qphot
      if(kvec(3) < 1.0e-5_dp) then
         inten = 0.0_dp
         return
      end if     

      if(all(abs(epsk + wphot - Epe) > 6 * eta)) then
         inten = 0.0_dp
         return
      end if

      allocate(matel(nbnd,3),matel_pol(nbnd))

      call PES_Slab_MatrixElements(wann,nlayer,scwfs,lmax,bessel_integ,kvec,vectk,lam,Matel)

      matel_pol = zero
      do idir=1,3
         matel_pol(:) = matel_pol(:) + pol(idir) * matel(:,idir)
      end do

      inten = 0.0_dp
      do ibnd=1,nbnd
         if(epsk(ibnd) > mu) cycle
         inten = inten + abs(matel_pol(ibnd))**2 * gauss(eta, epsk(ibnd) + wphot - Epe)   
      end do

      deallocate(matel,matel_pol)


   end function PES_Slab_Intensity_besselinteg
!--------------------------------------------------------------------------------------
   subroutine VectorPhase(norb,coords,kvec,lam,vectk,vectk_phase)
      real(dp),parameter :: rthresh=-20.0_dp
      integer,intent(in)  :: norb
      real(dp),intent(in) :: coords(:,:)
      real(dp),intent(in) :: kvec(3)
      real(dp),intent(in) :: lam
      complex(dp),intent(in) :: vectk(:,:)
      complex(dp),allocatable,intent(out) :: vectk_phase(:,:)
      integer :: j,ibnd
      real(dp) :: z0,phi,xlam

      allocate(vectk_phase(norb,norb)); vectk_phase = zero

      z0 = maxval(coords(:,3))
      do j=1,norb
         phi = dot_product(kvec,coords(j,1:3))
         xlam = lam * (coords(j,3) - z0)
         ! if(xlam < rthresh) cycle
         do ibnd=1,norb
            vectk_phase(j,ibnd) = exp(-iu * phi) * save_exp(xlam) * vectk(j,ibnd) 
         end do
      end do

   end subroutine VectorPhase
!--------------------------------------------------------------------------------------

!======================================================================================
end module pes_main
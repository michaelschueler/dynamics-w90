module Mpes_intensity
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp, iu, zero, one, gauss
   use Mutils,only: linspace
   use Mlinalg,only: get_large_size, util_matmul, util_zgemm, EigvalsHE, EigHE
   use Mbsplines,only: spline1d_t
   use Mvector_bsplines,only: cplx_vector_spline_t, cplx_matrix_spline_t
   use Mangcoeff,only: Transform_Y2X
   use Mradialwf,only: radialwf_t
   use Mscattwf,only: scattwf_t
   use Mradialintegral,only: radialinteg_t
   use Mmatrix_elements,only: ScattMatrixElement_Momentum, ScattMatrixElement_Length, &
      ApplyDipoleOp, Dipole_lambda_projection, ScattMatrixElement_Lambda, &
      dipole_lambda_BesselTransform
   use Mwannier_orbitals,only: wannier_orbs_t
   use Mham_w90,only: wann90_tb_t
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: PES_MatrixElements, PES_Intensity, PES_Intensity_escapedepth
   public :: PES_AtomicIntegrals_lambda, PES_Intensity_besselinteg

   interface PES_MatrixElements
      module procedure PES_MatrixElements_comp, PES_MatrixElements_precomp, PES_MatrixElements_besselinteg
   end interface PES_MatrixElements

   interface PES_Intensity
      module procedure PES_Intensity_comp, PES_Intensity_precomp, PES_Intensity_besselinteg
   end interface PES_Intensity   
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
   subroutine PES_AtomicIntegrals_lambda(orbs,scwfs,lam,lmax,kmin,kmax,bessel_integ,gauge,Nr,Nk)
      use Mtime,only: Timer_Tic,Timer_toc
      type(wannier_orbs_t),intent(in) :: orbs
      type(scattwf_t),intent(in)      :: scwfs(:)
      real(dp),intent(in)             :: lam
      integer,intent(in)              :: lmax
      real(dp),intent(in)             :: kmin,kmax
      type(cplx_matrix_spline_t),allocatable,intent(out) :: bessel_integ(:)
      integer,intent(in),optional     :: gauge
      integer,intent(in),optional     :: Nr
      integer,intent(in),optional     :: Nk
      integer :: gauge_,Nr_,nk_
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
               call ApplyDipoleOp(rwf,orbs%L_indx(iorb),mabs,dipwf_m1,dipwf_p1,gauge_)

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

                  call ApplyDipoleOp(rwf,orbs%L_indx(iorb),-mabs,dipwf_m1,dipwf_p1,gauge_)

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

               call ApplyDipoleOp(rwf,orbs%L_indx(iorb),orbs%M_indx(iorb),dipwf_m1,dipwf_p1,gauge_)
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
   subroutine PES_MatrixElements_comp(orbs,wann,scwfs,kvec,vectk,lam,Matel,gauge)
      type(wannier_orbs_t),intent(in) :: orbs
      type(wann90_tb_t),intent(in)    :: wann
      type(scattwf_t),intent(in)      :: scwfs(:)
      real(dp),intent(in)             :: kvec(3)
      complex(dp),intent(inout)       :: vectk(:,:)
      real(dp),intent(in)             :: lam
      complex(dp),intent(inout)       :: matel(:,:)
      integer,intent(in),optional     :: gauge
      integer :: gauge_
      logical :: large_size
      integer :: norb,nbnd,iorb,ibnd,mabs,idir
      real(dp) :: phi,z0
      complex(dp) :: mat_m(3),mat_mm(3)
      complex(dp),allocatable :: matomic(:,:)    
      type(radialwf_t) :: rwf

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      large_size = get_large_size(norb)

      norb = wann%num_wann
      nbnd = norb
      call assert(orbs%norb == norb, "PES_MatrixElements: orbs%norb == norb")
      call assert_shape(matel, [norb,3], "PES_MatrixElements", "matel")
      call assert_shape(vectk, [norb,norb], "PES_MatrixElements", "vectk")      

      allocate(matomic(norb,3)); matomic = zero
      do iorb=1,norb
         if(orbs%weight(iorb) < 1.0e-5_dp) cycle
         call WannOrb_to_RadialWF(orbs,iorb,rwf)

         select case(gauge_)
         case(gauge_len)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               mat_m(1:3) = ScattMatrixElement_Length(scwfs(iorb),rwf,orbs%L_indx(iorb),mabs,kvec)
               mat_mm(1:3) = ScattMatrixElement_Length(scwfs(iorb),rwf,orbs%L_indx(iorb),-mabs,kvec)
               do idir=1,3
                  Matomic(iorb,idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
               end do
            else
               matomic(iorb,1:3) = ScattMatrixElement_Length(scwfs(iorb),rwf,orbs%L_indx(iorb),&
                  orbs%M_indx(iorb),kvec)
            end if
         case(gauge_mom)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               mat_m(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),rwf,orbs%L_indx(iorb),mabs,kvec)
               mat_mm(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),rwf,orbs%L_indx(iorb),-mabs,kvec)
               do idir=1,3
                  Matomic(iorb,idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
               end do
            else
               matomic(iorb,1:3) = ScattMatrixElement_Momentum(scwfs(iorb),rwf,orbs%L_indx(iorb),&
                  orbs%M_indx(iorb),kvec)
            end if         
         end select

         matomic(iorb,1:3) = matomic(iorb,1:3) * orbs%weight(iorb)

         call rwf%Clean()

      end do

      z0 = maxval(wann%coords(:,3))

      do iorb=1,norb
         phi = dot_product(kvec,wann%coords(iorb,1:3))
         do ibnd=1,nbnd
            vectk(iorb,ibnd) = exp(-iu * phi) * exp(lam * (wann%coords(iorb,3) - z0) ) * vectk(iorb,ibnd) 
         end do
      end do

      matel = zero
      call util_zgemm(vectk, matomic, matel, transa_opt='T')
      ! call util_matmul(vectk(1:norb,1:norb), matomic(1:norb,1:3), matel(1:norb,1:3), &
      !    large_size=large_size)

      deallocate(matomic)

   end subroutine PES_MatrixElements_comp
!--------------------------------------------------------------------------------------
   subroutine PES_MatrixElements_precomp(orbs,wann,scwfs,radints,kvec,vectk,lam,Matel,gauge)
      type(wannier_orbs_t),intent(in)   :: orbs
      type(wann90_tb_t),intent(in)      :: wann
      type(scattwf_t),intent(in)        :: scwfs(:)
      type(radialinteg_t),intent(in) :: radints(:)
      real(dp),intent(in)               :: kvec(3)
      complex(dp),intent(inout)         :: vectk(:,:)
      real(dp),intent(in)               :: lam
      complex(dp),intent(inout)         :: matel(:,:)
      integer,intent(in),optional       :: gauge
      integer :: gauge_
      logical :: large_size
      integer :: norb,nbnd,ibnd,iorb,mabs,idir
      real(dp) :: phi,z0
      complex(dp) :: mat_m(3),mat_mm(3)
      complex(dp),allocatable :: matomic(:,:)    

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      large_size = get_large_size(norb)

      norb = wann%num_wann
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
               mat_m(1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),mabs,kvec)
               mat_mm(1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),-mabs,kvec)
               do idir=1,3
                  Matomic(iorb,idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
               end do
            else
               matomic(iorb,1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),&
                  orbs%M_indx(iorb),kvec)
            end if
         case(gauge_mom)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               mat_m(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),mabs,kvec)
               mat_mm(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),-mabs,kvec)
               do idir=1,3
                  Matomic(iorb,idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
               end do
            else
               matomic(iorb,1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),&
                  orbs%M_indx(iorb),kvec)
            end if         
         end select

         matomic(iorb,1:3) = matomic(iorb,1:3) * orbs%weight(iorb)

      end do

      z0 = maxval(wann%coords(:,3))

      do iorb=1,norb
         phi = dot_product(kvec,wann%coords(iorb,1:3))
         do ibnd=1,nbnd
            vectk(iorb,ibnd) = exp(-iu * phi) * exp(lam * (wann%coords(iorb,3) - z0) ) * vectk(iorb,ibnd) 
         end do
      end do

      matel = zero
      call util_zgemm(vectk, matomic, matel, transa_opt='T')


      deallocate(matomic)

   end subroutine PES_MatrixElements_precomp
!--------------------------------------------------------------------------------------
   subroutine PES_MatrixElements_besselinteg(wann,scwfs,lmax,bessel_integ,kvec,vectk,lam,Matel)
      type(wann90_tb_t),intent(in)           :: wann
      type(scattwf_t),intent(in)             :: scwfs(:)
      integer,intent(in)                     :: lmax
      type(cplx_matrix_spline_t),intent(in)  :: bessel_integ(:)
      real(dp),intent(in)                    :: kvec(3)
      complex(dp),intent(inout)              :: vectk(:,:)
      real(dp),intent(in)                    :: lam
      complex(dp),intent(inout)              :: matel(:,:)
      logical :: large_size
      integer :: norb,nbnd,ibnd,iorb,mabs,idir
      real(dp) :: phi,z0
      complex(dp),allocatable :: matomic(:,:)    

      large_size = get_large_size(norb)

      norb = wann%num_wann
      call assert_shape(matel, [norb,3], "PES_MatrixElements_besselinteg", "matel")
      call assert_shape(vectk, [norb,norb], "PES_MatrixElements_besselinteg", "vectk")      

      allocate(matomic(norb,3)); matomic = zero
      do iorb=1,norb
         matomic(iorb,1:3) = ScattMatrixElement_Lambda(lmax,scwfs(iorb),bessel_integ(iorb),kvec)
      end do

      z0 = maxval(wann%coords(:,3))

      do iorb=1,norb
         phi = dot_product(kvec,wann%coords(iorb,1:3))
         do ibnd=1,nbnd
            vectk(iorb,ibnd) = exp(-iu * phi) * exp(lam * (wann%coords(iorb,3) - z0) ) * vectk(iorb,ibnd) 
         end do
      end do

      matel = zero
      call util_zgemm(vectk, matomic, matel, transa_opt='T')

      deallocate(matomic)

   end subroutine PES_MatrixElements_besselinteg
!--------------------------------------------------------------------------------------
   subroutine MatrixElements_Atomic_escape(orbs,wann,scwfs,radints,kvec2d,kperps,lam,lam_order,matel_atomic,gauge)
      type(wannier_orbs_t),intent(in)   :: orbs
      type(wann90_tb_t),intent(in)      :: wann
      type(scattwf_t),intent(in)        :: scwfs(:)
      type(radialinteg_t),intent(in)    :: radints(:)
      real(dp),intent(in)               :: kvec2d(2)
      real(dp),intent(in)               :: kperps(:)
      real(dp),intent(in)               :: lam
      integer,intent(in)                :: lam_order
      complex(dp),intent(inout)         :: matel_atomic(:,:,:)
      integer,intent(in),optional       :: gauge
      integer :: gauge_
      logical :: large_size
      integer :: norb,nk,ik,iorb,mabs,idir,n
      real(dp) :: kvec(3),derr,deri
      complex(dp) :: xiln,mat_m(3),mat_mm(3),mat_r(3)
      real(dp),allocatable :: melr(:,:),meli(:,:)
      integer :: spl_order,iflag,inbvx(2)
      type(spline1d_t) :: splr,spli

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      large_size = get_large_size(norb)

      norb = wann%num_wann
      nk = size(kperps)
      call assert(orbs%norb == norb, "MatrixElements_Atomic_escape: orbs%norb == norb")
      call assert_shape(matel_atomic, [norb,3,nk], "MatrixElements_Atomic_escape", "matel_atomic")

      spl_order = lam_order + 2

      allocate(melr(nk,3),meli(nk,3))
      matel_atomic = zero

      do iorb=1,norb
         if(orbs%weight(iorb) < 1.0e-5_dp) cycle

         select case(gauge_)
         case(gauge_len)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               !$OMP PARALLEL DO PRIVATE(kvec,mat_m,mat_mm,mat_r)
               do ik=1,nk
                  kvec(1:2) = kvec2d(1:2)
                  kvec(3) = kperps(ik)
                  mat_m(1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),mabs,kvec)
                  mat_mm(1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),-mabs,kvec)
                  do idir=1,3
                     mat_r(idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
                     melr(ik,idir) = dble(mat_r(idir))
                     meli(ik,idir) = aimag(mat_r(idir))
                  end do       
               end do
               !$OMP END PARALLEL DO 
            else
               !$OMP PARALLEL DO PRIVATE(kvec,mat_m,mat_mm,mat_r)
               do ik=1,nk
                  kvec(1:2) = kvec2d(1:2)
                  kvec(3) = kperps(ik)
                  mat_r(1:3) = ScattMatrixElement_Length(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),&
                     orbs%M_indx(iorb),kvec)
                  do idir=1,3
                     melr(ik,idir) = dble(mat_r(idir))
                     meli(ik,idir) = aimag(mat_r(idir))
                  end do              
               end do    
               !$OMP END PARALLEL DO
            end if
         case(gauge_mom)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               !$OMP PARALLEL DO PRIVATE(kvec,mat_m,mat_mm,mat_r)
               do ik=1,nk
                  kvec(1:2) = kvec2d(1:2)
                  kvec(3) = kperps(ik)               
                  mat_m(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),mabs,kvec)
                  mat_mm(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),-mabs,kvec)
                  do idir=1,3
                     mat_r(idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
                     melr(ik,idir) = dble(mat_r(idir))
                     meli(ik,idir) = aimag(mat_r(idir))
                  end do
               end do
               !$OMP END PARALLEL DO
            else
               !$OMP PARALLEL DO PRIVATE(kvec,mat_m,mat_mm,mat_r)
               do ik=1,nk
                  kvec(1:2) = kvec2d(1:2)
                  kvec(3) = kperps(ik)               
                  mat_r(1:3) = ScattMatrixElement_Momentum(scwfs(iorb),radints(iorb),orbs%L_indx(iorb),&
                     orbs%M_indx(iorb),kvec)
                  do idir=1,3
                     melr(ik,idir) = dble(mat_r(idir))
                     meli(ik,idir) = aimag(mat_r(idir))
                  end do                     
               end do
               !$OMP END PARALLEL DO
            end if         
         end select

         do idir=1,3
            iflag = 0
            call splr%Init(kperps,melr(:,idir),spl_order,iflag)

            iflag = 0
            call spli%Init(kperps,meli(:,idir),spl_order,iflag)     

            xiln = one
            do n=0,lam_order
               if(n == 0) then
                  xiln = one
               else
                  xiln = xiln * iu * lam / dble(n)
               end if

               inbvx = 1
               do ik=1,nk
                  derr = splr%Eval(kperps(ik),n,iflag,inbvx(1))
                  deri = spli%Eval(kperps(ik),n,iflag,inbvx(2))
                  matel_atomic(iorb,idir,ik) = matel_atomic(iorb,idir,ik) + xiln * cmplx(derr,deri,kind=dp)
               end do

            end do

            call splr%Clean()
            call spli%Clean()

         end do
      end do   

      deallocate(melr,meli)  

   end subroutine MatrixElements_Atomic_escape
!--------------------------------------------------------------------------------------
   function PES_Intensity_comp(orbs,wann,scwfs,kpar,wphot,pol,Epe,epsk,vectk,mu,lam,eta,gauge) result(int)
      type(wannier_orbs_t),intent(in) :: orbs
      type(wann90_tb_t),intent(in)    :: wann
      type(scattwf_t),intent(in)      :: scwfs(:)
      real(dp),intent(in)             :: kpar(2)
      real(dp),intent(in)             :: wphot    
      complex(dp),intent(in)          :: pol(3)  
      real(dp),intent(in)             :: Epe
      real(dp),intent(in)             :: epsk(:)   
      complex(dp),intent(inout)       :: vectk(:,:)        
      real(dp),intent(in)             :: mu      
      real(dp),intent(in)             :: lam       
      real(dp),intent(in)             :: eta           
      integer,intent(in),optional     :: gauge
      real(dp)                        :: int
      integer :: gauge_    
      integer :: idir,nbnd,ibnd
      real(dp) :: Ez,kvec(3)
      complex(dp),allocatable :: matel(:,:),matel_pol(:)

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      nbnd = wann%num_wann
      call assert(orbs%norb == nbnd, "PES_Intensity: orbs%norb == nbnd")
      call assert_shape(vectk,[nbnd,nbnd],"PES_Intensity","vectk")

      Ez = Epe - 0.5_dp * (kpar(1)**2 + kpar(2)**2)
      if(Ez < 1.0e-5_dp) then
         int = 0.0_dp
         return
      end if

      kvec(1:2) = kpar
      kvec(3) = sqrt(2.0_dp * Ez)


      if(all(abs(epsk + wphot - Epe) > 6 * eta)) then
         int = 0.0_dp
         return
      end if

      allocate(matel(nbnd,3),matel_pol(nbnd))

      call PES_MatrixElements(orbs,wann,scwfs,kvec,vectk,lam,matel,gauge=gauge_)

      matel_pol = zero
      do idir=1,3
         matel_pol(:) = matel_pol(:) + pol(idir) * matel(:,idir)
      end do

      int = 0.0_dp
      do ibnd=1,nbnd
         if(epsk(ibnd) > mu) cycle
         int = int + abs(matel_pol(ibnd))**2 * gauss(eta, epsk(ibnd) + wphot - Epe)   
      end do

      deallocate(matel,matel_pol)


   end function PES_Intensity_comp
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
   function PES_Intensity_precomp(orbs,wann,scwfs,radints,kpar,wphot,pol,Epe,epsk,vectk,mu,lam,eta,gauge) result(inten)
      type(wannier_orbs_t),intent(in)   :: orbs
      type(wann90_tb_t),intent(in)      :: wann
      type(scattwf_t),intent(in)        :: scwfs(:)
      type(radialinteg_t),intent(in)    :: radints(:)      
      real(dp),intent(in)               :: kpar(2)
      real(dp),intent(in)               :: wphot    
      complex(dp),intent(in)            :: pol(3)  
      real(dp),intent(in)               :: Epe
      real(dp),intent(in)               :: epsk(:)   
      complex(dp),intent(inout)         :: vectk(:,:)        
      real(dp),intent(in)               :: mu      
      real(dp),intent(in)               :: lam       
      real(dp),intent(in)               :: eta           
      integer,intent(in),optional       :: gauge
      real(dp)                          :: inten
      integer :: gauge_    
      integer :: idir,nbnd,ibnd
      real(dp) :: Ez,kvec(3)
      complex(dp),allocatable :: matel(:,:),matel_pol(:)

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

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


      if(all(abs(epsk + wphot - Epe) > 6 * eta)) then
         inten = 0.0_dp
         return
      end if

      allocate(matel(nbnd,3),matel_pol(nbnd))

      call PES_MatrixElements(orbs,wann,scwfs,radints,kvec,vectk,lam,matel,gauge=gauge_)

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
   function PES_Intensity_besselinteg(wann,scwfs,lmax,bessel_integ,kpar,wphot,pol,Epe,epsk,&
      vectk,mu,lam,eta) result(inten)
      type(wann90_tb_t),intent(in)           :: wann
      type(scattwf_t),intent(in)             :: scwfs(:)
      integer,intent(in)                     :: lmax
      type(cplx_matrix_spline_t),intent(in)  :: bessel_integ(:)   
      real(dp),intent(in)                    :: kpar(2)
      real(dp),intent(in)                    :: wphot    
      complex(dp),intent(in)                 :: pol(3)  
      real(dp),intent(in)                    :: Epe
      real(dp),intent(in)                    :: epsk(:)   
      complex(dp),intent(inout)              :: vectk(:,:)        
      real(dp),intent(in)                    :: mu      
      real(dp),intent(in)                    :: lam       
      real(dp),intent(in)                    :: eta           
      real(dp)                               :: inten
      integer :: idir,nbnd,ibnd
      real(dp) :: Ez,kvec(3)
      complex(dp),allocatable :: matel(:,:),matel_pol(:)

      nbnd = wann%num_wann
      call assert_shape(vectk,[nbnd,nbnd],"PES_Intensity_besselinteg","vectk")

      Ez = Epe - 0.5_dp * (kpar(1)**2 + kpar(2)**2)
      if(Ez < 1.0e-5_dp) then
         inten = 0.0_dp
         return
      end if

      kvec(1:2) = kpar
      kvec(3) = sqrt(2.0_dp * Ez)

      if(all(abs(epsk + wphot - Epe) > 6 * eta)) then
         inten = 0.0_dp
         return
      end if

      allocate(matel(nbnd,3),matel_pol(nbnd))

      call PES_MatrixElements(wann,scwfs,lmax,bessel_integ,kvec,vectk,lam,Matel)

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

!--------------------------------------------------------------------------------------
   subroutine PES_Intensity_escapedepth(orbs,wann,scwfs,radints,kpar,wphot,pol,Epe,epsk,&
         vectk,mu,lam,eta,inten,gauge,lam_order)
      type(wannier_orbs_t),intent(in)   :: orbs
      type(wann90_tb_t),intent(in)      :: wann
      type(scattwf_t),intent(in)        :: scwfs(:)
      type(radialinteg_t),intent(in)    :: radints(:)      
      real(dp),intent(in)               :: kpar(2)
      real(dp),intent(in)               :: wphot    
      complex(dp),intent(in)            :: pol(3)  
      real(dp),intent(in)               :: Epe(:)
      real(dp),intent(in)               :: epsk(:)   
      complex(dp),intent(in)            :: vectk(:,:)        
      real(dp),intent(in)               :: mu      
      real(dp),intent(in)               :: lam       
      real(dp),intent(in)               :: eta  
      real(dp),intent(inout)            :: inten(:)         
      integer,intent(in),optional       :: gauge
      integer,intent(in),optional       :: lam_order      
      integer :: gauge_,lam_order_
      integer :: nk,ik,idir,nbnd,norb,ibnd,iorb
      real(dp) :: Ez,z0,phi,kvec(3)
      real(dp),allocatable :: kperps(:)
      complex(dp),allocatable :: matel_atomic(:,:,:),matel(:,:),matel_pol(:),vect_phase(:,:)

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      lam_order_ = 4
      if(present(lam_order)) lam_order_ = lam_order

      norb = wann%num_wann
      nbnd = wann%num_wann
      nk = size(epe)
      call assert(orbs%norb == nbnd, "PES_Intensity_escapedepth: orbs%norb == nbnd")
      call assert_shape(epsk,[nbnd],"PES_Intensity_escapedepth","epsk")
      call assert_shape(vectk,[nbnd,nbnd],"PES_Intensity_escapedepth","vectk")

      allocate(kperps(nk))
      do ik=1,nk
         Ez = Epe(ik) - 0.5_dp * (kpar(1)**2 + kpar(2)**2)
         if(Ez < 1.0e-5_dp) then
            kperps(ik) = -1.0_dp
         else
            kperps(ik) = sqrt(2.0_dp * Ez)
         end if
      end do

      allocate(matel_atomic(nbnd,3,nk))

      call MatrixElements_Atomic_escape(orbs,wann,scwfs,radints,kpar,kperps,lam,&
         lam_order_,matel_atomic,gauge=gauge_)

      allocate(matel(nbnd,3),matel_pol(nbnd),vect_phase(nbnd,nbnd))

      z0 = maxval(wann%coords(:,3))

      inten = 0.0_dp

      do ik=1,nk
         if(kperps(ik) > 0.0_dp) then
            kvec = [kpar(1), kpar(2), kperps(ik)]

            do iorb=1,norb
               phi = dot_product(kvec,wann%coords(iorb,1:3))
               do ibnd=1,nbnd
                  vect_phase(iorb,ibnd) = exp(-iu * phi) * exp(lam * (wann%coords(iorb,3) - z0) ) * vectk(iorb,ibnd) 
               end do
            end do

            call util_zgemm(vect_phase, matel_atomic(:,:,ik), matel, transa_opt='T')
         else
            matel = zero
         end if

         matel_pol = zero
         do idir=1,3
            matel_pol(1:nbnd) = matel_pol(1:nbnd) + pol(idir) * matel(:,idir)
         end do

         do ibnd=1,nbnd
            if(epsk(ibnd) > mu) cycle
            inten(ik) = inten(ik) + abs(matel_pol(ibnd))**2 * gauss(eta, epsk(ibnd) + wphot - Epe(ik))   
         end do         
      end do

      deallocate(kperps)
      deallocate(matel_atomic)
      deallocate(matel,matel_pol)

   end subroutine PES_Intensity_escapedepth
!--------------------------------------------------------------------------------------

!======================================================================================
end module Mpes_intensity
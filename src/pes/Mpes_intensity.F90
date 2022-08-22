module Mpes_intensity
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp, iu, zero, gauss
   use Mlinalg,only: get_large_size, util_matmul, EigvalsHE, EigHE
   use Mangcoeff,only: Transform_Y2X
   use Mradialwf,only: radialwf_t
   use Mscattwf,only: scattwf_t
   use Mmatrix_elements,only: ScattMatrixElement_Momentum, ScattMatrixElement_Length
   use Mwannier_orbitals,only: wannier_orbs_t
   use Mham_w90,only: wann90_tb_t
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: PES_MatrixElements, PES_Intensity
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
         call rwf%InitSlater(orbs%Zeff(iorb),orbs%N_indx(iorb),orbs%L_indx(iorb))
      else
         call rwf%InitGrid(orbs%rs,orbs%Rrad(:,iorb))
      end if

   end subroutine WannOrb_to_RadialWF
!--------------------------------------------------------------------------------------
   subroutine PES_MatrixElements(orbs,wann,scattwf,kvec,vectk,lam,Matel,gauge)
      type(wannier_orbs_t),intent(in) :: orbs
      type(wann90_tb_t),intent(in)    :: wann
      type(scattwf_t),intent(in)      :: scattwf
      real(dp),intent(in)             :: kvec(3)
      complex(dp),intent(inout)       :: vectk(:,:)
      real(dp),intent(in)             :: lam
      complex(dp),intent(inout)       :: matel(:,:)
      integer,intent(in),optional     :: gauge
      integer :: gauge_
      logical :: large_size
      integer :: norb,iorb,mabs,idir
      real(dp) :: phi,z0
      complex(dp) :: mat_m(3),mat_mm(3)
      complex(dp),allocatable :: matomic(:,:)    
      type(radialwf_t) :: rwf

      gauge_ = gauge_len
      if(present(gauge)) gauge_ = gauge

      large_size = get_large_size(norb)

      norb = wann%num_wann
      call assert(orbs%norb == norb, "PES_MatrixElements: orbs%norb == norb")
      call assert_shape(matel, [norb,3], "PES_MatrixElements", "matel")
      call assert_shape(vectk, [norb,norb], "PES_MatrixElements", "vectk")      

      allocate(matomic(norb,3)); matomic = zero
      do iorb=1,norb

         call WannOrb_to_RadialWF(orbs,iorb,rwf)

         select case(gauge_)
         case(gauge_len)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               mat_m(1:3) = ScattMatrixElement_Length(scattwf,rwf,orbs%L_indx(iorb),mabs,kvec)
               mat_mm(1:3) = ScattMatrixElement_Length(scattwf,rwf,orbs%L_indx(iorb),-mabs,kvec)
               do idir=1,3
                  Matomic(iorb,idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
               end do
            else
               matomic(iorb,1:3) = ScattMatrixElement_Length(scattwf,rwf,orbs%L_indx(iorb),&
                  orbs%M_indx(iorb),kvec)
            end if
         case(gauge_mom)
            if(orbs%real_lm) then
               mabs = abs(orbs%M_indx(iorb))
               mat_m(1:3) = ScattMatrixElement_Momentum(scattwf,rwf,orbs%L_indx(iorb),mabs,kvec)
               mat_mm(1:3) = ScattMatrixElement_Momentum(scattwf,rwf,orbs%L_indx(iorb),-mabs,kvec)
               do idir=1,3
                  Matomic(iorb,idir) = Transform_Y2X(orbs%M_indx(iorb),mat_m(idir),mat_mm(idir))
               end do
            else
               matomic(iorb,1:3) = ScattMatrixElement_Momentum(scattwf,rwf,orbs%L_indx(iorb),&
                  orbs%M_indx(iorb),kvec)
            end if         
         end select

         call rwf%Clean()

      end do

      z0 = maxval(wann%coords(:,3))

      do iorb=1,norb
         phi = dot_product(kvec,wann%coords(iorb,1:3))
         vectk(1:norb,iorb) = exp(-iu * phi) * exp(lam * (wann%coords(iorb,3) - z0) ) * vectk(1:norb,iorb) 
      end do

      matel = zero
      call util_matmul(vectk(1:norb,1:norb), matomic(1:norb,1:3), matel(1:norb,1:3), &
         large_size=large_size)

      deallocate(matomic)

   end subroutine PES_MatrixElements
!--------------------------------------------------------------------------------------
   function PES_Intensity(orbs,wann,scattwf,kpar,wphot,pol,Epe,epsk,vectk,mu,lam,eta,gauge) result(int)
      type(wannier_orbs_t),intent(in) :: orbs
      type(wann90_tb_t),intent(in)    :: wann
      type(scattwf_t),intent(in)      :: scattwf
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

      call PES_MatrixElements(orbs,wann,scattwf,kvec,vectk,lam,matel,gauge=gauge_)

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


   end function PES_Intensity
!--------------------------------------------------------------------------------------

!======================================================================================
end module Mpes_intensity
module wan_fourier
!! Provides discrete Fourier transforms to construct k-space operators from TB 
!! representation.
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero,one
   implicit none
   include '../units_inc.f90'
   include '../formats.h'
!--------------------------------------------------------------------------------------
   private 
   public :: &
      fourier_R_to_k, &
      fourier_R_to_k_deriv, &
      fourier_R_to_k_slab, &
      fourier_D2_R_to_k, &
      fourier_R_to_k_truevec,&
      fourier_k_to_R
!--------------------------------------------------------------------------------------
   integer,parameter :: blocksize=32
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k(kpt, irvec, ndegen, OO_R, OO)
      !! Performs the Fourier transformation R -> k
      !! For \(\alpha=0\): 
      !! \(O_{ij}(R) \rightarrow O_{ij}(k) = \sum_R e^{i k \cdot R} O_{ij}(R)\)
      !! For \(\alpha=1,2,3\):
      !! \(i \sum_R R_\alpha e^{i k \cdot R} O_{ij}(R) \)
      real(kind=dp)                                     :: kpt(3) !! k-point (reduced coordinates)
      integer,intent(in)                                :: irvec(:,:) !! hopping vectors
      integer,intent(in)                                :: ndegen(:) !! degeneracy factors
      complex(kind=dp), dimension(:, :, :), intent(in)  :: OO_R !! operator in real space O(R)
      complex(kind=dp), dimension(:, :), intent(inout)  :: OO !! operator in k-space O(k)

      integer          :: nrpts,num_wann,ir,nn,nr,m
      complex(kind=dp) :: phase_fac(blocksize)
      integer :: numblock,imin,imax

      nrpts = size(irvec, dim=1)
      num_wann = size(OO_R, dim=1)
      call assert_shape(OO_R, [num_wann, num_wann, nrpts], "fourier_R_to_k", "OO_R")
      call assert_shape(OO, [num_wann, num_wann], "fourier_R_to_k", "OO")

      ! compute the number of chuncks
      numblock  = (nrpts+blocksize-1)/blocksize
      nn = num_wann**2

      OO(:, :) = zero
      do m = 1, numblock
         imin = (m-1)*blocksize+1
         imax = min(m * blocksize, nrpts)
         nr = imax-imin+1
         call GetPhase(nr,kpt,irvec(imin:imax,:),ndegen(imin:imax),phase_fac)

         call ZGEMV("N",nn,nr,one,OO_R(1,1,imin),nn,phase_fac(1),1,one,OO(1,1),1)  
      end do

   end subroutine fourier_R_to_k
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k_deriv(kpt, irvec, crvec, ndegen, OO_R, OO)
      !! Performs the Fourier transformation R -> k
      !! For \(\alpha=0\): 
      !! \(O_{ij}(R) \rightarrow O_{ij}(k) = \sum_R e^{i k \cdot R} O_{ij}(R)\)
      !! For \(\alpha=1,2,3\):
      !! \(i \sum_R R_\alpha e^{i k \cdot R} O_{ij}(R) \)
      real(kind=dp)                                     :: kpt(3) !! k-point (reduced coordinates)
      integer,intent(in)                                :: irvec(:,:) !! hopping vectors
      real(dp),intent(in)                               :: crvec(:,:) !! cartesian hopping vectors
      integer,intent(in)                                :: ndegen(:) !! degeneracy factors
      complex(kind=dp), dimension(:, :, :), intent(in)  :: OO_R !! operator in real space O(R)
      complex(kind=dp), dimension(:, :, :), intent(inout)  :: OO !! operator in k-space O(k)

      integer          :: nrpts,num_wann,m,nn,nr,ir,i,j,idir
      integer :: numblock,imin,imax
      complex(kind=dp) :: phase_fac(blocksize),r_phase(blocksize)

      nrpts = size(irvec, dim=1)
      num_wann = size(OO_R, dim=1)
      call assert_shape(OO_R, [num_wann, num_wann, nrpts], "fourier_R_to_k_deriv", "OO_R")
      call assert_shape(OO, [num_wann, num_wann, 3], "fourier_R_to_k_deriv", "OO")

      ! compute the number of chuncks
      numblock  = (nrpts+blocksize-1)/blocksize
      nn = num_wann**2

      OO(:, :, :) = zero
      do m = 1, numblock
         imin = (m-1)*blocksize+1
         imax = min(m * blocksize, nrpts)
         nr = imax-imin+1
         call GetPhase(nr,kpt,irvec(imin:imax,:),ndegen(imin:imax),phase_fac)

         do idir=1,3
            r_phase(1:nr) = iu * crvec(imin:imax,idir) * phase_fac(1:nr)
            call ZGEMV("N",nn,nr,one,OO_R(1,1,imin),nn,r_phase(1),1,one,OO(1,1,idir),1)
         end do
      end do

   end subroutine fourier_R_to_k_deriv
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k_slab(ijmax, kpt, irvec, crvec, ndegen, OO_R, OO, alpha)
      !! 2D Fourier Transformation for the slab calculation
      !!
      !! For alpha=0:
      !! O_ij(R) --> O_ij(k) = sum_R e^{+ik.R}*O_ij(R)
      !!
      !! For alpha=1,2,3:
      !! sum_R [cmplx_i*R_alpha*e^{+ik.R}*O_ij(R)]
      !! where R_alpha is a Cartesian component of R
      !! ***REMOVE EVENTUALLY*** (replace with pw90common_fourier_R_to_k_new)
      !!
      !! TO DO:
      !! (1) transformation for random direction of irvec (U)
      !! (2) slab(logic) --> input parameters into the type w90
      ! Arguments
      !
      integer,intent(in)                                            :: ijmax
      real(kind=dp)                                                 :: kpt(2)
      integer,intent(in)                                            :: irvec(:,:) !! hopping vectors
      real(dp),intent(in)                                           :: crvec(:,:) !! cartesian hopping vectors
      integer,intent(in)                                            :: ndegen(:) !! degeneracy factors
      complex(kind=dp), dimension(:, :, :), intent(in)              :: OO_R
      complex(kind=dp), dimension(:, :, :), intent(inout)           :: OO
      integer                                                       :: alpha

      integer          :: nrpts, ir, i, j, ideg, i3
      real(kind=dp)    :: rdotk
      complex(kind=dp) :: phase_fac

      nrpts = size(irvec, dim=1)

      OO(:, :, :) = zero
      do ir = 1, nrpts
         rdotk = DPI*dot_product(kpt(:), irvec(ir, :))
         phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir), dp)
         i3 = irvec(ir, 3)
         if (abs(i3) < ijmax) then
             if (alpha == 0) then
                OO(:, :,i3+ijmax+1) = OO(:, :,i3+ijmax+1) + &
                    phase_fac*OO_R(:, :, ir)
             elseif (alpha == 1 .or. alpha == 2 .or. alpha == 3) then
                OO(:, :,i3+ijmax+1) = OO(:, :,i3+ijmax+1) + &
                   iu*crvec(ir, alpha)*phase_fac*OO_R(:, :, ir)
             else
                stop 'wrong value of alpha in fourier_R_to_k_2D'
             end if
         end if

      end do

   end subroutine fourier_R_to_k_slab
!--------------------------------------------------------------------------------------
   subroutine fourier_D2_R_to_k(kpt, irvec, crvec, ndegen, OO_R, OO, a, b)                                                     !
      !! sum_R [- R_a * R_b * e^{+ik.R}*O_ij(R)]
      !! where R_a is a Cartesian component of R

      ! Arguments
      !
      real(kind=dp)                                     :: kpt(3)
      integer,intent(in)                                :: irvec(:,:) !! hopping vectors
      real(dp),intent(in)                               :: crvec(:,:) !! cartesian hopping vectors
      integer,intent(in)                                :: ndegen(:) !! degeneracy factors
      complex(kind=dp), dimension(:, :, :), intent(in)  :: OO_R
      complex(kind=dp), dimension(:, :), intent(inout)  :: OO
      integer,intent(in)                                :: a,b

      integer          :: nrpts, ir, i, j, ideg
      real(kind=dp)    :: rdotk
      complex(kind=dp) :: phase_fac

      nrpts = size(irvec, dim=1)

      OO(:, :) = zero
      do ir = 1, nrpts
         rdotk = DPI*dot_product(kpt(:), irvec(ir, :))
         phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir), dp)
         OO(:, :) = OO(:, :) - &
            crvec(ir, a)*crvec(ir, b)*phase_fac*OO_R(:, :, ir)
      end do

   end subroutine fourier_D2_R_to_k
!--------------------------------------------------------------------------------------
   subroutine GetPhase(nr,kpt,irvecs,dgens,exp_iphase) 
      integer,intent(in)  :: nr
      real(dp),intent(in) :: kpt(3)
      integer,intent(in)  :: irvecs(:,:)
      integer,intent(in)  :: dgens(:)
      complex(dp),intent(inout) :: exp_iphase(:)
      real(dp),dimension(blocksize) :: s,c,rdotk

      rdotk(1:nr) = DPI*(kpt(1) * irvecs(1:nr,1) + kpt(2) * irvecs(1:nr,2) &
         + kpt(3) * irvecs(1:nr,3))
      c(1:nr) = cos(rdotk(1:nr))
      s(1:nr) = sin(rdotk(1:nr))
      exp_iphase(1:nr) = cmplx(c(1:nr), s(1:nr), kind=dp) / dgens(1:nr)

   end subroutine GetPhase
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k_truevec(kpt, irvec, ndegen, OO_R, OO_true)                                                              !
      !! For OO_true (true vector):
      !! $${\vec O}_{ij}(k) = \sum_R e^{+ik.R} {\vec O}_{ij}(R)$$
      ! Arguments
      !
      real(kind=dp)                                     :: kpt(3)
      integer,intent(in)                                :: irvec(:,:) !! hopping vectors
      integer,intent(in)                                :: ndegen(:) !! degeneracy factors
      complex(kind=dp), dimension(:, :, :, :), intent(in)  :: OO_R
      complex(kind=dp), dimension(:, :, :), intent(inout)   :: OO_true

      integer          :: nrpts,num_wann,nr,ir,m,nn
      complex(kind=dp) :: phase_fac(blocksize)
      integer :: numblock,imin,imax

      nrpts = size(irvec, dim=1)
      num_wann = size(OO_R, dim=1)
      call assert_shape(OO_R, [num_wann, num_wann, 3, nrpts], "fourier_R_to_k_truevec", "OO_R")
      call assert_shape(OO_true, [num_wann, num_wann, 3], "fourier_R_to_k_truevec", "OO_true")

      ! compute the number of chuncks
      numblock  = (nrpts+blocksize-1)/blocksize
      nn = 3 * (num_wann)**2

      OO_true = zero

      do m = 1, numblock
         imin = (m-1)*blocksize+1
         imax = min(m * blocksize, nrpts)
         nr = imax-imin+1
         call GetPhase(nr,kpt,irvec(imin:imax,:),ndegen(imin:imax),phase_fac)
         ! do ir=imin,imax
         !    OO_true(:,:,:) = OO_true(:,:,:) + phase_fac(ir-imin+1) * OO_R(:,:,:,ir)
         ! end do
         call ZGEMV("N",nn,nr,one,OO_R(1,1,1,imin),nn,phase_fac(1),1,one,OO_true(1,1,1),1)
      end do

   end subroutine fourier_R_to_k_truevec
!--------------------------------------------------------------------------------------
   subroutine fourier_R_to_k_vec(kpt, irvec, crvec, ndegen, OO_R, OO_true, OO_pseudo)                                                               !
      !! For OO_true (true vector):
      !! $${\vec O}_{ij}(k) = \sum_R e^{+ik.R} {\vec O}_{ij}(R)$$

      ! Arguments
      !
      real(kind=dp)                                     :: kpt(3)
      integer,intent(in)                                :: irvec(:,:) !! hopping vectors
      real(dp),intent(in)                               :: crvec(:,:) !! cartesian hopping vectors
      integer,intent(in)                                :: ndegen(:) !! degeneracy factors
      complex(kind=dp), dimension(:, :, :, :), intent(in)  :: OO_R
      complex(kind=dp), optional, dimension(:, :, :), intent(inout)   :: OO_true
      complex(kind=dp), optional, dimension(:, :, :), intent(inout)   :: OO_pseudo

      integer          :: nrpts, ir, i, j, ideg
      real(kind=dp)    :: rdotk
      complex(kind=dp) :: phase_fac

      nrpts = size(irvec, dim=1)

      if (present(OO_true)) OO_true = zero
      if (present(OO_pseudo)) OO_pseudo = zero

      do ir = 1, nrpts
         rdotk = DPI*dot_product(kpt(:), irvec(ir, :))
         phase_fac = cmplx(cos(rdotk), sin(rdotk), dp)/real(ndegen(ir), dp)

         if (present(OO_true)) then
            OO_true(:, :, 1) = OO_true(:, :, 1) + phase_fac*OO_R(:, :, 1, ir)
            OO_true(:, :, 2) = OO_true(:, :, 2) + phase_fac*OO_R(:, :, 2, ir)
            OO_true(:, :, 3) = OO_true(:, :, 3) + phase_fac*OO_R(:, :, 3, ir)
         end if
         if (present(OO_pseudo)) then
            OO_pseudo(:, :, 1) = OO_pseudo(:, :, 1) &
               + iu*crvec(ir, 2)*phase_fac*OO_R(:, :, 3, ir) &
               - iu*crvec(ir, 3)*phase_fac*OO_R(:, :, 2, ir)
            OO_pseudo(:, :, 2) = OO_pseudo(:, :, 2) &
               + iu*crvec(ir, 3)*phase_fac*OO_R(:, :, 1, ir) &
               - iu*crvec(ir, 1)*phase_fac*OO_R(:, :, 3, ir)
            OO_pseudo(:, :, 3) = OO_pseudo(:, :, 3) &
               + iu*crvec(ir, 1)*phase_fac*OO_R(:, :, 2, ir) &
               - iu*crvec(ir, 2)*phase_fac*OO_R(:, :, 1, ir)
         end if

      end do

   end subroutine fourier_R_to_k_vec
!--------------------------------------------------------------------------------------


   subroutine fourier_k_to_R(kpts, irvec, ndegen, OO_k, OO_R)
      !! Performs the Fourier transformation R -> k
      !! For \(\alpha=0\): 
      !! \(O_{ij}(R) \rightarrow O_{ij}(k) = \sum_R e^{i k \cdot R} O_{ij}(R)\)
      !! For \(\alpha=1,2,3\):
      !! \(i \sum_R R_\alpha e^{i k \cdot R} O_{ij}(R) \)
      real(kind=dp)                                     :: kpts(:,:) !! k-point grid (reduced coordinates)
      integer,intent(in)                                :: irvec(3) !! hopping vectors
      integer,intent(in)                                :: ndegen  !! degeneracy factors
      complex(kind=dp), dimension(:, :, :), intent(in)  :: OO_k !! operator in k-space O(k)
      complex(kind=dp), dimension(:, :), intent(inout)  :: OO_R !! operator in real space O(R)

      integer          :: nkpts,nrpts,num_wann,nn,nk,m
      complex(kind=dp) :: phase_fac(blocksize),degen_fac
      integer :: numblock,imin,imax

      nkpts = size(kpts, dim=1)
      nrpts = size(irvec, dim=1)
      num_wann = size(OO_k, dim=1)
      call assert_shape(OO_R, [num_wann, num_wann], "fourier_k_to_R", "OO_R")
      call assert_shape(OO_k, [num_wann, num_wann, nkpts], "fourier_k_to_R", "OO_k")

      ! compute the number of chuncks
      numblock  = (nkpts+blocksize-1)/blocksize
      nn = num_wann**2

      OO_R(:, :) = zero

      do m = 1, numblock
         imin = (m-1)*blocksize+1
         imax = min(m * blocksize, nkpts)
         nk = imax - imin + 1
         call GetInvPhase(nk,irvec,kpts(imin:imax,:),phase_fac)

         call ZGEMV("N",nn,nk,one,OO_k(1,1,imin),nn,phase_fac(1),1,one,OO_R(1,1),1)  
      end do

      degen_fac = one * ndegen / nkpts
      call ZSCAL(nn, degen_fac, OO_R(1,1), 1)


   end subroutine fourier_k_to_R
!--------------------------------------------------------------------------------------
   subroutine GetInvPhase(nk,irvec,kpts,exp_iphase) 
      integer,intent(in)   :: nk
      integer,intent(in)   :: irvec(3)
      real(dp),intent(in)  :: kpts(:,:)
      complex(dp),intent(inout) :: exp_iphase(:)
      real(dp),dimension(blocksize) :: s,c,rdotk

      rdotk(1:nk) = DPI*(kpts(1:nk,1) * irvec(1) + kpts(1:nk,2) * irvec(2) &
         + kpts(1:nk,3) * irvec(3))
      c(1:nk) = cos(rdotk(1:nk))
      s(1:nk) = -sin(rdotk(1:nk))
      exp_iphase(1:nk) = cmplx(c(1:nk), s(1:nk), kind=dp) 

   end subroutine GetInvPhase
!--------------------------------------------------------------------------------------

!======================================================================================
end module wan_fourier

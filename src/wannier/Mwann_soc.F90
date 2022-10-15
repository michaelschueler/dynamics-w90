module Mwann_soc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero
   use Mham_w90,only: wann90_tb_t
   implicit none
   include "../formats.h"
   include "../units_inc.f90"
!--------------------------------------------------------------------------------------
   private
   public :: ham_soc_t
   public :: AddSOC_Wannier
!--------------------------------------------------------------------------------------
   type ham_soc_t
      integer                                  :: ngroups,norb
      integer,allocatable,dimension(:)         :: ndim,Lorb
      complex(dp),allocatable,dimension(:,:,:) :: Lmat
   contains
      procedure,public  :: ReadFromTXT => ham_soc_ReadFromTXT
#ifdef WITHHDF5
      procedure,public  :: ReadFromHDF5 => ham_soc_ReadFromHDF5
#endif
   end type ham_soc_t
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine ham_soc_ReadFromTXT(me,fname)
      class(ham_soc_t) :: me
      character(len=*),intent(in) :: fname
      integer :: unit_inp
      logical :: file_ok
      integer :: i,j
      real(dp) :: rvec(6)


      inquire(file=trim(fname),exist=file_ok)
      if(.not.file_ok) then
         write(error_unit,fmt900) 'Input file does not exist: '//trim(fname)
         stop
      end if

      open(newunit=unit_inp, file=trim(fname), status='old', action='read')
      read(unit_inp,*) me%ngroups
      allocate(me%ndim(me%ngroups),me%Lorb(me%ngroups))
      read(unit_inp,*) me%ndim
      read(unit_inp,*) me%Lorb

      me%norb = sum(me%ndim)
      allocate(me%Lmat(me%norb,me%norb,3))
      do i=1,me%norb
         do j=1,me%norb
            read(unit_inp,*) rvec
            me%Lmat(i,j,1:3) = rvec(1:3) + iu * rvec(4:6)
         end do
      end do

      close(unit_inp)

   end subroutine ham_soc_ReadFromTXT
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine ham_soc_ReadFromHDF5(me,fname)
      use scitools_hdf5_utils
      class(ham_soc_t) :: me
      character(len=*),intent(in) :: fname
      integer(HID_t) :: file_id
      real(dp),allocatable :: rdata(:,:,:)
      logical :: file_ok

      inquire(file=trim(fname),exist=file_ok)
      if(.not.file_ok) then
         write(error_unit,fmt900) 'Input file does not exist: '//trim(fname)
         stop
      end if

      call hdf_open_file(file_id, trim(fname), STATUS='OLD', ACTION='READ')

      call hdf_read_attribute(file_id,'','ngroups', me%ngroups)
      allocate(me%ndim(me%ngroups),me%Lorb(me%ngroups))
      call hdf_read_dataset(file_id,'ndim',me%ndim)
      call hdf_read_dataset(file_id,'lorb',me%Lorb)

      me%norb = sum(me%ndim)

      allocate(me%Lmat(me%norb,me%norb,3))
      allocate(rdata(me%norb,me%norb,3))

      call hdf_read_dataset(file_id,'hsoc_real',rdata)
      me%Lmat = rdata
      call hdf_read_dataset(file_id,'hsoc_imag',rdata)
      me%Lmat = me%Lmat + iu * rdata

      deallocate(rdata)

      call hdf_close_file(file_id) 

   end subroutine ham_soc_ReadFromHDF5
#endif
!--------------------------------------------------------------------------------------
   subroutine AddSOC_Wannier(soc,lam,w90_nosoc,w90_soc)
   !! Extends the Hamiltonian \(H(\mathbf{k})\) without spin-orbit coupling (SOC) 
   !! to a new Hamiltonian \(H_s(\mathbf{k})\) that includes SOC. The new Hamiltonian
   !! is constructed by
   !! $$ [H_s(\mathbf{k})]_{\sigma \sigma^\prime} = H(\mathbf{k})\delta_{\sigma \sigma^\prime}
   !! + \sum^{n_g}_{i=1} \lambda_i \mathbf{s}_{\sigma \sigma^\prime} \cdot \mathbf{L}_i \ . $$
   !! Here, \(n_g\) is the number of groups of orbitals, where \(\mathbf{L}_i\) is the matrix
   !! representation of the orbital angular momentum in each orbital subspace with dimension \(N_i\).
   !! Typically \(N_i = 2l_i  + 1\) where \(l_i\) is the angular momentum of full shell of s,p,d,...
   !! orbitals. \(\mathbf{s}_{\sigma \sigma^\prime}\) is the vector of Pauli matrices.
   !! See the script `SOCInput.py` / `SOCInput_txt.py` for how to construct the matrices \(\mathbf{L}_i\).
      type(ham_soc_t),intent(in)    :: soc !! SOC class storing the matrices \(\mathbf{L}_i\) 
      real(dp),intent(in)           :: lam(:) !! SOC constants for each group
      type(wann90_tb_t),intent(in)  :: w90_nosoc !! Wannier Hamiltonian without SOC
      type(wann90_tb_t),intent(out) :: w90_soc !! Wannier Hamiltonian with SOC. We follow the storage format
                                               !! of `Wannier90` with `spinors=.true.`.
      integer :: norb,nbnd,i,j,idir,ir
      complex(dp),allocatable :: Hsoc(:,:)
      
      norb = soc%norb
      if(norb /= w90_nosoc%num_wann) then
         write(error_unit,fmt900) "AddSOC_Wannier: Number of orbitals does not match the number of Wannier functions."
         stop
      end if

      w90_soc%num_wann = 2 * w90_nosoc%num_wann
      w90_soc%real_lattice = w90_nosoc%real_lattice
      w90_soc%recip_lattice = w90_nosoc%recip_lattice
      w90_soc%recip_reduced = w90_nosoc%recip_reduced
      w90_soc%nrpts = w90_nosoc%nrpts
      w90_soc%coords_present = w90_nosoc%coords_present

      if(allocated(w90_soc%ndegen)) deallocate(w90_soc%ndegen)
      if(allocated(w90_soc%irvec)) deallocate(w90_soc%irvec)
      if(allocated(w90_soc%ham_r)) deallocate(w90_soc%ham_r)
      if(allocated(w90_soc%pos_r)) deallocate(w90_soc%pos_r)
      if(allocated(w90_soc%coords)) deallocate(w90_soc%coords)

      allocate(w90_soc%ndegen(w90_nosoc%nrpts))
      allocate(w90_soc%irvec(w90_nosoc%nrpts,3))
      allocate(w90_soc%ham_r(w90_soc%num_wann,w90_soc%num_wann,w90_nosoc%nrpts))
      allocate(w90_soc%pos_r(w90_soc%num_wann,w90_soc%num_wann,w90_nosoc%nrpts,3))
      allocate(w90_soc%coords(w90_soc%num_wann,3))

      w90_soc%coords = w90_nosoc%coords
      w90_soc%ndegen = w90_nosoc%ndegen
      w90_soc%irvec = w90_nosoc%irvec

      nbnd = 2 * norb
      allocate(Hsoc(nbnd,nbnd))
      call GetHam_SOC_onsite(soc,lam,Hsoc)

      do ir=1,w90_soc%nrpts
         do j=1,norb
            do i=1,norb
               w90_soc%ham_r(2*i-1,2*j-1,ir) = w90_nosoc%ham_r(i,j,ir)
               w90_soc%ham_r(2*i,2*j,ir) = w90_nosoc%ham_r(i,j,ir)
            end do
         end do
      end do

      do idir=1,3
         do ir=1,3
            do j=1,norb
               do i=1,norb
                  w90_soc%pos_r(2*i-1,2*j-1,ir,idir) = w90_nosoc%pos_r(i,j,ir,idir)
                  w90_soc%pos_r(2*i,2*j,ir,idir) = w90_nosoc%pos_r(i,j,ir,idir)
               end do
            end do
         end do
      end do

      ! reorder to Wannier90 indexing
      do ir=1,3      
         do j=1,norb
            do i=1,norb
               w90_soc%ham_r(2*i-1,2*j-1,ir) = w90_soc%ham_r(2*i-1,2*j-1,ir) + Hsoc(i,j)
               w90_soc%ham_r(2*i-1,2*j,ir) = w90_soc%ham_r(2*i-1,2*j,ir) + Hsoc(i,j+norb)
               w90_soc%ham_r(2*i,2*j-1,ir) = w90_soc%ham_r(2*i,2*j-1,ir) + Hsoc(i+norb,j)
               w90_soc%ham_r(2*i,2*j,ir) = w90_soc%ham_r(2*i,2*j,ir) + Hsoc(i+norb,j+norb)
            end do
         end do
      end do

      deallocate(Hsoc)

   end subroutine AddSOC_Wannier
!--------------------------------------------------------------------------------------
   subroutine GetHam_SOC_onsite(soc,lam,H0)
      !! Generates the on-site atomic SOC term 
      !! $$H^\mathrm{SOC}_{\sigma \sigma^\prime} = \sum^{n_g}_{i=1} \lambda_i \mathbf{s}_{\sigma \sigma^\prime} \cdot \mathbf{L}_i $$
      !! (see [[AddSOC_Wannier]]). 
      type(ham_soc_t),intent(in)    :: soc !! SOC class storing the matrices \(\mathbf{L}_i\) 
      real(dp),intent(in)           :: lam(:) !! SOC constants for each group
      complex(dp),intent(inout)     :: H0(:,:) !! atomic SOC term
      integer :: norb,nbnd,ig,imin,imax
      integer,allocatable :: orb_min(:)   

      norb = soc%norb
      nbnd = 2 * norb

      call assert_shape(lam, [soc%ngroups], "GetHam_SOC_onsite", "lam")
      call assert_shape(H0, [nbnd,nbnd], "GetHam_SOC_onsite", "H0")

      allocate(orb_min(soc%ngroups)); orb_min = 1
      do ig=1,soc%ngroups-1
         orb_min(ig+1) = orb_min(ig) + soc%ndim(ig)
      end do

      H0 = zero
      do ig=1,soc%ngroups
         imin = orb_min(ig)
         imax = orb_min(ig) + soc%ndim(ig) - 1

         ! sigma_x
         H0(imin:imax,imin+norb:imax+norb) = H0(imin:imax,imin+norb:imax+norb) &
            + lam(ig) * soc%Lmat(imin:imax,imin:imax,1)
         H0(imin+norb:imax+norb,imin:imax) = H0(imin+norb:imax+norb,imin:imax) &
            + lam(ig) * soc%Lmat(imin:imax,imin:imax,1)

         ! sigma_y
         H0(imin:imax,imin+norb:imax+norb) = H0(imin:imax,imin+norb:imax+norb) &
            - iu * lam(ig) * soc%Lmat(imin:imax,imin:imax,2)
         H0(imin+norb:imax+norb,imin:imax) = H0(imin+norb:imax+norb,imin:imax) &
            + iu * lam(ig) * soc%Lmat(imin:imax,imin:imax,2)

         ! sigma_z
         H0(imin:imax,imin:imax) = H0(imin:imax,imin:imax) &
            + lam(ig) * soc%Lmat(imin:imax,imin:imax,3)
         H0(imin+norb:imax+norb,imin+norb:imax+norb) = H0(imin+norb:imax+norb,imin+norb:imax+norb) &
            - lam(ig) * soc%Lmat(imin:imax,imin:imax,3)
      end do

   end subroutine GetHam_SOC_onsite
!--------------------------------------------------------------------------------------

!======================================================================================
end module Mwann_soc
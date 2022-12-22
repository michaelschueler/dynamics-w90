module io_obs
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use scitools_def,only: dp,iu,zero
   use scitools_utils,only: savetxt,str,save_griddata
#ifdef WITHHDF5
   use scitools_hdf5_utils
#endif
   implicit none
   include "../formats.h"  
!--------------------------------------------------------------------------------------
   type WannierCalcOutput_t
   !! Class for output of [[wann_calc]]. Controls which output is written to file.
      logical            :: calc_orbweight=.false.
      logical            :: calc_spin=.false.
      logical            :: calc_berry=.false.
      logical            :: calc_oam=.false.
      logical            :: calc_metric=.false.
      logical            :: calc_evecs=.false.
      logical            :: write_kpts=.false. 
      logical            :: write_velo=.false. 
      real(dp),pointer,dimension(:,:)        :: epsk => null()    
      real(dp),pointer,dimension(:,:)        :: kpts => null() 
      real(dp),pointer,dimension(:,:,:)      :: orb_weight => null()
      real(dp),pointer,dimension(:,:,:)      :: spin => null()
      real(dp),pointer,dimension(:,:,:)      :: berry => null()
      real(dp),pointer,dimension(:,:,:)      :: oam => null()
      real(dp),pointer,dimension(:,:,:,:)    :: metric => null()
      real(dp),pointer,dimension(:,:,:,:)    :: spin_berry => null()
      complex(dp),pointer,dimension(:,:,:)   :: evecs => null()
      complex(dp),pointer,dimension(:,:,:,:) :: velok => null()
   contains
      procedure, public  :: AddEpsk => wann_calc_AddEpsk
      procedure, public  :: AddKpts => wann_calc_AddKpts
      procedure, public  :: AddOrbWeight => wann_calc_AddOrbWeight
      procedure, public  :: AddSpin => wann_calc_AddSpin
      procedure, public  :: AddBerry => wann_calc_AddBerry
      procedure, public  :: AddSpinBerry => wann_calc_AddSpinBerry
      procedure, public  :: AddOAM => wann_calc_AddOAM
      procedure, public  :: AddMetric => wann_calc_AddMetric      
      procedure, public  :: AddEvecs => wann_calc_AddEvecs
      procedure, public  :: AddVelocity => wann_calc_AddVelocity
      procedure, public  :: SaveToFile => wann_calc_SaveToFile
      procedure, private :: SaveToFile_txt => wann_calc_SaveToFile_txt
#ifdef WITHHDF5
      procedure, private :: SaveToFile_hdf5 => wann_calc_SaveToFile_hdf5
#endif
   end type WannierCalcOutput_t
!--------------------------------------------------------------------------------------
   private
   public :: WannierCalcOutput_t, SaveTDObs, SaveTDOccupation, SaveSpinCurrent

   interface SaveTDObs
   !! Interface for saving observables from a time-dependent run with [[wann_evol]].
   !! Automatically selects the relevant observables based on the gauge of the light-matter coupling
      module procedure SaveTDObs_velo, SaveTDObs_dip
   end interface SaveTDObs
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine wann_calc_AddEpsk(me,epsk)
   !! Adds the band structure to the output
      class(WannierCalcOutput_t) :: me
      real(dp),target,intent(in) :: epsk(:,:) !! band structure, dimensions `[nbnd,nk]`

      me%epsk => epsk

   end subroutine wann_calc_AddEpsk
!--------------------------------------------------------------------------------------
   subroutine wann_calc_AddKpts(me,kpts)
   !! Adds the k-points used in the calculation to the output
      class(WannierCalcOutput_t) :: me
      real(dp),target,intent(in) :: kpts(:,:) !! k-points

      me%kpts => kpts

   end subroutine wann_calc_AddKpts
!--------------------------------------------------------------------------------------
   subroutine wann_calc_AddOrbWeight(me,orb_weight)
   !! Adds the orbital weight of each Wannier function to the output
      class(WannierCalcOutput_t) :: me
      real(dp),target,intent(in) :: orb_weight(:,:,:) !! orbital weight, dimension `[nwan,nwan,nk]`

      me%orb_weight => orb_weight

   end subroutine wann_calc_AddOrbWeight
!--------------------------------------------------------------------------------------
   subroutine wann_calc_AddSpin(me,spin)
   !! Adds the spin texture to the output
      class(WannierCalcOutput_t) :: me
      real(dp),target,intent(in) :: spin(:,:,:) !! spin texture, dimension `[nwan,3,nk]`

      me%spin => spin

   end subroutine wann_calc_AddSpin
!--------------------------------------------------------------------------------------
   subroutine wann_calc_AddBerry(me,berry)
   !! Adds the Berry curvature to the output
      class(WannierCalcOutput_t) :: me
      real(dp),target,intent(in) :: berry(:,:,:) !! Berry curvature, dimension `[nwan,3,nk]`

      me%berry => berry

   end subroutine wann_calc_AddBerry
!--------------------------------------------------------------------------------------
   subroutine wann_calc_AddSpinBerry(me,spin_berry)
   !! Adds the spin Berry curvature to the output
      class(WannierCalcOutput_t) :: me
      real(dp),target,intent(in) :: spin_berry(:,:,:,:) !! spin Berry curvature, dimension `[nwan,3,3,nk]`

      me%spin_berry => spin_berry

   end subroutine wann_calc_AddSpinBerry
!--------------------------------------------------------------------------------------
   subroutine wann_calc_AddOAM(me,oam)
   !! Adds the orbital angular momentum (OAM) texture to the output
      class(WannierCalcOutput_t) :: me
      real(dp),target,intent(in) :: oam(:,:,:) !! OAM, dimension `[nwan,3,nk]`

      me%oam => oam

   end subroutine wann_calc_AddOAM
!--------------------------------------------------------------------------------------
   subroutine wann_calc_AddMetric(me,metric)
   !! Adds the quantum metric to the output
      class(WannierCalcOutput_t) :: me
      real(dp),target,intent(in) :: metric(:,:,:,:) !! quantum metric, dimension `[nwan,3,3,nk]`

      me%metric => metric

   end subroutine wann_calc_AddMetric
!--------------------------------------------------------------------------------------
   subroutine wann_calc_AddEvecs(me,evecs)
   !! Adds the full eigenvectors at each k-point to the output
      class(WannierCalcOutput_t) :: me
      complex(dp),target,intent(in) :: evecs(:,:,:) !! eigenvectors, dimension `[nwan,nwan,nk]`

      me%evecs => evecs

   end subroutine wann_calc_AddEvecs
!--------------------------------------------------------------------------------------
   subroutine wann_calc_AddVelocity(me,velok)
   !! Adds velocity matrix elements at each k-point to the output
      class(WannierCalcOutput_t) :: me
      complex(dp),target,intent(in) :: velok(:,:,:,:) !! velocity matrix elements, dimension `[nwan,nwan,3,nk]`

      me%velok => velok

   end subroutine wann_calc_AddVelocity
!--------------------------------------------------------------------------------------
   subroutine wann_calc_SaveToFile_txt(me,prefix)
   !! Stores the output selected by the `wann_calc_AddXXX` routines to plain text file.
      class(WannierCalcOutput_t) :: me
      character(len=*),intent(in) :: prefix !! Prefix of output file name. The output files will
                                            !! be named `prefix_XXX.txt`, where `XXX` denotes the
                                            !! specific output quantity.
      integer :: nbnd,nwan,nk,i,j,ik
      integer :: unit_out
      real(dp),allocatable :: rdata(:,:)
      character(len=256) :: fout

      if(associated(me%epsk)) then
         call savetxt(trim(prefix)//'_epsk.txt', me%epsk, transp=.true.)
      end if

      if(associated(me%kpts)) then
         call savetxt(trim(prefix)//'_kpts.txt', me%kpts)
      end if

      if(associated(me%orb_weight)) then
         nwan = size(me%orb_weight, dim=1)
         do i=1,nwan
            fout = trim(prefix)//'_orbweight_'//str(i)//'.txt'
            call savetxt(trim(fout), me%orb_weight(i,:,:), transp=.true.)
         end do
      end if

      if(associated(me%spin)) then
         fout = trim(prefix)//'_spin_x.txt'
         call savetxt(trim(fout), me%spin(:,1,:), transp=.true.)
         fout = trim(prefix)//'_spin_y.txt'
         call savetxt(trim(fout), me%spin(:,2,:), transp=.true.)
         fout = trim(prefix)//'_spin_z.txt'
         call savetxt(trim(fout), me%spin(:,3,:), transp=.true.)
      end if

      if(associated(me%berry)) then
         fout = trim(prefix)//'_berry_x.txt'
         call savetxt(trim(fout), me%berry(:,1,:), transp=.true.)
         fout = trim(prefix)//'_berry_y.txt'
         call savetxt(trim(fout), me%berry(:,2,:), transp=.true.)
         fout = trim(prefix)//'_berry_z.txt'
         call savetxt(trim(fout), me%berry(:,3,:), transp=.true.)          
      end if

      if(associated(me%spin_berry)) then
         fout = trim(prefix)//'_spin_berry_x_sigx.txt'
         call savetxt(trim(fout), me%spin_berry(:,1,1,:), transp=.true.)
         fout = trim(prefix)//'_spin_berry_x_sigy.txt'
         call savetxt(trim(fout), me%spin_berry(:,1,2,:), transp=.true.)
         fout = trim(prefix)//'_spin_berry_x_sigz.txt'
         call savetxt(trim(fout), me%spin_berry(:,1,3,:), transp=.true.)     
         fout = trim(prefix)//'_spin_berry_y_sigx.txt'
         call savetxt(trim(fout), me%spin_berry(:,2,1,:), transp=.true.)
         fout = trim(prefix)//'_spin_berry_y_sigy.txt'
         call savetxt(trim(fout), me%spin_berry(:,2,2,:), transp=.true.)
         fout = trim(prefix)//'_spin_berry_y_sigz.txt'
         call savetxt(trim(fout), me%spin_berry(:,2,3,:), transp=.true.) 
         fout = trim(prefix)//'_spin_berry_z_sigx.txt'
         call savetxt(trim(fout), me%spin_berry(:,3,1,:), transp=.true.)
         fout = trim(prefix)//'_spin_berry_z_sigy.txt'
         call savetxt(trim(fout), me%spin_berry(:,3,2,:), transp=.true.)
         fout = trim(prefix)//'_spin_berry_z_sigz.txt'
         call savetxt(trim(fout), me%spin_berry(:,3,3,:), transp=.true.)        
      end if

      if(associated(me%oam)) then
         fout = trim(prefix)//'_oam_x.txt'
         call savetxt(trim(fout), me%oam(:,1,:), transp=.true.)
         fout = trim(prefix)//'_oam_y.txt'
         call savetxt(trim(fout), me%oam(:,2,:), transp=.true.)
         fout = trim(prefix)//'_oam_z.txt'
         call savetxt(trim(fout), me%oam(:,3,:), transp=.true.)   
      end if

      if(associated(me%metric)) then
         fout = trim(prefix)//'_metric_xx.txt'
         call savetxt(trim(fout), me%metric(:,1,1,:), transp=.true.)       
         fout = trim(prefix)//'_metric_xy.txt'
         call savetxt(trim(fout), me%metric(:,1,2,:), transp=.true.)   
         fout = trim(prefix)//'_metric_xz.txt'
         call savetxt(trim(fout), me%metric(:,1,3,:), transp=.true.) 
         fout = trim(prefix)//'_metric_yx.txt'
         call savetxt(trim(fout), me%metric(:,2,1,:), transp=.true.) 
         fout = trim(prefix)//'_metric_yy.txt'
         call savetxt(trim(fout), me%metric(:,2,2,:), transp=.true.) 
         fout = trim(prefix)//'_metric_yz.txt'
         call savetxt(trim(fout), me%metric(:,2,3,:), transp=.true.) 
         fout = trim(prefix)//'_metric_zx.txt'
         call savetxt(trim(fout), me%metric(:,3,1,:), transp=.true.)
         fout = trim(prefix)//'_metric_zy.txt'
         call savetxt(trim(fout), me%metric(:,3,2,:), transp=.true.) 
         fout = trim(prefix)//'_metric_zz.txt'
         call savetxt(trim(fout), me%metric(:,3,3,:), transp=.true.) 
      end if

      if(associated(me%evecs)) then
         nbnd = size(me%evecs, dim=1)
         nk = size(me%evecs, dim=3)
         allocate(rdata(2*nbnd,nk))
         do i=1,nbnd
            fout = trim(prefix)//'_evec_bnd_'//str(i)//'.txt' 
            rdata(1:nbnd,1:nk) = dble(me%evecs(1:nbnd,i,1:nk))
            rdata(nbnd+1:,1:nk) = aimag(me%evecs(1:nbnd,i,1:nk))
            call savetxt(trim(fout), rdata, transp=.true.)           
         end do
         deallocate(rdata)
      end if

      if(associated(me%velok)) then
         nbnd = size(me%velok, dim=1)
         nk = size(me%velok, dim=4)
         fout = trim(prefix)//'_velocity.txt' 
         open(newunit=unit_out,file=trim(fout),status="replace")
         do ik=1,nk
            do i=1,nbnd
               do j=1,nbnd
                  write(unit_out,'(3(I6,1x),6(F14.7,1x))') ik, i, j, dble(me%velok(i,j,1,ik)), &
                     aimag(me%velok(i,j,1,ik)), dble(me%velok(i,j,2,ik)), aimag(me%velok(i,j,2,ik)), &
                     dble(me%velok(i,j,3,ik)), aimag(me%velok(i,j,3,ik))
               end do
            end do
         end do
         close(unit_out)
      end if

   end subroutine wann_calc_SaveToFile_txt
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine wann_calc_SaveToFile_hdf5(me,prefix)
   !! Stores the output selected by the `wann_calc_AddXXX` routines to hdf5.
      class(WannierCalcOutput_t) :: me
      character(len=*),intent(in) :: prefix !! Prefix of output file name. Output will be
                                            !! written to `prefix_wann_calc.h5`.
      integer :: nbnd,nwan,nk,i
      integer(HID_t) :: file_id
      real(dp),allocatable :: rdata(:,:,:),rdata4(:,:,:,:)

      call hdf_open_file(file_id, trim(prefix)//'_wann_calc.h5', STATUS='NEW')

      if(associated(me%epsk)) then
         call hdf_write_dataset(file_id,'epsk',me%epsk)
      end if

      if(associated(me%kpts)) then
         call hdf_write_dataset(file_id,'kpts',me%kpts)
      end if

      if(associated(me%orb_weight)) then
         call hdf_write_dataset(file_id,'orbweight',me%orb_weight)
      end if

      if(associated(me%spin)) then
         call hdf_write_dataset(file_id,'spin',me%spin)
      end if

      if(associated(me%berry)) then
         call hdf_write_dataset(file_id,'berry',me%berry)       
      end if

      if(associated(me%spin_berry)) then
         call hdf_write_dataset(file_id,'spin_berry',me%spin_berry)       
      end if

      if(associated(me%oam)) then
         call hdf_write_dataset(file_id,'oam',me%oam)
      end if

      if(associated(me%metric)) then
         call hdf_write_dataset(file_id,'metric',me%metric)
      end if

      if(associated(me%evecs)) then
         nbnd = size(me%evecs, dim=1)
         nk = size(me%evecs, dim=3)
         allocate(rdata(nbnd,nbnd,nk))
         rdata = dble(me%evecs)
         call hdf_write_dataset(file_id,'evecs-real',rdata)
         rdata = aimag(me%evecs)
         call hdf_write_dataset(file_id,'evecs-imag',rdata)
         deallocate(rdata)
      end if

      if(associated(me%velok)) then
         nbnd = size(me%velok, dim=1)
         nk = size(me%velok, dim=3)
         allocate(rdata4(nbnd,nbnd,3,nk))
         rdata4 = dble(me%velok)
         call hdf_write_dataset(file_id,'velok-real',rdata4)
         rdata4 = aimag(me%velok)
         call hdf_write_dataset(file_id,'velok-imag',rdata4)
         deallocate(rdata4)
      end if      

      call hdf_close_file(file_id)

   end subroutine wann_calc_SaveToFile_hdf5
#endif
!--------------------------------------------------------------------------------------
   subroutine wann_calc_SaveToFile(me,prefix)
   !! High-level interface to store the output selected by the `wann_calc_AddXXX` routines
   !! to file. If `dynamics-w90` has been build with hdf5 support, the output will be 
   !! written in hdf5 format; otherwise, plain text format is used.
      class(WannierCalcOutput_t) :: me
      character(len=*),intent(in) :: prefix !! Prefix of output file name. Dependening
                                            !! on the output format, a single hdf5 file or
                                            !! multiple plain text files will be written.

#ifdef WITHHDF5
      call me%SaveToFile_hdf5(prefix)
#else
      call me%SaveToFile_txt(prefix)
#endif

   end subroutine wann_calc_SaveToFile
!--------------------------------------------------------------------------------------
   subroutine SaveTDObs_velo(prefix,Nt,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
      Jcurr_para,Jcurr_dia,Jcurr_intra)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt,output_step
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Etot(:)
      real(dp),intent(in) :: Ekin(:)
      real(dp),intent(in) :: BandOcc(:,:)
      real(dp),intent(in) :: Jcurr(:,:)
      real(dp),intent(in) :: Dip(:,:)
      real(dp),intent(in) :: Jcurr_para(:,:)
      real(dp),intent(in) :: Jcurr_dia(:,:)
      real(dp),intent(in) :: Jcurr_intra(:,:)

#ifdef WITHHDF5
      call SaveTDObs_velo_hdf5(prefix,Nt,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
         Jcurr_para,Jcurr_dia,Jcurr_intra)
#else
      call SaveTDObs_velo_txt(prefix,Nt,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
         Jcurr_para,Jcurr_dia,Jcurr_intra)
#endif

   end subroutine SaveTDObs_velo
!--------------------------------------------------------------------------------------
   subroutine SaveTDObs_dip(prefix,Nt,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
      Jcurr_hk,Jcurr_pol)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt,output_step
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Etot(:)
      real(dp),intent(in) :: Ekin(:)
      real(dp),intent(in) :: BandOcc(:,:)
      real(dp),intent(in) :: Jcurr(:,:)
      real(dp),intent(in) :: Dip(:,:)
      real(dp),intent(in) :: Jcurr_hk(:,:)
      real(dp),intent(in) :: Jcurr_pol(:,:)

#ifdef WITHHDF5
      call SaveTDObs_dip_hdf5(prefix,Nt,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
         Jcurr_hk,Jcurr_pol)
#else
      call SaveTDObs_dip_txt(prefix,Nt,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
         Jcurr_hk,Jcurr_pol)
#endif

   end subroutine SaveTDObs_dip
!--------------------------------------------------------------------------------------
   subroutine SaveTDOccupation(prefix,Nt,output_step,dt,Occk)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt,output_step
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Occk(:,:,0:)

#ifdef WITHHDF5
      call SaveTDOccupation_hdf5(prefix,Nt,output_step,dt,Occk)
#else
      call SaveTDOccupation_txt(prefix,Nt,Occk)
#endif

   end subroutine SaveTDOccupation
!--------------------------------------------------------------------------------------
   subroutine SaveSpinCurrent(prefix,Nt,output_step,dt,Jspin)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt,output_step
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Jspin(:,:,0:)

#ifdef WITHHDF5
      call SaveSpinCurrent_hdf5(prefix,Nt,output_step,dt,Jspin)
#else
      call SaveSpinCurrent_txt(prefix,Nt,Jspin)
#endif

   end subroutine SaveSpinCurrent
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine SaveTDObs_velo_hdf5(prefix,Nt,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
      Jcurr_para,Jcurr_dia,Jcurr_intra)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt,output_step
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Etot(:)
      real(dp),intent(in) :: Ekin(:)
      real(dp),intent(in) :: BandOcc(:,:)
      real(dp),intent(in) :: Jcurr(:,:)
      real(dp),intent(in) :: Dip(:,:)
      real(dp),intent(in) :: Jcurr_para(:,:)
      real(dp),intent(in) :: Jcurr_dia(:,:)
      real(dp),intent(in) :: Jcurr_intra(:,:)
      integer(HID_t) :: file_id
      character(len=255) :: Flname
      integer  :: tstp
      real(dp),allocatable :: ts(:)

      Flname = trim(prefix)//'_observables.h5'
      call hdf_open_file(file_id, trim(Flname), STATUS='NEW')

      allocate(ts(0:Nt))
      forall(tstp=0:Nt) ts(tstp) = dt * tstp * output_step
      call hdf_write_dataset(file_id,'time',ts)      

      call hdf_write_dataset(file_id,'etot',Etot)
      call hdf_write_dataset(file_id,'ekin',Ekin)
      call hdf_write_dataset(file_id,'occ',BandOcc)
      call hdf_write_dataset(file_id,'current',Jcurr)
      call hdf_write_dataset(file_id,'dipole',Dip)

      call hdf_write_dataset(file_id,'current_para',Jcurr_para)
      call hdf_write_dataset(file_id,'current_dia',Jcurr_dia)
      call hdf_write_dataset(file_id,'current_intra',Jcurr_intra)

      call hdf_close_file(file_id)
      deallocate(ts)

   end subroutine SaveTDObs_velo_hdf5
#endif
!--------------------------------------------------------------------------------------
   subroutine SaveTDObs_velo_txt(prefix,Nt,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
      Jcurr_para,Jcurr_dia,Jcurr_intra)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt,output_step
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Etot(:)
      real(dp),intent(in) :: Ekin(:)
      real(dp),intent(in) :: BandOcc(:,:)
      real(dp),intent(in) :: Jcurr(:,:)
      real(dp),intent(in) :: Dip(:,:)
      real(dp),intent(in) :: Jcurr_para(:,:)
      real(dp),intent(in) :: Jcurr_dia(:,:)
      real(dp),intent(in) :: Jcurr_intra(:,:)
      integer  :: tstp
      real(dp),allocatable :: ts(:)

      allocate(ts(0:Nt))
      forall(tstp=0:Nt) ts(tstp) = dt * tstp * output_step

      call save_griddata(trim(prefix)//'_etot.txt', ts, Etot)
      call save_griddata(trim(prefix)//'_ekin.txt', ts, Ekin)
      call save_griddata(trim(prefix)//'_occ.txt', ts, BandOcc, transp=.true.)
      call save_griddata(trim(prefix)//'_curr.txt', ts, Jcurr, transp=.true.)
      call save_griddata(trim(prefix)//'_dip.txt', ts, Dip, transp=.true.)     
      
      call save_griddata(trim(prefix)//'_curr_para.txt', ts, Jcurr_para, transp=.true.)
      call save_griddata(trim(prefix)//'_curr_dia.txt', ts, Jcurr_dia, transp=.true.)
      call save_griddata(trim(prefix)//'_curr_intra.txt', ts, Jcurr_intra, transp=.true.) 

      deallocate(ts)

   end subroutine SaveTDObs_velo_txt
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine SaveTDObs_dip_hdf5(prefix,Nt,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
      Jcurr_hk,Jcurr_pol)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt,output_step
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Etot(:)
      real(dp),intent(in) :: Ekin(:)
      real(dp),intent(in) :: BandOcc(:,:)
      real(dp),intent(in) :: Jcurr(:,:)
      real(dp),intent(in) :: Dip(:,:)
      real(dp),intent(in) :: Jcurr_hk(:,:)
      real(dp),intent(in) :: Jcurr_pol(:,:)
      integer(HID_t) :: file_id
      character(len=255) :: Flname
      integer  :: tstp
      real(dp),allocatable :: ts(:)

      Flname = trim(prefix)//'_observables.h5'
      call hdf_open_file(file_id, trim(Flname), STATUS='NEW')

      allocate(ts(0:Nt))
      forall(tstp=0:Nt) ts(tstp) = dt * tstp * output_step
      call hdf_write_dataset(file_id,'time',ts)      

      call hdf_write_dataset(file_id,'etot',Etot)
      call hdf_write_dataset(file_id,'ekin',Ekin)
      call hdf_write_dataset(file_id,'occ',BandOcc)
      call hdf_write_dataset(file_id,'current',Jcurr)
      call hdf_write_dataset(file_id,'dipole',Dip)

      call hdf_write_dataset(file_id,'current_hk',Jcurr_hk)
      call hdf_write_dataset(file_id,'current_pol',Jcurr_pol)

      call hdf_close_file(file_id)
      deallocate(ts)

   end subroutine SaveTDObs_dip_hdf5
#endif
!--------------------------------------------------------------------------------------
   subroutine SaveTDObs_dip_txt(prefix,Nt,output_step,dt,Etot,Ekin,BandOcc,Jcurr,Dip,&
      Jcurr_hk,Jcurr_pol)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt,output_step
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Etot(:)
      real(dp),intent(in) :: Ekin(:)
      real(dp),intent(in) :: BandOcc(:,:)
      real(dp),intent(in) :: Jcurr(:,:)
      real(dp),intent(in) :: Dip(:,:)
      real(dp),intent(in) :: Jcurr_hk(:,:)
      real(dp),intent(in) :: Jcurr_pol(:,:)
      integer  :: tstp
      real(dp),allocatable :: ts(:)

      allocate(ts(0:Nt))
      forall(tstp=0:Nt) ts(tstp) = dt * tstp * output_step

      call save_griddata(trim(prefix)//'_etot.txt', ts, Etot)
      call save_griddata(trim(prefix)//'_ekin.txt', ts, Ekin)
      call save_griddata(trim(prefix)//'_occ.txt', ts, BandOcc, transp=.true.)
      call save_griddata(trim(prefix)//'_curr.txt', ts, Jcurr, transp=.true.)
      call save_griddata(trim(prefix)//'_dip.txt', ts, Dip, transp=.true.)     
      
      call save_griddata(trim(prefix)//'_curr_hk.txt', ts, Jcurr_hk, transp=.true.)
      call save_griddata(trim(prefix)//'_curr_pol.txt', ts, Jcurr_pol, transp=.true.)

      deallocate(ts)

   end subroutine SaveTDObs_dip_txt
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine SaveTDOccupation_hdf5(prefix,Nt,output_step,dt,Occk)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt,output_step
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Occk(:,:,0:)
      integer(HID_t) :: file_id
      character(len=255) :: Flname
      integer  :: tstp
      real(dp),allocatable :: ts(:)

      Flname = trim(prefix)//'_occupation.h5'
      call hdf_open_file(file_id, trim(Flname), STATUS='NEW')

      allocate(ts(0:Nt))
      forall(tstp=0:Nt) ts(tstp) = dt * tstp * output_step
      call hdf_write_dataset(file_id,'time',ts)   

      call hdf_write_dataset(file_id,'occupation',Occk)

      call hdf_close_file(file_id)
      deallocate(ts)

   end subroutine SaveTDOccupation_hdf5
#endif
!--------------------------------------------------------------------------------------
   subroutine SaveTDOccupation_txt(prefix,Nt,Occk)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt
      real(dp),intent(in) :: Occk(:,:,0:)
      integer  :: tstp
      character(len=256) :: file_out

      do tstp=0,Nt
         file_out = trim(prefix)//'_occupation_step'//str(tstp)//'.txt'
         call savetxt(trim(file_out),Occk(:,:,tstp))
      end do

   end subroutine SaveTDOccupation_txt
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine SaveSpinCurrent_hdf5(prefix,Nt,output_step,dt,Jspin)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt,output_step
      real(dp),intent(in) :: dt
      real(dp),intent(in) :: Jspin(:,:,0:)
      integer(HID_t) :: file_id
      character(len=255) :: Flname
      integer  :: tstp
      real(dp),allocatable :: ts(:)

      Flname = trim(prefix)//'_occupation.h5'
      call hdf_open_file(file_id, trim(Flname), STATUS='NEW')

      allocate(ts(0:Nt))
      forall(tstp=0:Nt) ts(tstp) = dt * tstp * output_step
      call hdf_write_dataset(file_id,'time',ts)   

      call hdf_write_dataset(file_id,'spin-current',Jspin)

      call hdf_close_file(file_id)
      deallocate(ts)

   end subroutine SaveSpinCurrent_hdf5
#endif
!--------------------------------------------------------------------------------------
   subroutine SaveSpinCurrent_txt(prefix,Nt,Jspin)
      character(len=*),intent(in) :: prefix
      integer,intent(in)  :: Nt
      real(dp),intent(in) :: Jspin(:,:,0:)
      integer  :: tstp
      integer :: unit_out
      integer :: i,j
      character(len=256) :: file_out

      file_out = trim(prefix)//'_jspin.h5'
      open(newunit=unit_out,file=trim(file_out),STATUS='REPLACE')
      do tstp=0,Nt
         write(unit_out,'(9(F14.7,1x))') ((Jspin(i,j,tstp),i=1,3),j=1,3)
      end do
      close(unit_out)

   end subroutine SaveSpinCurrent_txt
!--------------------------------------------------------------------------------------


!======================================================================================
end module io_obs
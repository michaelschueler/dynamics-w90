program wann_calc
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdebug
   use Mdef,only: dp,zero
   use Mtime,only: Timer_Act, Timer_SetName, Timer_Run, Timer_stop, Timer_DTShow
   use Mwannier_calc,only: wannier_calc_t
   implicit none
   include '../formats.h'
!--------------------------------------------------------------------------------------
   ! -- constants --
   character(len=*),parameter :: fmt_time='(a,"/",a,"/",a," at ",a,":",a,":",a)'
   character(len=*),parameter :: fmt_info='(" Info: ",a)'
   character(len=*),parameter :: fmt_input='(" Input: [",a,"]")'
   ! -- timing --
   character(len=8) :: date
   character(len=10) :: time
   ! -- for reading i/o --
   logical :: PrintToFile = .false.
   integer :: Narg,unit_inp
   character(len=255) :: FlIn,FlOutPref
   ! -- input variables --
   logical            :: calc_orbweight=.false.
   logical            :: calc_spin=.false.
   logical            :: calc_berry=.false.
   logical            :: calc_oam=.false.
   logical            :: calc_evecs=.false.
   logical            :: write_kpts=.false.
   integer            :: gauge=0
   namelist/CALCOPT/calc_orbweight,calc_spin,calc_berry,calc_oam,calc_evecs,write_kpts,gauge
   ! -- internal variables --
   real(dp),allocatable,dimension(:,:,:) :: orb_weight,spin,berry,oam
   type(wannier_calc_t)                  :: wann
!--------------------------------------------------------------------------------------
   write(output_unit,'(A)') '+------------------------------------------------------+'
   write(output_unit,'(A)') '|             Wannier90 post-processing                |'
   write(output_unit,'(A)') '+------------------------------------------------------+'
   write(output_unit,*)

   call date_and_time(date=date,time=time)
   write(output_unit,'(a)',advance='no') '  Calculation started on '
   write(output_unit,fmt_time) date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
   write(output_unit,*)

   call Timer_Act
   call Timer_SetName('total',1)
   call Timer_Run(N=1)
!--------------------------------------------------------------------------------------
!                               ++  Read input ++
!--------------------------------------------------------------------------------------
   call Timer_SetName('Initialize Wannier', 2); call Timer_Run(N=2)

   Narg=command_argument_count()
   if(Narg>=1) then
      call get_command_argument(1,FlIn)
      open(newunit=unit_inp,file=trim(FlIn),status='OLD',action='READ')
      read(unit_inp,nml=CALCOPT)
      close(unit_inp)
   else
      write(error_unit,fmt900) 'Please provide a namelist input file.'
      stop
   end if

   if(Narg>=2) then
      call get_command_argument(2,FlOutPref)
      PrintToFile=.true.
   end if
   if(.not.PrintToFile) then
      write(output_unit,fmt_info) 'No output prefix given. No output will be produced.'
   end if

   call wann%Init(FlIn)

   write(output_unit,*)
   call Timer_Stop(N=2); call Timer_DTShow(N=2)
!--------------------------------------------------------------------------------------
!                               ++  Calculation ++
!--------------------------------------------------------------------------------------
   if(calc_orbweight) then
      call Timer_SetName('orbital weight', 2); call Timer_Run(N=2)
      call wann%GetOrbitalWeight(orb_weight)
      write(output_unit,*)
      call Timer_Stop(N=2); call Timer_DTShow(N=2)
   end if

   if(calc_spin) then
      call Timer_SetName('spin', 2); call Timer_Run(N=2)
      call wann%GetSpin(spin)
      write(output_unit,*)
      call Timer_Stop(N=2); call Timer_DTShow(N=2)
   end if

   if(calc_berry) then
      call Timer_SetName('Berry curvature', 2); call Timer_Run(N=2)
      call wann%GetBerryCurvature(berry)
      write(output_unit,*)
      call Timer_Stop(N=2); call Timer_DTShow(N=2)
   end if

   if(calc_oam) then
      call Timer_SetName('OAM', 2); call Timer_Run(N=2)
      call wann%GetOAM(oam)
      write(output_unit,*)
      call Timer_Stop(N=2); call Timer_DTShow(N=2)
   end if

   if(PrintToFile) call SaveOutput(FlOutPref)
!--------------------------------------------------------------------------------------
   write(output_unit,*)
   call Timer_Stop(N=1)
   call Timer_DTShow(N=1)
   write(output_unit,'(A)') '+------------------------------------------------------+'
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine SaveOutput(pref)
      character(len=*),intent(in) :: pref

#ifdef WITHHDF5
      call SaveOutput_hdf5(trim(pref)//'_wann_calc.h5')
#else
      call SaveOutput_txt(pref)
#endif

   end subroutine SaveOutput
!--------------------------------------------------------------------------------------
   subroutine SaveOutput_txt(fname)
      use Mutils,only: savetxt,str
      character(len=*),intent(in) :: fname
      character(len=256) :: fout
      integer :: i
      real(dp),allocatable :: rdata(:,:)

      call savetxt(trim(fname)//'_epsk.txt', wann%epsk, transp=.true.)
      if(write_kpts) call savetxt(trim(fname)//'_kpts.txt', wann%kpts)
      if(calc_orbweight) then
         do i=1,wann%nwan
            fout = trim(fname)//'_orbweight_'//str(i)//'.txt'
            call savetxt(trim(fout), orb_weight(i,:,:), transp=.true.)
         end do
      end if
      if(calc_spin) then
         fout = trim(fname)//'_spin_x.txt'
         call savetxt(trim(fout), spin(1,:,:), transp=.true.)
         fout = trim(fname)//'_spin_y.txt'
         call savetxt(trim(fout), spin(2,:,:), transp=.true.)
         fout = trim(fname)//'_spin_z.txt'
         call savetxt(trim(fout), spin(3,:,:), transp=.true.)
      end if
      if(calc_berry) then
         fout = trim(fname)//'_berry_x.txt'
         call savetxt(trim(fout), berry(:,1,:), transp=.true.)
         fout = trim(fname)//'_berry_y.txt'
         call savetxt(trim(fout), berry(:,2,:), transp=.true.)
         fout = trim(fname)//'_berry_z.txt'
         call savetxt(trim(fout), berry(:,3,:), transp=.true.)         
      end if
      if(calc_oam) then
         fout = trim(fname)//'_oam_x.txt'
         call savetxt(trim(fout), oam(:,1,:), transp=.true.)
         fout = trim(fname)//'_oam_y.txt'
         call savetxt(trim(fout), oam(:,2,:), transp=.true.)
         fout = trim(fname)//'_oam_z.txt'
         call savetxt(trim(fout), oam(:,3,:), transp=.true.)         
      end if
      if(calc_evecs) then
         allocate(rdata(2*wann%nbnd,wann%Nk))
         do i=1,wann%nbnd
            fout = trim(fname)//'_evec_bnd_'//str(i)//'.txt' 
            rdata(1:wann%nbnd,1:wann%Nk) = dble(wann%vectk(1:wann%nbnd,i,1:wann%Nk))
            rdata(wann%nbnd+1:,1:wann%Nk) = aimag(wann%vectk(1:wann%nbnd,i,1:wann%Nk))
            call savetxt(trim(fout), rdata, transp=.true.)     
         end do
         deallocate(rdata)
      end if

   end subroutine SaveOutput_txt
!--------------------------------------------------------------------------------------
#ifdef WITHHDF5
   subroutine SaveOutput_hdf5(fname)
      use Mhdf5_utils
      character(len=*),intent(in) :: fname
      integer(HID_t) :: file_id
      real(dp),allocatable :: rdata(:,:,:)
  
      call hdf_open_file(file_id, trim(fname), STATUS='NEW')
      call hdf_write_attribute(file_id,'','nk', wann%Nk)
      call hdf_write_attribute(file_id,'','nbnd', wann%nbnd)
      call hdf_write_attribute(file_id,'','nwan', wann%nwan)

      call hdf_write_dataset(file_id,'epsk',wann%epsk)

      if(write_kpts) call hdf_write_dataset(file_id,'kpts',wann%kpts)
      if(calc_orbweight) call hdf_write_dataset(file_id,'orbweight',orb_weight)
      if(calc_spin) call hdf_write_dataset(file_id,'spin',spin)
      if(calc_berry) call hdf_write_dataset(file_id,'berry',berry)
      if(calc_oam) call hdf_write_dataset(file_id,'oam',oam)

      if(calc_evecs) then
         allocate(rdata(wann%nbnd,wann%nbnd,wann%Nk))
         rdata = dble(wann%vectk)
         call hdf_write_dataset(file_id,'evecs-real',rdata)
         rdata = aimag(wann%vectk)
         call hdf_write_dataset(file_id,'evecs-imag',rdata)
         deallocate(rdata)
      end if

      call hdf_close_file(file_id)

   end subroutine SaveOutput_hdf5
#endif
!--------------------------------------------------------------------------------------

!======================================================================================
end program wann_calc
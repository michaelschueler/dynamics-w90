!======================================================================================
!****m* src/Mtime.f90
!
! NAME
!  Mtime.f90
!
! AUTHOR
!  Yaroslav Pavlyukh
!
! DESCRIPTION
!  This module provides various tools for measuring and printing
!  the execution time of parts of a source code.
!
! CONTAINS
!   Mtime.f90/Tact
!   Mtime.f90/TRun
!   Mtime.f90/TGet
!   Mtime.f90/TmsGet
!   Mtime.f90/Diff
!   Mtime.f90/Diffms
!   Mtime.f90/SetName
!   Mtime.f90/DTPrint
!   Mtime.f90/DTmsPrint
!   Mtime.f90/DTmsShow
!   Mtime.f90/Tstop
!   Mtime.f90/Int2Char
!
!****
!======================================================================================
module Mtime
  type NM
     character(40)::Name
  end type Nm
  type (nm),private,dimension(100)::Tmnm
  logical,private,dimension(100)::TmSt
  integer,private,dimension(8,100)::TmStart,TmStop,DT
  integer,private,dimension(100)::DTms
  integer,private,dimension(12),parameter::NDinM=[31,28,31,30,31,30,31,31,30,31,30,31]
  integer,private,dimension(8),parameter::C=[1,12,0,0,24,60,60,1000]

  private
  public &
       Timer_act,&
       Timer_Run,&
       TGet,&
       TmsGet,&
       Diff,&
       Diffms,&
       Timer_SetName,&
       DTPrint,&
       DTmsPrint,&
       DTmsShow,&
       Timer_DTShow,&
       Timer_stop,&
       PrintTime,&
       Timer_tic,&
       Timer_toc

Contains
!******************
!======================================================================================
!****u* Mtime.f90/Tact
!
! NAME 
!  Tact
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Activates the timer.
!
! SOURCE
! 
subroutine Timer_act
  TmSt=.FALSE.
  DT=0
  DTms=0
  Tmnm(:)%Name='     EMPTY TIMER    '
end subroutine Timer_Act
!******************
!****u* Mtime.f90/TRun
!
! NAME 
!  Tact
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Starts timer N
!
! SOURCE
! 
subroutine Timer_Run(N,NV,SNV,RunAll)
logical,intent(in)::RunAll
integer,intent(in)::N,SNV,NV(:)
optional::N,NV,SNV,RunAll
integer::CrTm(8)
call Date_And_Time(Values=CrTm)
if (Present(RunAll)) then
   if (RunAll) then
      TmSt=.True.
      do i=1,100 
         TmStart(:,i)=CrTm
      end do
   end if
   Return
end if
if (Present(SNV)) then
   do i=1,SNV
      if (.not.TmSt(NV(i))) then
         TmStart(:,NV(i))=CrTm
         TmSt(NV(i))=.True.
      else
         print*,'Timer Thread #',NV(i),' is already running'
      end if
    end do
   Return
end if
if (Present(N)) then
   if (.not.TmSt(N)) then
      TmStart(:,N)=CrTm
      TmSt(N)=.True.
   else
      print*,'Timer Thread #',N,' is already running'
   end if
end if
end subroutine Timer_Run
!******************
!****f* Mtime.f90/TGet
!
! NAME 
!  TGet
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Gives the time measured by timer N.
!
! SOURCE
! 
function TGet(N)
integer,dimension(5)::TGet
integer::N
  if (TmSt(N)) then
     call Date_And_Time(Values=TmStop(:,N))
     TGet=Diff(TmStart(:,N),TmStop(:,N))
  else
     print*,'Timer Thread #',N,' is not running'
  end if   
end function TGet
!******************
!****f* Mtime.f90/TmsGet
!
! NAME 
!  TmsGet
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Gives the time measured by timer N in milliseconds.
!
! SOURCE
! 
function TmsGet(N)
integer::TmsGet
integer::N
  if (TmSt(N)) then
     call Date_And_Time(Values=TmStop(:,N))
     TmsGet=Diffms(TmStart(:,N),TmStop(:,N))
  else
     print*,'Timer Thread #',N,' is not running'
  end if   
end function TmsGet
!******************
!****f* Mtime.f90/Diff
!
! NAME 
!  Diff
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Computes the time difference between two dates.
!
! SOURCE
! 
function Diff(T1,T2)
integer,dimension(5)::Diff
integer,dimension(8)::T1,T2
integer::Ndays1,Ndays2
if (T1(1).ne.T2(1)) then
   print*,'Time interval that covers more than one year not implemented yet'
   Diff=0
   Return
end if
Ndays1=Sum(NDinM(1:T1(2)-1))+T1(3)
Ndays2=Sum(NDinM(1:T2(2)-1))+T2(3)
Diff(1)=Ndays2-Ndays1
do i=2,5
   Diff(i)=T2(i+3)-T1(i+3)
end do
do i=0,3
   j1=8-i
   j2=5-i
   if (Diff(j2)<0) then
      Diff(j2)=Diff(j2)+C(j1)
      Diff(j2-1)=Diff(j2-1)-1
   end if
end do
end function Diff
!**************************
!****f* Mtime.f90/Diffms
!
! NAME 
!  Diffms
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Computes the time difference between two dates in 
!  milliseconds.
!
! SOURCE
! 
integer function Diffms(T1,T2)
integer,dimension(8)::T1,T2
integer,dimension(5)::dT
dT=Diff(T1,T2)
jjj=C(8)*(C(7)*(c(6)*(C(5)*dT(1)+dT(2))+dT(3))+dT(4))+dT(5)
Diffms=jjj
end function Diffms
!**************************
!****u* Mtime.f90/SetName
!
! NAME 
!  SetName
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Assigns a name string to timer N.
!
! SOURCE
! 
subroutine Timer_SetName(A,N)
character(*)::A
integer::N
Tmnm(N)%Name=A
end subroutine Timer_SetName
!**************************
!****u* Mtime.f90/DTPrint
!
! NAME 
!  DTPrint
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Printe the full time of timer N.
!
! SOURCE
! 
subroutine DTPrint(N,FlUnt)
integer::N,dT(5),LN
integer,optional::FlUnt
character(100)::FM
character(61)::FM2=',''] = '',I5,'' d; '',I2,'' hr; '',I2,'' mn; '',I2,'' sc; '',I3,'' msc'')'
dT=TGet(N)
LN=Len_Trim(Tmnm(N)%Name)
FM='(''Time ['',A'//Int2Char(LN,2)//FM2
if (Present(FlUnt)) then
   write(UNIT=FlUnt,FMT=FM) Tmnm(N)%Name(1:LN),dT
else
   write(UNIT=6,FMT=FM) Tmnm(N)%Name(1:LN),dT
end if
end subroutine DtPrint
!**************************
!****f* Mtime.f90/DTmsPrint
!
! NAME 
!  DTmsPrint
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Printe the full time of timer N in milliseconds.
!
! SOURCE
! 
subroutine DTmsPrint(N,FlUnt)
integer::N,dT,LN
integer,optional::FlUnt
character(100)::FM
dT=TmsGet(N)
LN=Len_Trim(Tmnm(N)%Name)
FM='(''Time ['',A'//Int2Char(LN,2)//',''] = '',I8,'' msc'')'
if (Present(FlUnt)) then
   write(FlUnt,FMT=FM) Tmnm(N)%Name(1:LN),dT
else
   write(*,FMT=FM) Tmnm(N)%Name(1:LN),dT
end if
end subroutine DTmsPrint
!**************************
!****u* Mtime.f90/Tstop
!
! NAME 
!  Tstop
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Stops the timer N.
!
! SOURCE
! 
subroutine Timer_stop(N)
if (TmSt(n)) then
  call Date_And_Time(Values=TmStop(:,N))
  TmSt(N)=.false.
else
   print*,'Timer Thread #',N,' is not running'
end if
end subroutine  Timer_stop
!**************************
subroutine DTmsShow(N,FlUnt)
!****u* Mtime.f90/DTmsShow
!
! NAME 
!  DTmsShow
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Prints the time of timer N in milliseconds together with the 
!  name assigned by SetName.
!
! SOURCE
! 
integer::N,dT,LN
integer,optional::FlUnt
character(100)::FM
dT=Diffms(TmStart(:,N),TmStop(:,N))
LN=Len_Trim(Tmnm(N)%Name)
FM='(''Time ['',A'//Int2Char(LN,2)//',''] = '',I8,'' msc'')'
if (Present(FlUnt)) then
   write(FlUnt,FMT=FM) Tmnm(N)%Name(1:LN),dT
else
   write(*,FMT=FM) Tmnm(N)%Name(1:LN),dT
end if
end subroutine DTmsShow
!**************************
subroutine Timer_DTShow(N,FlUnt)
!****u* Mtime.f90/DTmsShow
!
! NAME 
!  DTmsShow
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Prints the time of timer N in milliseconds/seconds/minutes
!  together with the name assigned by SetName.
!
! SOURCE
! 
integer::N,dT,LN
integer,optional::FlUnt
character(100)::FM
dT=Diffms(TmStart(:,N),TmStop(:,N))
LN=Len_Trim(Tmnm(N)%Name)
if(dT<=1.0D3) then
   FM='(''Time ['',A'//Int2Char(LN,2)//',''] = '',I8,'' msc'')'
end if
if((dT>1E3).and.(dT<1E6)) then
   dT=nint(dble(dT)/1.0D3)
   FM='(''Time ['',A'//Int2Char(LN,2)//',''] = '',I8,'' sc'')'
end if
if(dT>1E6) then
   dT=nint(dble(dT)/6.0D4)
   FM='(''Time ['',A'//Int2Char(LN,2)//',''] = '',I8,'' min'')'
end if

if (Present(FlUnt)) then
   write(FlUnt,FMT=FM) Tmnm(N)%Name(1:LN),dT
else
   write(*,FMT=FM) Tmnm(N)%Name(1:LN),dT
end if
end subroutine Timer_DTShow


!**************************
!****f* Mtime.f90/Int2Char
!
! NAME 
!  Int2Char
!
! AUTHOR
!  Yaroslav Pavlukh
! 
! DESCRIPTION
!  Converts an integer into a set of characters suitable
!  for printing.
!
! SOURCE
! 
function Int2Char(Iy,fieldw)
integer::Ix,fieldw,Iy
character(fieldw)::Int2Char
character(10)::dig='0123456789'
integer,dimension(10)::A
Ix=Iy
if ((Ix<10**fieldw).and.(fieldw<=10)) then
do i1=1,fieldw
i2=fieldw-i1
A(i1)=int(Ix/10**i2)
Ix=Ix-A(i1)*10**i2
end do
A=A+1
do i1=1,fieldw
Int2Char(i1:i1)=dig(A(i1):A(i1))
end do
else 
Int2Char=repeat(' ',fieldw)
end if
end function Int2Char
!**************************

subroutine Timer_tic(A,N)
character(*)::A
integer::N
call Timer_SetName(A,N)
call Timer_Run(N=N)
end subroutine Timer_tic

subroutine Timer_toc(N)
integer::N
call Timer_stop(N=N)
call Timer_DTShow(N=N)
end subroutine Timer_toc

subroutine PrintTime(tag,dT,FlUnt)
character(len=*),intent(in)::tag
double precision,intent(in)::dT
integer,optional::FlUnt
character(100)::FM
integer :: ndT
LN=Len_Trim(tag)
if(dT<=1.0D0) then
   ndT=nint(1D3*dT)
   FM='(''Time ['',A'//Int2Char(LN,2)//',''] = '',I8,'' msc'')'
end if
if((dT>1D0).and.(dT<1D3)) then
   ndT=nint(dT)
   FM='(''Time ['',A'//Int2Char(LN,2)//',''] = '',I8,'' sc'')'
end if
if(dT>1D3) then
   ndT=nint(dT/1.0D3)
   FM='(''Time ['',A'//Int2Char(LN,2)//',''] = '',I8,'' min'')'
end if

if (Present(FlUnt)) then
   write(FlUnt,FMT=FM) tag(1:LN),ndT
else
   write(*,FMT=FM) tag(1:N),ndT
end if
end subroutine PrintTime

!**************************
end module Mtime
















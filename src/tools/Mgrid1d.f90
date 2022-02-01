module Mgrid1D
!==================================================================
!****m* src/Mgrid1d.f03
!
! NAME
!  Mgrid1d
!
! AUTHOR
!  Michael Schueler
!
! DESCRIPTION
!  Contains the type mesh_t (one-dimensional grid) along with
!  some basic integration routines.  
!
! CONTAINS
!   Mgrid1d.f03/mesh_t
!   Mgrid1d.f03/Simpson
!   Mgrid1d.f03/Trapz
! 
!****  
!==================================================================
  use Mdef,only:dp
  use Mutils,only:savetxt,loadtxt
  implicit none
!------------------------------------------------------------------
  private
  public &
       mesh_t,&
       Simpson,&
       Trapz
!------------------------------------------------------------------
  type mesh_t
     integer :: Npts
     real(dp) :: Xmin,Xmax,dX
     real(dp),allocatable :: Xpts(:) 
   contains 
     procedure :: Init => mesh_Init
     procedure :: GridInit => mesh_GridInit
     procedure :: Clean => mesh_Clean
     procedure :: SaveToFile => mesh_savetofile
     procedure :: ReadFromFile => mesh_readfromfile
  end type mesh_t
!------------------------------------------------------------------
  interface Simpson
     module procedure DSimpson,ZSimpson
  end interface Simpson

  interface Trapz
     module procedure DTrapz,ZTrapz
  end interface Trapz
!------------------------------------------------------------------
contains
!------------------------------------------------------------------
  subroutine mesh_GridInit(Mesh,Npts,Xmin,Xmax)
    integer::Npts
    real(dp)::Xmin,Xmax
    class(mesh_t)::Mesh

    Mesh%XMIN=Xmin
    Mesh%XMAX=Xmax
    Mesh%Npts=Npts

    call Mesh%Init
    
  end subroutine Mesh_GridInit
!------------------------------------------------------------------  
  subroutine mesh_Init(Mesh)
    class(mesh_t) :: Mesh
    integer :: i

    if(.not.allocated(Mesh%Xpts)) allocate(Mesh%Xpts(0:Mesh%Npts))
    Mesh%dX=(Mesh%Xmax-Mesh%Xmin)/Mesh%Npts
    do i=0,Mesh%Npts
       Mesh%Xpts(i)=Mesh%Xmin+Mesh%dX*i
    end do

  end subroutine Mesh_Init
!------------------------------------------------------------------
  subroutine mesh_Clean(Mesh)
    class(mesh_t) :: Mesh

    if(allocated(Mesh%Xpts)) deallocate(Mesh%Xpts)
    
  end subroutine Mesh_Clean
!------------------------------------------------------------------
  subroutine mesh_savetofile(mesh,Flname)
    class(mesh_t),intent(in) :: mesh
    character(len=*),intent(in) :: Flname
    real(dp),allocatable :: d(:,:)

    allocate(d(mesh%npts+1,1))
    d(1:mesh%npts+1,1) = mesh%Xpts(0:mesh%npts)
    
    call savetxt(Flname,d)

    deallocate(d)
    
  end subroutine mesh_savetofile
!------------------------------------------------------------------
  subroutine mesh_readfromfile(mesh,Flname)
    real(dp),parameter :: eps10=1.0e-10_dp
    class(mesh_t),intent(inout) :: mesh
    character(len=*),intent(in) :: Flname
    real(dp),allocatable :: d(:,:)

    call loadtxt(Flname,d)

    mesh%Npts = size(d) - 1
    mesh%Xmin = minval(d(:,1))
    mesh%Xmax = maxval(d(:,1))
    if(.not.allocated(mesh%Xpts)) allocate(mesh%Xpts(0:mesh%Npts))
    mesh%Xpts(0:mesh%npts) = d(1:mesh%npts+1,1)
    
    deallocate(d)
    
  end subroutine mesh_readfromfile
!------------------------------------------------------------------
  pure function DSimpson(Mesh,y,Nmax) result(Dint) 
    type(mesh_t),intent(in) :: Mesh
    integer,intent(in),optional :: Nmax
    real(dp),intent(in) :: y(0:)
    real(dp) :: Dint
    !...............
    integer :: m,n,j
    real(dp) :: h
    real(dp), parameter :: c1=2D0/3D0, c2=4D0/3D0
    !...............

    if(present(Nmax)) then
       n=Nmax
    else
       n=Mesh%Npts
    end if

    m = n
    if((mod(n,2)>0).or.(mod(n,2)<0)) m = n-1      
    h = Mesh%dX
    Dint = h/3 * (y(0) + y(m))
    do j=1,m-1
       if(mod(j,2)==0) then
          Dint = Dint + h*c1*y(j)
       else
          Dint = Dint + h*c2*y(j)
       end if
    end do

  end function DSimpson
!------------------------------------------------------------------
  pure function ZSimpson(Mesh,y,Nmax) result(Zint) 
    type(mesh_t),intent(in) :: Mesh
    integer,intent(in), optional :: Nmax
    complex(dp),intent(in) :: y(0:)
    complex(dp) :: Zint
    !...............
    integer :: m,n,j
    real(dp) :: h
    real(dp), parameter :: c1=2D0/3D0, c2=4D0/3D0
    !...............

    if(present(Nmax)) then
       n=Nmax
    else
       n=Mesh%Npts
    end if

    m = n
    if((mod(n,2)>0).or.(mod(n,2)<0)) m = n-1      
    h = Mesh%dX 
    Zint = h/3 * (y(0) + y(m))
    do j=1,m-1
       if(mod(j,2)==0) then
          Zint = Zint + h*c1*y(j)
       else
          Zint = Zint + h*c2*y(j)
       end if
    end do

  end function ZSimpson
!------------------------------------------------------------------
  function DTrapz(Mesh,y,Nmax) result(Dint)
    type(mesh_t) :: Mesh
    integer, optional :: Nmax
    real(dp) :: y(0:),Dint
    !...............
    integer :: n,i
    real(dp) :: h
    !...............
    
    h=Mesh%dX

    if(present(Nmax)) then
       n=min(Nmax,Mesh%Npts)
    else
       n=Mesh%Npts
    end if
    
    if(n==0) then
       Dint=h*y(0)
       return
    else
       Dint=0.5D0*h*(y(0)+y(n))
       do i=1,n-1
          Dint=Dint+h*y(i)
       end do
    end if

  end function DTrapz
!------------------------------------------------------------------
  function ZTrapz(Mesh,y,Nmax) result(Zint)
    type(mesh_t) :: Mesh
    integer, optional :: Nmax
    complex(dp) :: y(0:),Zint
    !...............
    integer :: n,i
    real(dp) :: h
    !...............
    
    h=Mesh%dX

    if(present(Nmax)) then
       n=min(Nmax,Mesh%Npts)
    else
       n=Mesh%Npts
    end if
    
    if(n==0) then
       Zint=h*y(0)
       return
    else
       Zint=0.5D0*h*(y(0)+y(n))
       do i=1,n-1
          Zint=Zint+h*y(i)
       end do
    end if

  end function ZTrapz
!------------------------------------------------------------------

!==================================================================
end module Mgrid1D

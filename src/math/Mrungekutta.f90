module Mrungekutta
!! 5th-order Runga-Kutta-Fehlbert integrator
!======================================================================================
  use Mdef,only: dp
  implicit none
!--------------------------------------------------------------------------------------
  private
  public &
       ODE_step_rk5
!--------------------------------------------------------------------------------------
  real(dp) :: ARK5(5),BRK5(6),CRK5(5,5)  
!--------------------------------------------------------------------------------------
  data ARK5(1),ARK5(2),ARK5(3),ARK5(4),ARK5(5) / &
       0.25D0,0.375D0,0.92307692307692307692D0,1.0D0,0.5D0 /

  data BRK5(1),BRK5(2),BRK5(3),BRK5(4),BRK5(5),BRK5(6) / &
       0.11851851851851851851D0,&
       0.0D0,&
       0.51898635477582846003D0,&
       0.50613149034201665780D0,&
       -0.18D0,&
       0.03636363636363636363D0 /
       
  data CRK5(1,1),CRK5(1,2),CRK5(1,3),CRK5(1,4),CRK5(1,5) / &
       0.25D0,0.0D0,0.0D0,0.0D0,0.0D0 /
  data CRK5(2,1),CRK5(2,2),CRK5(2,3),CRK5(2,4),CRK5(2,5) / &
       0.09375D0,0.28125D0,0.0D0,0.0D0,0.0D0 /
  data CRK5(3,1),CRK5(3,2),CRK5(3,3),CRK5(3,4),CRK5(3,5) / &
       0.87938097405553026854D0,&
       -3.27719617660446062812D0,&
       3.32089212562585343650D0,&
       0.0D0,0.0D0 /
  data CRK5(4,1),CRK5(4,2),CRK5(4,3),CRK5(4,4),CRK5(4,5) / &
       2.03240740740740740740D0,&
       -8.0D0,&
       7.17348927875243664717D0,&
       -0.20589668615984405458D0,&
       0.0D0 /
  data CRK5(5,1),CRK5(5,2),CRK5(5,3),CRK5(5,4),CRK5(5,5) / &
       -0.29629629629629629629D0,&
       2.0D0,&
       -1.38167641325536062378D0,&
       0.45297270955165692007D0,&
       -0.275D0 /
!--------------------------------------------------------------------------------------
  abstract interface
     function deriv_func_dscalar(t,y) result(dydt)
       import :: dp
       real(dp),intent(in) :: t
       real(dp),intent(in) :: y
       real(dp)            :: dydt
     end function deriv_func_dscalar

     function deriv_func_zscalar(t,y) result(dydt)
       import :: dp
       real(dp),intent(in)    :: t
       complex(dp),intent(in) :: y
       complex(dp)            :: dydt
     end function deriv_func_zscalar

     function deriv_func_dvector(size,t,y) result(dydt)
       import :: dp
       integer,intent(in)  :: size
       real(dp),intent(in) :: t
       real(dp),intent(in) :: y(:)
       real(dp)            :: dydt(size)
     end function deriv_func_dvector

     function deriv_func_zvector(size,t,y) result(dydt)
       import :: dp
       integer,intent(in)     :: size
       real(dp),intent(in)    :: t
       complex(dp),intent(in) :: y(:)
       complex(dp)            :: dydt(size)
     end function deriv_func_zvector

     function deriv_func_dmatrix(size,t,y) result(dydt)
       import :: dp
       integer,intent(in)     :: size
       real(dp),intent(in)    :: t
       real(dp),intent(in)    :: y(:,:)
       real(dp)               :: dydt(size,size)
     end function deriv_func_dmatrix
     
     function deriv_func_zmatrix(size,t,y) result(dydt)
       import :: dp
       integer,intent(in)     :: size
       real(dp),intent(in)    :: t
       complex(dp),intent(in) :: y(:,:)
       complex(dp)            :: dydt(size,size)
     end function deriv_func_zmatrix
     
  end interface
!--------------------------------------------------------------------------------------
  interface ODE_step_RK5
  !! generic interface for performing a Runge-Kutta time step.
     module procedure &
          ODE_step_RK5_dscalar,&
          ODE_step_RK5_zscalar,&
          ODE_step_RK5_dvector,&
          ODE_step_RK5_zvector,&
          ODE_step_RK5_dmatrix,&
          ODE_step_RK5_zmatrix
  end interface ODE_step_RK5
!-------------------------------------------------------------------------------------- 
contains
!------------------------------------------------------------------------------------
  function ODE_step_RK5_dscalar(n,dt,deriv_func,xn) result(xn1)
    !! Runge-Kutta step \( x(n \Delta t) \rightarrow x((n+1)\Delta t)\) for a real scalar 
    !! function and right-hand side \( \dot{x}(t) = f(t, x(t)) \).
    integer,intent(in)            :: n !! index of the current time step
    real(dp),intent(in)           :: dt !! time step size \(\Delta t \)
    procedure(deriv_func_dscalar) :: deriv_func !! right-hand side \( f(t,x) \)
    real(dp),intent(in)           :: xn !! current value \(x(n \Delta t)\)
    real(dp)                      :: xn1 !! new value \(x((n+1) \Delta t)\)
    integer  :: k
    real(dp) :: dxk(6),xk

    dxk(1) = deriv_func(n*dt,xn)
    do k=1,5
       xk = xn + dt*sum(CRK5(k,1:k)*dxk(1:k))
       dxk(k+1) = deriv_func((n+ARK5(k))*dt,xk)
    end do

    xn1 = xn + dt*sum(BRK5*dxk)

  end function ODE_step_RK5_dscalar
!--------------------------------------------------------------------------------------  
  function ODE_step_RK5_zscalar(n,dt,deriv_func,xn) result(xn1)
    !! Runge-Kutta step \( x(n \Delta t) \rightarrow x((n+1)\Delta t)\) for a complex scalar 
    !! function and right-hand side \( \dot{x}(t) = f(t, x(t)) \).
    integer,intent(in)            :: n !! index of the current time step
    real(dp),intent(in)           :: dt !! time step size \(\Delta t \)
    procedure(deriv_func_zscalar) :: deriv_func !! right-hand side \( f(t,x) \)
    complex(dp),intent(in)        :: xn !! current value \(x(n \Delta t)\)
    complex(dp)                   :: xn1 !! new value \(x((n+1) \Delta t)\)
    integer  :: k
    complex(dp) :: dxk(6),xk

    dxk(1) = deriv_func(n*dt,xn)
    do k=1,5
       xk = xn + dt*sum(CRK5(k,1:k)*dxk(1:k))
       dxk(k+1) = deriv_func((n+ARK5(k))*dt,xk)
    end do

    xn1 = xn + dt*sum(BRK5*dxk)

  end function ODE_step_RK5_zscalar
!--------------------------------------------------------------------------------------  
  function ODE_step_RK5_dvector(size,n,dt,deriv_func,xn) result(xn1)
    integer,intent(in)            :: size
    integer,intent(in)            :: n
    real(dp),intent(in)           :: dt
    procedure(deriv_func_dvector) :: deriv_func
    real(dp),intent(in)           :: xn(:)
    real(dp)                      :: xn1(size)
    integer  :: k,j
    real(dp) :: dxk(size,6),xk(size)

    dxk(:,1) = deriv_func(size,n*dt,xn)
    do k=1,5
       xk = xn
       do j=1,k
          xk(:) = xk(:) + dt*CRK5(k,j)*dxk(:,j)
       end do
       dxk(:,k+1) = deriv_func(size,(n+ARK5(k))*dt,xk)
    end do

    xn1 = xn
    do j=1,6
       xn1(:) = xn1(:) + dt*BRK5(j)*dxk(:,j)
    end do

  end function ODE_step_RK5_dvector
!------------------------------------------------------------------------------------
  function ODE_step_RK5_zvector(size,n,dt,deriv_func,xn) result(xn1)
    integer,intent(in)            :: size
    integer,intent(in)            :: n
    real(dp),intent(in)           :: dt
    procedure(deriv_func_zvector) :: deriv_func
    complex(dp),intent(in)        :: xn(:)
    complex(dp)                   :: xn1(size)
    integer     :: k,j
    complex(dp) :: dxk(size,6),xk(size)

    dxk(:,1) = deriv_func(size,n*dt,xn)
    do k=1,5
       xk = xn
       do j=1,k
          xk(:) = xk(:) + dt*CRK5(k,j)*dxk(:,j)
       end do
       dxk(:,k+1) = deriv_func(size,(n+ARK5(k))*dt,xk)
    end do

    xn1 = xn
    do j=1,6
       xn1(:) = xn1(:) + dt*BRK5(j)*dxk(:,j)
    end do

  end function ODE_step_RK5_zvector
!------------------------------------------------------------------------------------
  function ODE_step_RK5_dmatrix(size,n,dt,deriv_func,xn) result(xn1)
    integer,intent(in)            :: size
    integer,intent(in)            :: n
    real(dp),intent(in)           :: dt
    procedure(deriv_func_dmatrix) :: deriv_func
    real(dp),intent(in)           :: xn(:,:)
    real(dp)                      :: xn1(size,size)
    integer  :: k,j
    real(dp) :: dxk(size,size,6),xk(size,size)

    dxk(:,:,1) = deriv_func(size,n*dt,xn)
    do k=1,5
       xk = xn(:,:)
       do j=1,k
          xk(:,:) = xk(:,:) + dt*CRK5(k,j)*dxk(:,:,j)
       end do
       dxk(:,:,k+1) = deriv_func(size,(n+ARK5(k))*dt,xk)
    end do

    xn1 = xn(:,:)
    do j=1,6
       xn1(:,:) = xn1(:,:) + dt*BRK5(j)*dxk(:,:,j)
    end do

  end function ODE_step_RK5_dmatrix
!------------------------------------------------------------------------------------
   function ODE_step_RK5_zmatrix(size,n,dt,deriv_func,xn) result(xn1)
    integer,intent(in)            :: size
    integer,intent(in)            :: n
    real(dp),intent(in)           :: dt
    procedure(deriv_func_zmatrix) :: deriv_func
    complex(dp),intent(in)        :: xn(:,:)
    complex(dp)                   :: xn1(size,size)
    integer     :: k,j
    complex(dp) :: dxk(size,size,6),xk(size,size)

    dxk(:,:,1) = deriv_func(size,n*dt,xn)
    do k=1,5
       xk = xn(:,:)
       do j=1,k
          xk(:,:) = xk(:,:) + dt*CRK5(k,j)*dxk(:,:,j)
       end do
       dxk(:,:,k+1) = deriv_func(size,(n+ARK5(k))*dt,xk)
    end do

    xn1 = xn(:,:)
    do j=1,6
       xn1(:,:) = xn1(:,:) + dt*BRK5(j)*dxk(:,:,j)
    end do

  end function ODE_step_RK5_zmatrix
!------------------------------------------------------------------------------------
  
!====================================================================================== 
end module Mrungekutta

module scitools_laserpulse
!======================================================================================
!!  Laserpulse class, constructed in time interval [t1,t2] from 4th-order B-splines 
  use scitools_def,only: dp,iu
  use scitools_bsplines,only: spline1d_t
  use scitools_rungekutta,only: ODE_step_RK5
  use scitools_utils,only: loadtxt
!--------------------------------------------------------------------------------------
  private
  public &
       scalarfunc_spline_t,&
       LaserPulse_spline_t
!--------------------------------------------------------------------------------------
  integer,parameter :: kx=4
!--------------------------------------------------------------------------------------
  type scalarfunc_spline_t
  !! scalar time-depndent function \(y(t)\)
     !.............................................. 
     integer,private                           :: Npts !! number of sample points
     real(dp),public                           :: Tmin,Tmax !! time interval
     real(dp),allocatable,dimension(:),public  :: tpts !! time samples
     type(spline1d_t),private                  :: fspl !! b-spline object
     !..............................................
   contains
     !..............................................
     procedure, public  :: Init => scalarfunc_Init
     procedure, public  :: fval => scalarfunc_fval
     procedure, public  :: Load_function
     !..............................................
  end type scalarfunc_spline_t 

  type LaserPulse_spline_t
  !! Laser pulse class describing electric field \(E(t)\) and vector potential \(A(t)\)
     !..............................................
     logical,private                           :: E_set=.false.,A_set=.false.
     integer,private                           :: Npts
     real(dp),public                           :: Tmin=0.0_dp
     real(dp),public                           :: Tmax=-huge(1.0_dp)
     real(dp),allocatable,dimension(:),public  :: tpts
     type(spline1d_t),private                  :: Aspl
     type(spline1d_t),private                  :: Espl
     !..............................................
   contains
     !..............................................
     procedure, public  :: Init => Laserpulse_Init
     procedure, public  :: Afield => Laserpulse_Afield
     procedure, public  :: Efield => Laserpulse_Efield
     procedure, public  :: Energy => Laserpulse_Energy
     procedure, public  :: Load_vectorpotential
     procedure, public  :: Load_electricfield
     !..............................................
  end type LaserPulse_spline_t 
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
  subroutine scalarfunc_Init(self,Npts,Tmin,Tmax)
  !! Initializes the scalar function
    class(scalarfunc_spline_t),intent(inout)::self
    integer,intent(in)::Npts
    real(dp),intent(in)::Tmin,Tmax
    integer::it,iflag

    self%Npts=Npts
    if(allocated(self%tpts)) deallocate(self%tpts)
    allocate(self%tpts(Npts))
    
    self%Tmin=Tmin
    self%Tmax=Tmax
    forall(it=1:Npts) self%tpts(it)=Tmin+(Tmax-Tmin)*(it-1)/dble(Npts-1)
    
  end subroutine scalarfunc_Init
!--------------------------------------------------------------------------------------
  subroutine Load_function(self,filename,usecol)
  !! Reads the sampled function \( (t_i, y(t_i)) \) from file and 
  !! constructs the b-spline cofficients for interpolation.
    class(scalarfunc_spline_t),intent(inout) :: self
    character(len=*),intent(in)              :: filename !! file name of text file
    integer,intent(in),optional              :: usecol !! which column in the text file to read \(y(t_i)\)
    integer :: it,col,iflag,inbvx
    real(dp),allocatable,dimension(:,:) :: Adata

    col = 2
    if(present(usecol)) col = usecol

    call loadtxt(filename,Adata)
    call self%Init(size(Adata,1),minval(Adata(:,1)),maxval(Adata(:,1)))

    iflag = 0
    call self%fspl%Init(self%tpts,Adata(:,col),kx,iflag)

  end subroutine Load_function
!--------------------------------------------------------------------------------------
  real(dp) function scalarfunc_fval(self,t)
  !! Evaluates \(y(t)\) from the b-spline cofficients. If \(t\) is out of bounds, 
  !! a constant value corresponding to the left or right boundary, respectively, is returned.
    real(dp),parameter :: eps=1.0e-12_dp
    class(scalarfunc_spline_t),intent(in)::self
    real(dp),intent(in)::t
    integer::iflag,inbvx
    real(dp)::EF
    
    inbvx=1

    if(t < self%Tmin) then
      scalarfunc_fval = self%fspl%eval(self%Tmin + eps,0,iflag,inbvx) 
      return
    end if

    if(t > self%Tmax) then
      scalarfunc_fval = self%fspl%eval(self%Tmax - eps,0,iflag,inbvx) 
      return
    end if

    scalarfunc_fval = self%fspl%eval(t,0,iflag,inbvx) 

  end function scalarfunc_fval
!--------------------------------------------------------------------------------------
  subroutine Laserpulse_Init(self,Npts,Tmin,Tmax)
  !! Initializes the laser pulse class and allocates internal storage.
    class(LaserPulse_spline_t),intent(inout)::self
    integer,intent(in)::Npts
    real(dp),intent(in)::Tmin,Tmax
    integer::it,iflag

    self%Npts=Npts
    if(allocated(self%tpts)) deallocate(self%tpts)
    allocate(self%tpts(Npts))
    
    self%Tmin=Tmin
    self%Tmax=Tmax
    forall(it=1:Npts) self%tpts(it)=Tmin+(Tmax-Tmin)*(it-1)/dble(Npts-1)
    
  end subroutine Laserpulse_Init
!--------------------------------------------------------------------------------------   
  real(dp) function Laserpulse_Afield(self,t)
  !! Returns the vector potential \(A(t)\) from the b-spline cofficients.
  !! If \(t < t_\mathrm{min} \) (left boundary), 0 is returned;
  !! If \(t > t_\mathrm{max} \) (right boundary), \(A(t_\mathrm{max})\) is returned.
    class(LaserPulse_spline_t),intent(in)::self
    real(dp),intent(in)::t !! time at which the vector potential is evaluated
    integer::iflag,inbvx
    real(dp)::AF

    if((t>self%Tmin).and.(t<self%Tmax)) then
       inbvx=1
       AF = self%Aspl%eval(t,0,iflag,inbvx) 
       Laserpulse_Afield=AF
    elseif(t <= self%Tmin) then
       Laserpulse_Afield=0.0_dp
    else 
       inbvx=1
       Laserpulse_Afield=self%Aspl%eval(self%Tmax,0,iflag,inbvx) 
    end if
       
  end function Laserpulse_Afield
!--------------------------------------------------------------------------------------
  real(dp) function Laserpulse_Efield(self,t)
  !! Returns the electric field \(E(t)\) from the b-spline cofficients. If \(t\) is 
  !! out of bounds, 0 is returned. 
    class(LaserPulse_spline_t),intent(in)::self
    real(dp),intent(in)::t  !! time at which the field is evaluated
    integer::iflag,inbvx
    real(dp)::EF
    
    if((t>self%Tmin).and.(t<self%Tmax)) then 
       inbvx=1
       EF = self%Espl%eval(t,0,iflag,inbvx) 
       Laserpulse_Efield=EF
    else
       Laserpulse_Efield=0.0_dp
    end if

  end function Laserpulse_Efield
!--------------------------------------------------------------------------------------
  real(dp) function Laserpulse_Energy(self)
  !! Computes the energy \(\int dt\, E(t)^2  \).
    class(LaserPulse_spline_t),intent(inout)::self
    integer,parameter::np=400
    integer::k
    real(dp)::dt,ct

    dt=(self%Tmax-self%Tmin)/dble(self%Npts-1)
    
    Laserpulse_Energy=0.5_dp*dt*(self%Efield(self%tpts(1))**2+self%Efield(self%tpts(self%Npts))**2)
    do k=2,self%Npts-1
       ct=self%Tmin+dt*(k-1)
       Laserpulse_Energy=Laserpulse_Energy+dt*self%Efield(ct)**2
    end do
    
  end function Laserpulse_Energy
!--------------------------------------------------------------------------------------
  subroutine Load_VectorPotential(self,filename,usecol,CalcEfield)
  !! Reads the vector potential from a plain text file. We assume that 
  !! the time sample points \(t_i\) are stored in the first column, while the
  !! sampled values \(A(t_i)\) are stored in another column specified by `usecol`. 
  !! After reading from file, the b-spline coefficients for interpolating are constructed.
  !! The electric field is also computed by \(E(t) = - \dot{A}(t)\) by default.
    class(LaserPulse_spline_t),intent(inout) :: self
    character(len=*),intent(in)              :: filename !! file name for text file (column structure)
    integer,intent(in),optional              :: usecol !! which column to use to read \(A(t_i\); default: `usecol=2`
    logical,intent(in),optional              :: CalcEfield !! Whether to compute the electric field. Default is `.true.`
    logical :: calcefield_
    integer :: it,col,iflag,inbvx
    real(dp),allocatable,dimension(:,:) :: Adata
    real(dp),allocatable :: Edata(:)

    col = 2
    if(present(usecol)) col = usecol

    calcefield_ = .true.
    if(present(CalcEfield)) calcefield_ = CalcEfield
       
    call loadtxt(filename,Adata)
    call self%Init(size(Adata,1),minval(Adata(:,1)),maxval(Adata(:,1)))

    iflag = 0
    call self%Aspl%Init(self%tpts,Adata(:,col),kx,iflag)
    self%A_set = .true.
    
    deallocate(Adata)
    
    if(.not.self%E_set) then
       allocate(Edata(self%Npts)); Edata = 0.0_dp
    
       if(calcefield_) then       
          inbvx = 1
          do it=1,self%Npts
             Edata(it) = -self%Aspl%Eval(self%tpts(it),1,iflag,inbvx)
          end do
       end if
       
       iflag = 0
       call self%Espl%Init(self%tpts,Edata,kx,iflag)
       self%E_set = .true.
       deallocate(Edata)
    end if
             
  end subroutine Load_VectorPotential
!--------------------------------------------------------------------------------------
  subroutine Load_ElectricField(self,filename,usecol,CalcAfield)
  !! Reads the electric field from a plain text file. We assume that 
  !! the time sample points \(t_i\) are stored in the first column, while the
  !! sampled values \(E(t_i)\) are stored in another column specified by `usecol`. 
  !! After reading from file, the b-spline coefficients for interpolating are constructed.
  !! The vector potential is also computed by \(A(t) = - \int^t dt\, E(t) \) by default.
    class(LaserPulse_spline_t),intent(inout) :: self
    character(len=*),intent(in)              :: filename !! file name for text file (column structure)
    integer,intent(in),optional              :: usecol !! which column to use to read \(A(t_i\); default: `usecol=2`
    logical,intent(in),optional              :: CalcAfield !! Whether to compute the vector potential. Default is `.true.`
    logical :: calcafield_
    integer :: it,col,iflag,inbvx
    real(dp) :: dt
    real(dp),allocatable,dimension(:,:) :: Edata
    real(dp),allocatable :: Adata(:)

    col = 2
    if(present(usecol)) col = usecol

    calcafield_ = .true.
    if(present(CalcAfield)) calcafield_ = CalcAfield
       
    call loadtxt(filename,Edata)
    call self%Init(size(Edata,1),minval(Edata(:,1)),maxval(Edata(:,1)))

    iflag = 0
    call self%Espl%Init(self%tpts,Edata(:,col),kx,iflag)
    self%E_set = .true.
    deallocate(Edata)

    if(.not.self%A_set) then
       allocate(Adata(self%Npts)); Adata = 0.0_dp
    
       if(calcafield_) then
          dt = self%tpts(2)-self%tpts(1)
          do it=1,self%Npts-1
             Adata(it+1) = ODE_step_RK5(it-1,dt,deriv_func,Adata(it))
          end do          
       end if
       
       iflag = 0
       call self%Aspl%Init(self%tpts,Adata,kx,iflag)
       self%A_set = .true.
       deallocate(Adata)
    end if
  !===================================
  contains
  !===================================  
    function deriv_func(t,y) result(dydt)
      real(dp),intent(in) :: t
      real(dp),intent(in) :: y
      real(dp) :: dydt
      integer :: iflag,inbvx

      inbvx = 1
      dydt = -self%Espl%Eval(t+self%Tmin,0,iflag,inbvx)
      
    end function deriv_func
  !===================================  
  end subroutine Load_ElectricField
!--------------------------------------------------------------------------------------

  
!======================================================================================    
end module scitools_laserpulse
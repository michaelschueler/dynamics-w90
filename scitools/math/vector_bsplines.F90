module scitools_vector_bsplines
!! Vector and matrix-valued B-splines, based on [[scitools_bsplines]]
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use scitools_debug
   use scitools_def,only: dp,iu,zero,one
   use scitools_bsplines,only: db1ink, db1val
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: &
      real_vector_spline_t, &
      real_matrix_spline_t, &
      cplx_vector_spline_t, &
      cplx_matrix_spline_t
!-------------------------------------------------------------------------------------- 
   type :: real_vector_spline_t
   !! class for real vector-valued splines
      integer :: kx=4
      integer :: nc,nx
      real(dp) :: xlim(2)
      real(dp),allocatable,dimension(:)   :: tx
      real(dp),allocatable,dimension(:,:) :: cff
   contains
      procedure :: Init => real_vector_init
      procedure :: Clean => real_vector_clean
      procedure :: Eval => real_vector_Eval    
      procedure :: Eval_component => real_vector_Eval_comp
   end type real_vector_spline_t

   type :: real_matrix_spline_t
   !! class for real matrix-valued splines
      integer :: kx=4
      integer :: nc,mc,nx
      real(dp) :: xlim(2)
      real(dp),allocatable,dimension(:)     :: tx
      real(dp),allocatable,dimension(:,:,:) :: cff
   contains
      procedure :: Init => real_matrix_init
      procedure :: Clean => real_matrix_clean
      procedure :: Eval => real_matrix_Eval  
      procedure :: Eval_component => real_matrix_Eval_comp     
   end type real_matrix_spline_t

   type :: cplx_vector_spline_t
   !! class for complex vector-valued splines
      integer :: kx=4
      integer :: nc,nx
      real(dp) :: xlim(2)
      real(dp),allocatable,dimension(:)   :: tx
      real(dp),allocatable,dimension(:,:) :: cffr,cffi
   contains
      procedure :: Init => cplx_vector_init
      procedure :: Clean => cplx_vector_clean
      procedure :: Eval => cplx_vector_Eval   
      procedure :: Eval_component => cplx_vector_Eval_comp    
   end type cplx_vector_spline_t

   type :: cplx_matrix_spline_t
   !! class for complex matrix-valued splines
      integer :: kx=4
      integer :: nc,mc,nx
      real(dp) :: xlim(2)
      real(dp),allocatable,dimension(:)     :: tx
      real(dp),allocatable,dimension(:,:,:) :: cffr,cffi
   contains
      procedure :: Init => cplx_matrix_init
      procedure :: Clean => cplx_matrix_clean
      procedure :: Eval => cplx_matrix_Eval       
      procedure :: Eval_component => cplx_matrix_Eval_comp
   end type cplx_matrix_spline_t
!-------------------------------------------------------------------------------------- 
contains
!--------------------------------------------------------------------------------------   
   subroutine real_vector_init(me,x,y,nc,kx) 
   !! initializes a real vector-valued B-spline from sampled values \( (x_i, y_i) \)
      class(real_vector_spline_t) :: me
      real(dp),intent(in) :: x(:) !! sample points \(x_i\), `nx` values
      real(dp),intent(in) :: y(:,:) !! sampled values \(y_i\), dimension `(nx,nc)` 
      integer,intent(in)  :: nc !! vector size of \(y(x)\)
      integer,intent(in),optional :: kx !! order of B-spline interpolant
      integer :: i,iflag

      if(present(kx)) me%kx = kx

      me%nx = size(x)
      me%nc = nc
      me%xlim = [minval(x), maxval(x)]

      if(.not.allocated(me%tx)) allocate(me%tx(me%nx + me%kx))
      if(.not.allocated(me%Cff)) allocate(me%Cff(me%nx,me%nc))      

      iflag = 0
      call db1ink(x,me%nx,y(:,1),me%kx,me%tx,me%Cff(:,1),iflag)

      do i=2,me%nc
         iflag = 1
         call db1ink(x,me%nx,y(:,i),me%kx,me%tx,me%Cff(:,i),iflag)
      end do

   end subroutine real_vector_init
!--------------------------------------------------------------------------------------
   subroutine real_vector_clean(me)
      class(real_vector_spline_t) :: me

      if(allocated(me%tx)) deallocate(me%tx)
      if(allocated(me%cff)) deallocate(me%cff)

   end subroutine real_vector_clean
!--------------------------------------------------------------------------------------
   function real_vector_Eval(me,x,idx,iflag) result(y)
   !! Evaluates the B-spline interpolant \(y(x)\) or a derivative. Returns the whole vector.
      class(real_vector_spline_t) :: me
      real(dp),intent(in) :: x !! point where the interpolant is evaluated
      integer,intent(in),optional  :: idx !! order of derivative
      integer,intent(out),optional :: iflag !! info flag
      real(dp) :: y(me%nc)
      integer :: idx_,iflag_
      integer :: i
      integer :: inbvx

      idx_ = 0
      if(present(idx)) idx_ = idx

      if(x < me%xlim(1) .or. x > me%xlim(2)) then
         write(output_unit,fmt700) "real_vector_Eval: out of bounds"
         y = 0.0_dp
         return
      end if       

      do i=1,me%nc
         inbvx = 1
         call db1val(x,idx_,me%tx,me%nx,me%kx,me%cff(:,i),y(i),iflag_,inbvx)
      end do

      if(present(iflag)) iflag = iflag_

   end function real_vector_Eval
!--------------------------------------------------------------------------------------
   function real_vector_Eval_comp(me,x,ic,idx,iflag,inbvx) result(y)
   !! evaluates a component of the vector \(y(x)\) or a derivative
      class(real_vector_spline_t) :: me
      real(dp),intent(in) :: x !! point where the interpolant is evaluated
      integer,intent(in)  :: ic !! index of vector component
      integer,intent(in),optional  :: idx !! order of derivative
      integer,intent(out),optional :: iflag !! info flag
      integer,intent(inout),optional :: inbvx !! initialization parameter which must be set to 
                                              !! 1 the first time this routine is called, and 
                                              !! must not be changed by the user.
      real(dp) :: y
      integer :: idx_,iflag_
      integer :: inbvx_

      idx_ = 0
      if(present(idx)) idx_ = idx

      if(x < me%xlim(1) .or. x > me%xlim(2)) then
         write(output_unit,fmt700) "real_vector_Eval: out of bounds"
         y = 0.0_dp
         return
      end if       

      inbvx_ = 1
      if(present(inbvx)) inbvx_ = inbvx
      call db1val(x,idx_,me%tx,me%nx,me%kx,me%cff(:,ic),y,iflag_,inbvx_)

      if(present(iflag)) iflag = iflag_
      if(present(inbvx)) inbvx = inbvx_

   end function real_vector_Eval_comp
!--------------------------------------------------------------------------------------
   subroutine real_matrix_init(me,x,y,nc,mc,kx)
   !! initializes a real matrix-valued B-spline from sampled values \( (x_i, y_i) \)
      class(real_matrix_spline_t) :: me
      real(dp),intent(in) :: x(:) !! sample points \(x_i\), dimension `nx`
      real(dp),intent(in) :: y(:,:,:) !! sampled values \(y_i\), dimension `(nx,nc,mc)` 
      integer,intent(in)  :: nc,mc !! matrix ranks 
      integer,intent(in),optional :: kx !! order of B-spline interpolant
      integer :: i,j,iflag

      if(present(kx)) me%kx = kx

      me%nx = size(x)
      me%nc = nc
      me%mc = mc
      me%xlim = [minval(x), maxval(x)]

      if(.not.allocated(me%tx)) allocate(me%tx(me%nx + me%kx))
      if(.not.allocated(me%Cff)) allocate(me%Cff(me%nx,me%nc,me%mc))      

      iflag = 0
      call db1ink(x,me%nx,y(:,1,1),me%kx,me%tx,me%Cff(:,1,1),iflag)

      do i=1,me%nc
         do j=1,me%mc
            iflag = 1
            call db1ink(x,me%nx,y(:,i,j),me%kx,me%tx,me%Cff(:,i,j),iflag)
         end do
      end do

   end subroutine real_matrix_init
!--------------------------------------------------------------------------------------
   subroutine real_matrix_clean(me)
      class(real_matrix_spline_t) :: me

      if(allocated(me%tx)) deallocate(me%tx)
      if(allocated(me%cff)) deallocate(me%cff)

   end subroutine real_matrix_clean
!--------------------------------------------------------------------------------------
   function real_matrix_Eval(me,x,idx,iflag) result(y)
   !! Evaluates the B-spline interpolant \(y(x)\) or a derivative. Returns the whole matrix.
      class(real_matrix_spline_t) :: me
      real(dp),intent(in) :: x !! point where the interpolant is evaluated
      integer,intent(in),optional  :: idx !! order of derivative
      integer,intent(out),optional :: iflag !! info flag
      real(dp) :: y(me%nc,me%mc)
      integer :: idx_,iflag_
      integer :: i,j
      integer :: inbvx

      idx_ = 0
      if(present(idx)) idx_ = idx

      if(x < me%xlim(1) .or. x > me%xlim(2)) then
         write(output_unit,fmt700) "real_matrix_Eval: out of bounds"
         y = 0.0_dp
         return
      end if       

      do i=1,me%nc
         do j=1,me%mc
            inbvx = 1
            call db1val(x,idx_,me%tx,me%nx,me%kx,me%cff(:,i,j),y(i,j),iflag_,inbvx)
         end do
      end do

      if(present(iflag)) iflag = iflag_

   end function real_matrix_Eval
!--------------------------------------------------------------------------------------
   function real_matrix_Eval_comp(me,x,ic,jc,idx,iflag,inbvx) result(y)
   !! evaluates an element of the matrix \(y(x)\) or a derivative
      class(real_matrix_spline_t) :: me
      real(dp),intent(in) :: x !! point where the interpolant is evaluated
      integer,intent(in)  :: ic,jc !! indices of matrix elements
      integer,intent(in),optional  :: idx !! order of derivative
      integer,intent(out),optional :: iflag !! info flag
      integer,intent(inout),optional :: inbvx !! initialization parameter which must be set to 
                                              !! 1 the first time this routine is called, and 
                                              !! must not be changed by the user.
      real(dp) :: y
      integer :: idx_,iflag_
      integer :: inbvx_

      idx_ = 0
      if(present(idx)) idx_ = idx

      if(x < me%xlim(1) .or. x > me%xlim(2)) then
         write(output_unit,fmt700) "real_matrix_Eval: out of bounds"
         y = 0.0_dp
         return
      end if         

      inbvx_ = 1
      if(present(inbvx)) inbvx_ = inbvx
      call db1val(x,idx_,me%tx,me%nx,me%kx,me%cff(:,ic,jc),y,iflag_,inbvx_)

      if(present(iflag)) iflag = iflag_
      if(present(inbvx)) inbvx = inbvx_

   end function real_matrix_Eval_comp
!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------   
   subroutine cplx_vector_init(me,x,y,nc,kx)
   !! initializes a complex vector-valued B-spline from sampled values \( (x_i, y_i) \)
      class(cplx_vector_spline_t) :: me
      real(dp),intent(in) :: x(:) !! sample points \(x_i\), dimension `nx`
      complex(dp),intent(in) :: y(:,:) !! sampled values \(y_i\), dimension `(nx,nc)` 
      integer,intent(in)  :: nc !! vector size of \(y(x)\)
      integer,intent(in),optional :: kx !! order of B-spline interpolant
      integer :: i,iflag

      if(present(kx)) me%kx = kx

      me%nx = size(x)
      me%nc = nc
      me%xlim = [minval(x), maxval(x)]

      if(.not.allocated(me%tx)) allocate(me%tx(me%nx + me%kx))
      if(.not.allocated(me%Cffr)) allocate(me%Cffr(me%nx,me%nc))     
      if(.not.allocated(me%Cffi)) allocate(me%Cffi(me%nx,me%nc))  

      iflag = 0
      call db1ink(x,me%nx,dble(y(:,1)),me%kx,me%tx,me%Cffr(:,1),iflag)

      do i=i,me%nc
         iflag = 1
         call db1ink(x,me%nx,dble(y(:,i)),me%kx,me%tx,me%Cffr(:,i),iflag)
         call db1ink(x,me%nx,aimag(y(:,i)),me%kx,me%tx,me%Cffi(:,i),iflag)
      end do

   end subroutine cplx_vector_init
!--------------------------------------------------------------------------------------
   subroutine cplx_vector_clean(me)
      class(cplx_vector_spline_t) :: me

      if(allocated(me%tx)) deallocate(me%tx)
      if(allocated(me%cffr)) deallocate(me%cffr)
      if(allocated(me%cffi)) deallocate(me%cffi)

   end subroutine cplx_vector_clean
!--------------------------------------------------------------------------------------
   function cplx_vector_Eval(me,x,idx,iflag) result(y)
   !! Evaluates the B-spline interpolant \(y(x)\) or a derivative. Returns the whole vector.
      class(cplx_vector_spline_t) :: me
      real(dp),intent(in) :: x !! point where the interpolant is evaluated
      integer,intent(in),optional  :: idx !! order of derivative
      integer,intent(out),optional :: iflag !! info flag
      complex(dp) :: y(me%nc)
      integer :: idx_,iflag_
      integer :: i
      integer :: inbvx
      real(dp) :: yr,yi

      idx_ = 0
      if(present(idx)) idx_ = idx

      if(x < me%xlim(1) .or. x > me%xlim(2)) then
         write(output_unit,fmt700) "cplx_vector_Eval: out of bounds"
         y = zero
         return
      end if         

      do i=1,me%nc
         inbvx = 1
         call db1val(x,idx_,me%tx,me%nx,me%kx,me%cffr(:,i),yr,iflag_,inbvx)
         inbvx = 1
         call db1val(x,idx_,me%tx,me%nx,me%kx,me%cffi(:,i),yi,iflag_,inbvx)
         y(i) = cmplx(yr, yi, kind=dp)
      end do

      if(present(iflag)) iflag = iflag_

   end function cplx_vector_Eval
!--------------------------------------------------------------------------------------
   function cplx_vector_Eval_comp(me,x,ic,idx,iflag,inbvx) result(y)
   !! evaluates a component of the vector \(y(x)\) or a derivative
      class(cplx_vector_spline_t) :: me
      real(dp),intent(in) :: x !! point where the interpolant is evaluated
      integer,intent(in)  :: ic !! index of vector component
      integer,intent(in),optional  :: idx !! order of derivative
      integer,intent(out),optional :: iflag !! info flag
      integer,intent(inout),optional :: inbvx(2) !! initialization parameter which must be set to 
                                              !! 1 the first time this routine is called, and 
                                              !! must not be changed by the user.
      complex(dp) :: y
      integer :: idx_,iflag_
      integer :: inbvx_(2)
      real(dp) :: yr,yi

      idx_ = 0
      if(present(idx)) idx_ = idx

      if(x < me%xlim(1) .or. x > me%xlim(2)) then
         write(output_unit,fmt700) "cplx_vector_Eval: out of bounds"
         y = zero
         return
      end if         

      inbvx_ = 1
      if(present(inbvx)) inbvx_ = inbvx
      call db1val(x,idx_,me%tx,me%nx,me%kx,me%cffr(:,ic),yr,iflag_,inbvx_(1))
      call db1val(x,idx_,me%tx,me%nx,me%kx,me%cffi(:,ic),yi,iflag_,inbvx_(2))
      y = cmplx(yr, yi, kind=dp)

      if(present(iflag)) iflag = iflag_
      if(present(inbvx)) inbvx = inbvx_

   end function cplx_vector_Eval_comp
!--------------------------------------------------------------------------------------
   subroutine cplx_matrix_init(me,x,y,nc,mc,kx)
   !! initializes a complex matrix-valued B-spline from sampled values \( (x_i, y_i) \)
      class(cplx_matrix_spline_t) :: me
      real(dp),intent(in) :: x(:) !! sample points \(x_i\), dimension `nx`
      complex(dp),intent(in) :: y(:,:,:) !! sampled values \(y_i\), dimension `(nx,nc,mc)` 
      integer,intent(in)  :: nc,mc !! matrix ranks 
      integer,intent(in),optional :: kx !! order of B-spline interpolant
      integer :: i,j,iflag

      if(present(kx)) me%kx = kx

      me%nx = size(x)
      me%nc = nc
      me%mc = mc
      me%xlim = [minval(x), maxval(x)]

      if(.not.allocated(me%tx)) allocate(me%tx(me%nx + me%kx))
      if(.not.allocated(me%Cffr)) allocate(me%Cffr(me%nx,me%nc,me%mc))    
      if(.not.allocated(me%Cffi)) allocate(me%Cffi(me%nx,me%nc,me%mc))     

      iflag = 0
      call db1ink(x,me%nx,dble(y(:,1,1)),me%kx,me%tx,me%Cffr(:,1,1),iflag)

      do i=1,me%nc
         do j=1,me%mc
            iflag = 1
            call db1ink(x,me%nx,dble(y(:,i,j)),me%kx,me%tx,me%Cffr(:,i,j),iflag)
            iflag = 1
            call db1ink(x,me%nx,aimag(y(:,i,j)),me%kx,me%tx,me%Cffi(:,i,j),iflag)
         end do
      end do

   end subroutine cplx_matrix_init
!--------------------------------------------------------------------------------------
   subroutine cplx_matrix_clean(me)
      class(cplx_matrix_spline_t) :: me

      if(allocated(me%tx)) deallocate(me%tx)
      if(allocated(me%cffr)) deallocate(me%cffr)
      if(allocated(me%cffi)) deallocate(me%cffi)

   end subroutine cplx_matrix_clean
!--------------------------------------------------------------------------------------
   function cplx_matrix_Eval(me,x,idx,iflag) result(y)
   !! Evaluates the B-spline interpolant \(y(x)\) or a derivative. Returns the whole matrix.
      class(cplx_matrix_spline_t) :: me
      real(dp),intent(in) :: x !! point where the interpolant is evaluated
      integer,intent(in),optional  :: idx !! order of derivative
      integer,intent(out),optional :: iflag !! info flag
      complex(dp) :: y(me%nc,me%mc)
      integer :: idx_,iflag_
      integer :: i,j
      integer :: inbvx
      real(dp) :: yr,yi

      idx_ = 0
      if(present(idx)) idx_ = idx

      if(x < me%xlim(1) .or. x > me%xlim(2)) then
         write(output_unit,fmt700) "cplx_matrix_Eval: out of bounds"
         y = zero
         return
      end if         

      do i=1,me%nc
         do j=1,me%mc
            inbvx = 1
            call db1val(x,idx_,me%tx,me%nx,me%kx,me%cffr(:,i,j),yr,iflag_,inbvx)
            inbvx = 1
            call db1val(x,idx_,me%tx,me%nx,me%kx,me%cffi(:,i,j),yi,iflag_,inbvx)
            y(i,j) = cmplx(yr, yi, kind=dp)
         end do
      end do

      if(present(iflag)) iflag = iflag_

   end function cplx_matrix_Eval
!--------------------------------------------------------------------------------------
   function cplx_matrix_Eval_comp(me,x,ic,jc,idx,iflag,inbvx) result(y)
   !! evaluates an element of the matrix \(y(x)\) or a derivative
      class(cplx_matrix_spline_t) :: me
      real(dp),intent(in) :: x !! point where the interpolant is evaluated
      integer,intent(in)  :: ic,jc !! indices of matrix elements
      integer,intent(in),optional  :: idx !! order of derivative
      integer,intent(out),optional :: iflag !! info flag
      integer,intent(inout),optional :: inbvx(2) !! initialization parameter which must be set to 
                                              !! 1 the first time this routine is called, and 
                                              !! must not be changed by the user.
      complex(dp) :: y
      integer :: idx_,iflag_
      integer :: inbvx_(2)
      real(dp) :: yr,yi

      idx_ = 0
      if(present(idx)) idx_ = idx

      if(x < me%xlim(1) .or. x > me%xlim(2)) then
         write(output_unit,fmt700) "cplx_matrix_Eval: out of bounds"
         y = zero
         return
      end if         

      inbvx_ = 1
      if(present(inbvx)) inbvx_ = inbvx

      call db1val(x,idx_,me%tx,me%nx,me%kx,me%cffr(:,ic,jc),yr,iflag_,inbvx_(1))
      call db1val(x,idx_,me%tx,me%nx,me%kx,me%cffi(:,ic,jc),yi,iflag_,inbvx_(2))
      y = cmplx(yr, yi, kind=dp)

      if(present(iflag)) iflag = iflag_
      if(present(inbvx)) inbvx = inbvx_

   end function cplx_matrix_Eval_comp
!--------------------------------------------------------------------------------------

!======================================================================================
end module scitools_vector_bsplines
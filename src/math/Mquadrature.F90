module Mquadrature
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit, error_unit
   use Mdebug
   use Mdef,only: dp
   use Mquadpack,only: dqags,dqagi
   use Mgausslegendre,only: integration_class_1d
   implicit none
   include "../formats.h"
!--------------------------------------------------------------------------------------
   private
   public :: integral_1d, integral_inf_1d
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   subroutine integral_1d(f,a,b,tol,res,meth,neval,abserr)
      real(dp),intent(in)           :: a,b
      real(dp),intent(in)           :: tol
      real(dp),intent(out)          :: res
      integer,intent(in),optional   :: meth
      integer,intent(out),optional  :: neval
      real(dp),intent(out),optional :: abserr
      integer :: meth_
      integer :: ier,neval_
      real(dp) :: abserr_
      type(integration_class_1d) :: my_int
      interface
          real(dp) function f(x)
          use Mdef, only: dp
          implicit none
          real(dp), intent(in) :: x
          end function
      end interface

      meth_ = 8 
      if(present(meth)) meth_ = meth
      
#ifdef MODERNQUAD
      call my_int%initialize(fx=wrapper, xl=a, xu=b, tolx=tol, methodx=meth_)
      call my_int%integrate(res, ier, abserr_)
      call quadrature_err("integration_class_1d", ier, tol, abserr_)
      neval_ = 0
#else
      call dqags(f,a,b,tol,tol,res,abserr_,neval_,ier)
      call quadpack_err("dqags", ier)
#endif

      if(present(neval)) neval = neval_
      if(present(abserr)) abserr = abserr_

   !.................................................
   contains
   !.................................................
   function wrapper(me,x) result(w)
      class(integration_class_1d),intent(inout) :: me
      real(dp),intent(in) :: x
      real(dp) :: w

      w = f(x)
   end function wrapper
   !.................................................

   end subroutine integral_1d
!--------------------------------------------------------------------------------------
   subroutine integral_inf_1d(f,bound,inf,tol,res,neval,abserr)
      real(dp),intent(in)           :: bound
      integer,intent(in)            :: inf
      real(dp),intent(in)           :: tol
      real(dp),intent(out)          :: res
      integer,intent(out),optional  :: neval
      real(dp),intent(out),optional :: abserr
      integer :: meth_
      integer :: ier,neval_
      real(dp) :: abserr_
      interface
          real(dp) function f(x)
          use Mdef, only: dp
          implicit none
          real(dp), intent(in) :: x
          end function
      end interface


      call dqagi(f,bound,inf,tol,tol,res,abserr_,neval_,ier)
      call quadpack_err("dqagi", ier)

      if(present(neval)) neval = neval_
      if(present(abserr)) abserr = abserr_

   end subroutine integral_inf_1d
!--------------------------------------------------------------------------------------
   subroutine quadpack_err(tag,ier)   
      use Mutils,only: str
      character(len=*),intent(in) :: tag
      integer,intent(in) :: ier

      if(ier /= 0) then
         write(error_unit,fmt901) trim(tag) // " exit code = "//str(ier)
         stop
      end if

   end subroutine quadpack_err
!--------------------------------------------------------------------------------------
   subroutine quadpack_info(tag,abserr,neval)   
      use Mutils,only: str
      character(len=*),intent(in) :: tag
      real(dp),intent(in) :: abserr
      integer,intent(in) :: neval

      write(output_unit,fmt602) trim(tag)//": err = ",abserr
      write(output_unit,fmt600) trim(tag)//": neval = "//str(neval)

   end subroutine quadpack_info
!--------------------------------------------------------------------------------------
   subroutine quadrature_err(tag,ier,tol,err)   
      use Mutils,only: str
      character(len=*),intent(in) :: tag
      integer,intent(in) :: ier
      real(dp),intent(in) :: tol,err

      if(ier == -1) then
         write(error_unit,fmt901) trim(tag) // " exit code = -1"
         stop
      end if

      if((ier == 2) .and. (abs(err) > tol)) then
         write(output_unit,fmt700) trim(tag) // " Error to large!"
      end if

   end subroutine quadrature_err
!--------------------------------------------------------------------------------------

!======================================================================================
end module Mquadrature
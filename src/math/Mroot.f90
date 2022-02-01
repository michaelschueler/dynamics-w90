module Mroot

  ! Optimization algorithms

  use Mdef, only: dp
  use Mutils, only: stop_error
  implicit none
  private
  public bisect,secant,newton,brent

  interface
     real(dp) function func(x)
       import :: dp
       implicit none
       real(dp), intent(in) :: x
     end function func
  end interface

contains

  real(dp) function bisect(f, a, b, tol) result(c)
    ! Solves f(x) = 0 on the interval [a, b] using the bisection method
    procedure(func) :: f
    real(dp), intent(in) :: a, b, tol
    real(dp) :: a_, b_, fa, fb, fc
    a_ = a; b_ = b
    fa = f(a_)
    fb = f(b_)
    if (fa * fb >= 0) then
       call stop_error("bisect: f(a) and f(b) must have opposite signs")
    end if
    do while (b_ - a_ > tol)
       c = (a_ + b_) / 2
       fc = f(c)
       if (abs(fc) < tiny(1.0_dp)) return   ! We need to make sure f(c) is not zero below
       if (fa * fc < 0) then
          b_ = c
          fb = fc
       else
          a_ = c
          fa = fc
       end if
    end do
    c = (a_ + b_)/2
  end function bisect

  real(dp) function secant(f, a,b,tol, maxiter) result (c)
    !Solves f(x) = 0 on using the a and b as starting values
    procedure(func) :: f
    real(dp),intent(in) :: a,b,tol
    integer,optional :: maxiter

    integer :: maxiter_
    real(dp) :: xstart, xnext
    real(dp) :: fstart, fnext
    integer ::  i

    if(present(maxiter)) then
       maxiter_ = maxiter
    else
       maxiter_ = 100
    endif


    xstart = a
    xnext = b
    fstart = f(xstart)
    if(abs(fstart ) < tiny(1.0_dp)) then
       c = xstart
       return
    endif
    fnext = f(xnext)
    if(abs(fnext ) < tiny(1.0_dp)) then
       c = xnext
       return
    endif

    do i = 1, maxiter_

       if( abs(fnext - fstart) < tiny(1.0_dp)) then
          call stop_error("secant: division by zero")
       endif
       c = (xstart*fnext - xnext*fstart)/(fnext - fstart)
       if (abs(c - xnext) < tol) return ! root found
       !update variables
       xstart = xnext
       fstart = fnext
       xnext = c
       fnext = f(c)
    enddo
    !max iterations number reached
    call stop_error("secant: method did not converge")
  end function secant

  real(dp) function newton(f, df, a,b,tol, maxiter) result (c)
    !Solves f(x) = 0 on using the a and b as starting values
    procedure(func) :: f,df
    real(dp),intent(in) :: a,b,tol
    integer,optional :: maxiter

    integer :: maxiter_
    real(dp) :: xstart, xnext
    real(dp) :: fx,dfx,fstart,fnext
    integer ::  i

    if(present(maxiter)) then
       maxiter_ = maxiter
    else
       maxiter_ = 100
    endif


    xstart = a
    xnext = b
    fstart = f(xstart)
    if(abs(fstart ) < tiny(1.0_dp)) then
       c = xstart
       return
    endif
    fnext = f(xnext)
    if(abs(fnext ) < tiny(1.0_dp)) then
       c = xnext
       return
    endif

    c = 0.5_dp*(b + a)
    xstart = c

    do i = 1, maxiter_

       fx = f(xstart)
       dfx = df(xstart)
       if( abs(dfx) < tiny(1.0_dp)) then
          call stop_error("newton: division by zero")
       endif

       c = c - fx/dfx
       if (abs(c - xstart) < tol) return ! root found
       !update variables
       xstart = c
    enddo
    !max iterations number reached
    call stop_error("newton: method did not converge")
  end function newton

   function brent(f,a,b,tol,maxiter) result(xzero)
    use iso_fortran_env, only: output_unit
    procedure(func) :: f
    real(dp),intent(in) :: a,b,tol
    integer,optional :: maxiter
    real(dp) :: xzero,fzero
    integer :: maxiter_,iflag

    if(present(maxiter)) then
       maxiter_ = maxiter
    else
       maxiter_ = 100
    endif

    call zeroin(f,a,b,tol,xzero,fzero,iflag,maxiter_)

    if(iflag == -2) write(output_unit,'(a)') "brent: exceeded number of iterations"
    
    if(iflag == -1) then
       call stop_error("brent: division by zero")
    end if


  end function brent
  
  subroutine zeroin(f,ax,bx,tol,xzero,fzero,iflag,maxiter,fax,fbx)
    use iso_fortran_env, only: error_unit
    real(dp),parameter :: one=1.0_dp,zero=0.0_dp,two=2.0_dp,three=3.0_dp
    procedure(func) :: f
    real(dp),intent(in)              :: ax      !! left endpoint of initial interval
    real(dp),intent(in)              :: bx      !! right endpoint of initial interval
    real(dp),intent(in)              :: tol     !! desired length of the interval of uncertainty of the final result (>=0)
    integer,intent(in)               :: maxiter !! max. number of iterations
    real(dp),intent(out)             :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(dp),intent(out)             :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)              :: iflag   !! status flag (`-1`=error, `0`=root found)
    real(dp),intent(in),optional     :: fax     !! if `f(ax)` is already known, it can be input here
    real(dp),intent(in),optional     :: fbx     !! if `f(ax)` is already known, it can be input here

    real(dp),parameter :: eps   = epsilon(one)  !! original code had d1mach(4)
    real(dp) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s
    integer :: iter

    tol1 = eps+one

    a=ax
    b=bx

    if (present(fax)) then
        fa = fax
    else
        fa=f(a)
    end if
    if (present(fbx)) then
        fb = fbx
    else
        fb=f(b)
    end if

    !check trivial cases first:
    if (fa==zero) then

        iflag = 0
        xzero = a
        fzero = fa

    elseif (fb==zero) then

        iflag = 0
        xzero = b
        fzero = fb

    elseif (fa*(fb/abs(fb))<zero) then  ! check that f(ax) and f(bx) have different signs

        c=a
        fc=fa
        d=b-a
        e=d

        do iter=1,maxiter

            if (abs(fc)<abs(fb)) then
                a=b
                b=c
                c=a
                fa=fb
                fb=fc
                fc=fa
            end if

            tol1=two*eps*abs(b)+0.5_dp*tol
            xm = 0.5_dp*(c-b)
            if ((abs(xm)<=tol1).or.(fb==zero)) exit

            ! see if a bisection is forced
            if ((abs(e)>=tol1).and.(abs(fa)>abs(fb))) then
                s=fb/fa
                if (a/=c) then
                    ! inverse quadratic interpolation
                    q=fa/fc
                    r=fb/fc
                    p=s*(two*xm*q*(q-r)-(b-a)*(r-one))
                    q=(q-one)*(r-one)*(s-one)
                else
                    ! linear interpolation
                    p=two*xm*s
                    q=one-s
                end if
                if (p<=zero) then
                    p=-p
                else
                    q=-q
                end if
                s=e
                e=d
                if (((two*p)>=(three*xm*q-abs(tol1*q))) .or. &
                    (p>=abs(0.5_dp*s*q))) then
                    d=xm
                    e=d
                else
                    d=p/q
                end if
            else
                d=xm
                e=d
            end if

            a=b
            fa=fb
            if (abs(d)<=tol1) then
                if (xm<=zero) then
                    b=b-tol1
                else
                    b=b+tol1
                end if
            else
                b=b+d
            end if
            fb=f(b)
            if ((fb*(fc/abs(fc)))>zero) then
                c=a
                fc=fa
                d=b-a
                e=d
            end if

        end do

        iflag = 0
        xzero = b
        fzero = fb

        if(iter >= maxiter) iflag = -2
        
    else

        iflag = -1
        write(error_unit,'(A)')&
            'Error in zeroin: f(ax) and f(bx) do not have different signs.'

    end if

    end subroutine zeroin

end module Mroot

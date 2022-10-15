module scitools_filon
!! Provides Filon quadrature for integrals of the type 
!! \(I(\omega) = \int^b_a dt\, f(t)\cos(\omega t)\) or \(I(\omega) = \int^b_a dt\, f(t)\sin(\omega t)\).
!======================================================================================
   use scitools_def,only: dp
   implicit none
!--------------------------------------------------------------------------------------  
   private
   public :: Filon_quad
!--------------------------------------------------------------------------------------  
contains
!------------------------------------------------------------------------------------
   subroutine Filon_quad(f,a,b,omega,sc_flag,nmax,epsabs,res,abserr,neval)
   !! Computes the integral
   !! $$I(\omega) = \int^b_a dt\, f(t)\mathrm{osc}(\omega t) , $$
   !! where \(\mathrm{osc}(x) = \cos(x)\) or \(\mathrm{osc}(x) = \sin(x)\),
   !! using the Filon quadrature method.
      interface 
      !! interface for the the function scalar real \(f(t)\)
         function f(t) result(ft)
            use scitools_def,only: dp
            real(dp),intent(in) :: t
            real(dp) :: ft
         end function f
      end interface
      real(dp),intent(in) :: a,b  !! integration bounds
      real(dp),intent(in) :: omega !! the frequency \(\omega\)
      integer,intent(in) :: sc_flag !! `sc_flag = 1`: cos-integral, `sc_flag = 2`: sin-integral
      integer,intent(in) :: nmax !! maximum number of sample points
      real(dp),intent(in) :: epsabs !! absolute error threshold
      real(dp),intent(out) :: res !! the integral value \(I\)
      real(dp),intent(out) :: abserr !! absolute error
      integer,intent(out) :: neval !! number of function evaluations needed
      integer :: n,i
      real(dp) :: x,res_old,err
      real(dp),allocatable :: ftab(:)

      allocate(ftab(nmax))

      res_old = HUGE(1.0_dp)
      neval = 0

      NTABLOOP:do n=51,nmax,50
         do i=1,n
            x = a + (b-a)*(i-1)/dble(n-1)
            ftab(i) = f(x)
            neval = neval + 1
         end do
         if(sc_flag == 1) then
            call filon_cos(n,ftab,a,b,omega,res)
         else 
            call filon_sin(n,ftab,a,b,omega,res)
         end if
         err = abs(res_old - res)
         if(err < epsabs) then
            exit NTABLOOP
         else
            res_old = res
         end if
      end do NTABLOOP

      abserr = abs(res_old - res)

      deallocate(ftab)

   end subroutine Filon_quad
!--------------------------------------------------------------------------------------  
  subroutine filon_cos ( ntab, ftab, a, b, t, result )
   !!  FILON_COS uses the Filon method on integrals with a cosine factor.
   !!
   !!  Discussion
   !!  ----------
   !!   The integral to be approximated has the form:
   !!
   !!    Integral ( A <= X <= B ) F(X) * COS(T*X) dX
   !!
   !!   where T is user specified.
   !!
   !!  The function is interpolated over each subinterval by
   !!  a parabolic arc.
   !!
   !!  Modified:
   !!
   !!    10 February 2006
   !!
   !!  Author:
   !!
   !!    John Burkardt
   !!
   !!  Reference
   !!  ---------
   !!
   !!  Milton Abramowitz, Irene Stegun,
   !!  Handbook of Mathematical Functions,
   !!  National Bureau of Standards, 1964,
   !!  ISBN: 0-486-61272-4,
   !!  LC: QA47.A34.
   !!
   !!  Stephen Chase, Lloyd Fosdick,
   !!  An Algorithm for Filon Quadrature,
   !!  Communications of the Association for Computing Machinery,
   !!  Volume 12, Number 8, August 1969, pages 453-457.
   !!
   !!  Stephen Chase, Lloyd Fosdick,
   !!  Algorithm 353:
   !!  Filon Quadrature,
   !!  Communications of the Association for Computing Machinery,
   !!  Volume 12, Number 8, August 1969, pages 457-458.
   !!
   !!  Philip Davis, Philip Rabinowitz,
   !!  Methods of Numerical Integration,
   !!  Second Edition,
   !!  Dover, 2007,
   !!  ISBN: 0486453391,
   !!  LC: QA299.3.D28.
   !!
   !!  Parameters
   !!  ----------
   !!    Input, integer ( kind = 4 ) NTAB, the number of data points.
   !!    NTAB must be odd, and greater than 1.
   !!
   !!    Input, real ( kind = 8 ) FTAB(NTAB), contains the value of the function
   !!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
   !!
   !!    Input, real ( kind = 8 ) A, B, the limits of integration.
   !!
   !!    Input, real ( kind = 8 ) T, the multiplier of the X argument of the cosine.
   !!
   !!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
     implicit none

     integer ( kind = 4 ) ntab

     real ( kind = 8 ) a
     real ( kind = 8 ) alpha
     real ( kind = 8 ) b
     real ( kind = 8 ) beta
     real ( kind = 8 ) c2n
     real ( kind = 8 ) c2nm1
     real ( kind = 8 ) cost
     real ( kind = 8 ) ftab(ntab)
     real ( kind = 8 ) gamma
     real ( kind = 8 ) h
     real ( kind = 8 ) result
     real ( kind = 8 ) sint
     real ( kind = 8 ) t
     real ( kind = 8 ) theta
     real ( kind = 8 ) xtab(ntab)

     if ( a == b ) then
       result = 0.0D+00
       return
     end if
    
     if ( ntab <= 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FILON_COS - Fatal error!'
       write ( *, '(a)' ) '  NTAB < 2'
       write ( *, '(a,i8)' ) '  NTAB = ', ntab
       stop 1
     end if
    
     if ( mod ( ntab, 2 ) /= 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FILON_COS - Fatal error!'
       write ( *, '(a)' ) '  NTAB must be odd.'
       write ( *, '(a,i8)' ) '  NTAB = ', ntab
       stop 1
     end if
   !
   !  Set up a vector of the NTAB X values.
   ! 
     call r8vec_even ( ntab, a, b, xtab )

     h = ( b - a ) / real ( ntab - 1, kind = 8 )

     theta = t * h
     sint = sin ( theta )
     cost = cos ( theta )

     if ( 6.0D+00 * abs ( theta ) <= 1.0D+00 ) then

       alpha = 2.0D+00 * theta**3 /   45.0D+00 &
             - 2.0D+00 * theta**5 /  315.0D+00 &
             + 2.0D+00 * theta**7 / 4725.0D+00
     
       beta =  2.0D+00            /     3.0D+00 &
             + 2.0D+00 * theta**2 /    15.0D+00 &
             - 4.0D+00 * theta**4 /   105.0D+00 &
             + 2.0D+00 * theta**6 /   567.0D+00 &
             - 4.0D+00 * theta**8 / 22275.0D+00

       gamma = 4.0D+00            /      3.0D+00 &
             - 2.0D+00 * theta**2 /     15.0D+00 &
             +           theta**4 /    210.0D+00 &
             -           theta**6 /  11340.0D+00

     else

       alpha = ( theta**2 + theta * sint * cost &
         - 2.0D+00 * sint**2 ) / theta**3

       beta = ( 2.0D+00 * theta + 2.0D+00 * theta * cost**2 &
         - 4.0D+00 * sint * cost ) / theta**3

       gamma = 4.0D+00 * ( sint - theta * cost ) / theta**3
     
     end if

     c2n = sum ( ftab(1:ntab:2) * cos ( t * xtab(1:ntab:2) ) ) &
       - 0.5D+00 * ( ftab(ntab) * cos ( t * xtab(ntab) ) &
                   + ftab(1) * cos ( t * xtab(1) ) )

     c2nm1 = sum ( ftab(2:ntab-1:2) * cos ( t * xtab(2:ntab-1:2) ) )
    
     result = h * ( &
         alpha * ( ftab(ntab) * sin ( t * xtab(ntab) ) & 
                 - ftab(1)    * sin ( t * xtab(1) ) ) &
       + beta * c2n &
       + gamma * c2nm1 )

   end subroutine filon_cos



   subroutine filon_sin ( ntab, ftab, a, b, t, result )
   !!  FILON_SIN uses the Filon method on integrals with a sine factor.
   !!
   !!  Discussion
   !!  ----------
   !!   The integral to be approximated has the form:
   !!
   !!    Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
   !!
   !!   where T is user specified.
   !!
   !!  The function is interpolated over each subinterval by
   !!  a parabolic arc.
   !!
   !!  Modified:
   !!
   !!    10 February 2006
   !!
   !!  Author:
   !!
   !!    John Burkardt
   !!
   !!  Reference
   !!  ---------
   !!
   !!  Milton Abramowitz, Irene Stegun,
   !!  Handbook of Mathematical Functions,
   !!  National Bureau of Standards, 1964,
   !!  ISBN: 0-486-61272-4,
   !!  LC: QA47.A34.
   !!
   !!  Stephen Chase, Lloyd Fosdick,
   !!  An Algorithm for Filon Quadrature,
   !!  Communications of the Association for Computing Machinery,
   !!  Volume 12, Number 8, August 1969, pages 453-457.
   !!
   !!  Stephen Chase, Lloyd Fosdick,
   !!  Algorithm 353:
   !!  Filon Quadrature,
   !!  Communications of the Association for Computing Machinery,
   !!  Volume 12, Number 8, August 1969, pages 457-458.
   !!
   !!  Philip Davis, Philip Rabinowitz,
   !!  Methods of Numerical Integration,
   !!  Second Edition,
   !!  Dover, 2007,
   !!  ISBN: 0486453391,
   !!  LC: QA299.3.D28.
   !!
   !!  Parameters
   !!  ----------
   !!    Input, integer ( kind = 4 ) NTAB, the number of data points.
   !!    NTAB must be odd, and greater than 1.
   !!
   !!    Input, real ( kind = 8 ) FTAB(NTAB), contains the value of the function
   !!    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(NTAB-1).
   !!
   !!    Input, real ( kind = 8 ) A, B, the limits of integration.
   !!
   !!    Input, real ( kind = 8 ) T, the multiplier of the X argument of the cosine.
   !!
   !!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.
     implicit none

     integer ( kind = 4 ) ntab

     real ( kind = 8 ) a
     real ( kind = 8 ) alpha
     real ( kind = 8 ) b
     real ( kind = 8 ) beta
     real ( kind = 8 ) cost
     real ( kind = 8 ) ftab(ntab)
     real ( kind = 8 ) gamma
     real ( kind = 8 ) h
     real ( kind = 8 ) result
     real ( kind = 8 ) s2n
     real ( kind = 8 ) s2nm1
     real ( kind = 8 ) sint
     real ( kind = 8 ) t
     real ( kind = 8 ) theta
     real ( kind = 8 ) xtab(ntab)

     if ( a == b ) then
       result = 0.0D+00
       return
     end if
    
     if ( ntab <= 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
       write ( *, '(a)' ) '  NTAB < 2'
       write ( *, '(a,i8)' ) '  NTAB = ',ntab
       stop 1
     end if
    
     if ( mod ( ntab, 2 ) /= 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'FILON_SIN - Fatal error!'
       write ( *, '(a)' ) '  NTAB must be odd.'
       write ( *, '(a,i8)' ) '  NTAB = ',ntab
       stop 1
     end if
   !
   !  Set up a vector of the NTAB X values.
   ! 
     call r8vec_even ( ntab, a, b, xtab )

     h = ( b - a ) / real ( ntab - 1, kind = 8 )
     theta = t * h

     sint = sin ( theta )
     cost = cos ( theta )

     if ( 6.0D+00 * abs ( theta ) <= 1.0D+00 ) then

       alpha = 2.0D+00 * theta**3 /   45.0D+00 &
             - 2.0D+00 * theta**5 /  315.0D+00 &
             + 2.0D+00 * theta**7 / 4725.0D+00
     
       beta =  2.0D+00            /     3.0D+00 &
             + 2.0D+00 * theta**2 /    15.0D+00 &
             - 4.0D+00 * theta**4 /   105.0D+00 &
             + 2.0D+00 * theta**6 /   567.0D+00 &
             - 4.0D+00 * theta**8 / 22275.0D+00

       gamma = 4.0D+00            /      3.0D+00 &
             - 2.0D+00 * theta**2 /     15.0D+00 &
             +           theta**4 /    210.0D+00 &
             -           theta**6 /  11340.0D+00

     else
    
       alpha = ( theta**2 + theta * sint * cost &
         - 2.0D+00 * sint**2 ) / theta**3

       beta = ( 2.0D+00 * theta + 2.0D+00 * theta * cost**2 &
         - 4.0D+00 * sint * cost ) / theta**3

       gamma = 4.0D+00 * ( sint - theta * cost ) / theta**3
    
     end if
     
     s2n = sum ( ftab(1:ntab:2) * sin ( t * xtab(1:ntab:2) ) ) &
       - 0.5D+00 * ( ftab(ntab) * sin ( t * xtab(ntab) ) &
                   + ftab(1) * sin ( t * xtab(1) ) )

     s2nm1 = sum ( ftab(2:ntab-1:2) * sin ( t * xtab(2:ntab-1:2) ) )

     result = h * ( &
         alpha * ( ftab(1) * cos ( t * xtab(1) ) &
                 - ftab(ntab) * cos ( t * xtab(ntab) ) ) &
       + beta * s2n &
       + gamma * s2nm1 )

   end subroutine filon_sin


   subroutine r8vec_even ( n, alo, ahi, a )
   !! Returns N values, evenly spaced between ALO and AHI.
     implicit none

     integer ( kind = 4 ) n !! The number of values. Normally, `A(1) = ALO` and `A(N) = AHI`.
                            !! However, if N = 1, then `A(1) = 0.5*(ALO+AHI)`.
     real ( kind = 8 ) a(n) !! N evenly spaced values. 
     real ( kind = 8 ) ahi !! Input, high value
     real ( kind = 8 ) alo !! Input, low value
     integer ( kind = 4 ) i

     if ( n == 1 ) then

       a(1) = 0.5D+00 * ( alo + ahi )

     else

       do i = 1, n
         a(i) = ( real ( n - i,     kind = 8 ) * alo   &
                + real (     i - 1, kind = 8 ) * ahi ) &
                / real ( n     - 1, kind = 8 )
       end do

     end if

   end subroutine r8vec_even
!-------------------------------------------------------------------------------------- 


!======================================================================================
end module scitools_filon
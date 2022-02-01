module MMatrixExp
!----------------------------------------------------------------------|
  use Mdef,only:&
       dp,one,zero
  implicit none
!----------------------------------------------------------------------|
  private
  public &
       InitExpM,&
       CleanExpM,&
       ExpM,&
       ExpM_ThrSv
!----------------------------------------------------------------------|
  integer, parameter :: ideg = 6
  integer,allocatable :: iwsp(:)
  complex(dp),allocatable :: wsp1(:)
!----------------------------------------------------------------------|
contains
!----------------------------------------------------------------------|
  pure elemental complex(dp) function my_cosh(z)
    complex(dp),intent(in) :: z

#ifdef COSHFIX
    my_cosh = 0.5_dp * (exp(z) + exp(-z))
#else
    my_cosh = cosh(z)
#endif

  end function my_cosh
!----------------------------------------------------------------------|
  pure elemental complex(dp) function my_sinh(z)
    complex(dp),intent(in) :: z

#ifdef COSHFIX
    my_sinh = 0.5_dp * (exp(z) - exp(-z))
#else
    my_sinh = sinh(z)
#endif

  end function my_sinh
!----------------------------------------------------------------------|
!----------------------------------------------------------------------|
  subroutine InitExpM(n)
    integer,intent(in)::n

    allocate(iwsp(N),wsp1(4*N**2+ideg+1))
    
  end subroutine InitExpM
!----------------------------------------------------------------------|
  subroutine CleanExpM
    deallocate(iwsp,wsp1)
  end subroutine CleanExpM
!----------------------------------------------------------------------|
  function expm(t, H) result(expH)

    real(dp), intent(in) :: t
    complex(dp), dimension(:,:), intent(in) :: H
    complex(dp), dimension(size(H,1),size(H,2)) :: expH

    ! Expokit variables
    ! external :: ZGPADM
    integer :: iexp, ns, iflag, n

    n = size(H,1)
    
    call ZGPADM(ideg, n, t, H, n, wsp1, size(wsp1,1), iwsp, iexp, ns, iflag)
    expH = reshape(wsp1(iexp:iexp+n*n-1), shape(expH))

  end function expm
!----------------------------------------------------------------------|
  function expm_thrsv(t, H) result(expH)
    real(dp),parameter :: eps=1.0e-6_dp
    real(dp), intent(in) :: t
    complex(dp), dimension(:,:), intent(in) :: H
    complex(dp), dimension(size(H,1),size(H,2)) :: expH
  
    ! Expokit variables
    ! external :: ZGPADM
    integer :: iexp, ns, iflag, n
    integer, parameter :: ideg = 6
    integer,allocatable :: iwsp(:)
    complex(dp),allocatable :: wsp1(:)
    complex(dp) :: deth,ct,st,exad

    n = size(H,1)

    if(n==1) then
       expH=exp(t*h)
       return
    end if
    if(n==2) then
       deth = t * sqrt((h(1,1)-h(2,2))**2 + 4.0D0 * h(1,2)*h(2,1))
       exad = exp(t * 0.5D0*(h(1,1)+h(2,2)))
       if(abs(deth) < eps) then
          ct = 1D0 + (0.5D0*deth)**2/2D0 + (0.5D0*deth)**4/24D0 
          st = 0.5D0 + deth**2/48D0 + deth**4/3840D0
       else
          ct = my_cosh(0.5D0*deth) ; st = my_sinh(0.5D0*deth)/deth
       end if
       expH(1,1) = ct + (h(1,1)-h(2,2))*st*t
       expH(1,2) = 2D0*h(1,2) * st * t
       expH(2,1) = 2D0*h(2,1) * st * t
       expH(2,2) = ct - (h(1,1)-h(2,2))*st * t
       expH = exad * expH
       return
    end if
    
    allocate(iwsp(N),wsp1(4*N**2+ideg+1))
    
    call ZGPADM(ideg, n, t, H, n, wsp1, size(wsp1,1), iwsp, iexp, ns, iflag)
    expH = reshape(wsp1(iexp:iexp+n*n-1), shape(expH))

    deallocate(iwsp,wsp1)
    
  end function expm_thrsv
!----------------------------------------------------------------------|
  subroutine ZGPADM(ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)

    real(dp),intent(in) :: t
    integer,intent(in) :: ideg, m, ldh, lwsp
    integer,intent(inout) :: ipiv(m)
    integer,intent(out) :: iexph,ns,iflag
    complex(dp),intent(in)::H(:,:)
    complex(dp),intent(inout)::wsp(:)

!-----Purpose----------------------------------------------------------|
!
!     Computes exp(t*H), the matrix exponential of a general complex 
!     matrix in full, using the irreducible rational Pade approximation
!     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
!     combined with scaling-and-squaring.
!
!-----Arguments--------------------------------------------------------|
!
!     ideg      : (input) the degre of the diagonal Pade to be used.
!                 a value of 6 is generally satisfactory.
!
!     m         : (input) order of H.
!
!     H(ldh,m)  : (input) argument matrix.
!
!     t         : (input) time-scale (can be < 0).
!                  
!     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
!
!     ipiv(m)   : (workspace)
!
!>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
!                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
!                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!                 NOTE: if the routine was called with wsp(iptr), 
!                       then exp(tH) will start at wsp(iptr+iexph-1).
!
!     ns        : (output) number of scaling-squaring used.
!
!     iflag     : (output) exit flag.
!                       0 - no problem
!                      <0 - problem
!
!----------------------------------------------------------------------|
!     Roger B. Sidje (rbs@maths.uq.edu.au)
!     EXPOKIT: Software Package for Computing Matrix Exponentials.
!     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
!----------------------------------------------------------------------|
      integer ::i,j,k,icoef,mm,ih2,iodd,iused,ifree,iq,ip,iput,iget
      real(dp):: hnorm
      complex(dp):: cp, cq, scale, scale2


      intrinsic ABS, CMPLX, DBLE, INT, LOG, MAX
!---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
     
!
!---  initialise pointers ...
!
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
!
!---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
!     and set scale = t/2^ns ...
!
      do i = 1,m
         wsp(i) = ZERO
      enddo
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         enddo
      enddo
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,DBLE(wsp(i)) )
      enddo
      hnorm = ABS( t*hnorm )
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale =  CMPLX( t/DBLE(2**ns),0.0d0 )
      scale2 = scale*scale
!
!---  compute Pade coefficients ...
!
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = ONE
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      enddo
!
!---  H2 = scale2*H*H ...
!
      call ZGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,ZERO,wsp(ih2),m )
!
!---  initialise p (numerator) and q (denominator) ...
!
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = ZERO
            wsp(iq + (j-1)*m + i-1) = ZERO
         enddo
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      enddo
!
!---  Apply Horner rule ...
!
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call ZGEMM( 'n','n',m,m,m, ONE,wsp(iused),m,wsp(ih2),m, ZERO,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      enddo
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
!
!---  Obtain (+/-)(I + 2*(p\q)) ...
!
      if ( iodd.ne.0 ) then
         call ZGEMM( 'n','n',m,m,m, scale,wsp(iq),m,H,ldh, ZERO,wsp(ifree),m )
         iq = ifree
      else
         call ZGEMM( 'n','n',m,m,m, scale,wsp(ip),m,H,ldh, ZERO,wsp(ifree),m )
         ip = ifree
      endif
      call ZAXPY( mm, -ONE,wsp(ip),1, wsp(iq),1 )
      call ZGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
      call ZDSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + ONE
      enddo
      iput = ip
      if ( ns.eq.0 .and. iodd.ne.0 ) then
         call ZDSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
!
!--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
!
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call ZGEMM( 'n','n',m,m,m, ONE,wsp(iget),m, wsp(iget),m,ZERO,wsp(iput),m )
         iodd = 1-iodd
      enddo
 200  continue
      iexph = iput
      END SUBROUTINE ZGPADM
!----------------------------------------------------------------------|



      
end module MMatrixExp

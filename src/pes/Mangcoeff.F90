module Mangcoeff
!======================================================================================
   use scitools_def,only: dp,iu,zero,one
   implicit none
   include "../units_inc.f90"
!--------------------------------------------------------------------------------------
   private 
   public :: &
      TransformationMatrix_Y2X,&
      Transform_Y2X,&
      Log3J,&
      ClebGord,&
      ThreeYlm,&
      ThreeXlm
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   pure elemental integer function minus_one_pwr(m)
      integer,intent(in) :: m

      minus_one_pwr = -2 * mod(m,2) + 1

   end function minus_one_pwr
!--------------------------------------------------------------------------------------
   subroutine TransformationMatrix_Y2X(l,R)
      real(dp),parameter :: sq2=sqrt(2.0_dp)
      integer,intent(in) :: l
      complex(dp),intent(inout) :: R(:,:)
      integer :: i,m,i_neg,nl

      nl = 2 * l + 1
      R = zero
      do i=1,nl
         R(i,i) = one
      end do

      do m=-l,l
         i = m + l + 1
         i_neg = nl - i + 1
         if(m < 0) then
            R(i_neg,i) = iu/sq2 * minus_one_pwr(m+1)
            R(i,i) = iu/sq2
         elseif(m > 0) then
            R(i_neg,i) = 1.0_dp/sq2
            R(i,i) = minus_one_pwr(m)/sq2
         end if
      end do

   end subroutine TransformationMatrix_Y2X
!--------------------------------------------------------------------------------------
   pure elemental function Transform_Y2X(M,Cy_M,Cy_mm) result(Cx)
      integer,intent(in)     :: M
      complex(dp),intent(in) :: Cy_m,Cy_mm
      complex(dp)            :: Cx
      integer :: s

      s = 1
      if(mod(M,2) .ne. 0) s = -1

      if(M == 0) then
         Cx = Cy_m
      elseif(M < 0) then
         Cx = iu/sqrt(2.0_dp) * (Cy_mm - s * Cy_m)
      else
         Cx = 1.0_dp/sqrt(2.0_dp) * (Cy_mm + s * Cy_m)
      end if

   end function Transform_Y2X
!--------------------------------------------------------------------------------------
   pure elemental real(dp) function LogFactorial(n)
      integer,intent(in) :: n
      integer  :: i
      real(dp) :: z

      if(n <= 0) then
         LogFactorial = 0.0_dp
         return
      end if

      z = 0.0_dp
      do i=1,n
         z = z + log(dble(i))
      end do

      LogFactorial = z

   end function LogFactorial
!--------------------------------------------------------------------------------------
   pure logical function triagineq(j1,j2,j3)
      integer,intent(in) :: j1,j2,j3

      triagineq = (abs(j1-j2)<=j3).and.(j3<=j1+j2)
      triagineq = triagineq.and.(abs(j2-j3)<=j1).and.(j1<=j2+j3)
      triagineq = triagineq.and.(abs(j1-j3)<=j2).and.(j2<=j1+j3)

   end function triagineq
!--------------------------------------------------------------------------------------
   pure real(dp) function Log3J(j1,m1,j2,m2,j3,m3)
      integer,intent(in) :: j1,m1,j2,m2,j3,m3
      integer :: s,sup(3),smax,p(5),k
      real(dp) :: z(10),s3j,delta,beta
      logical :: neg

      If((m1+m2+m3>0).or.(m1+m2+m3)<0) Then
         Log3J = 0D0
         Return
      End If

      If(.not.triagineq(j1,j2,j3)) Then
         Log3J = 0D0
         Return
      End If

      sup(1)=j1+j2-j3
      sup(2)=j1-m1
      sup(3)=j2+m2
      smax=Maxval(sup)
      s3j=0D0

      z(1)=LogFactorial(j1+j2-j3)
      z(2)=LogFactorial(j1-j2+j3)
      z(3)=LogFactorial(j2+j3-j1)
      z(4)=LogFactorial(j1+j2+j3+1)
      delta=0.5D0*(z(1)+z(2)+z(3)-z(4))

      z(1)=LogFactorial(j1+m1)
      z(2)=LogFactorial(j1-m1)
      z(3)=LogFactorial(j2+m2)
      z(4)=LogFactorial(j2-m2)
      z(5)=LogFactorial(j3+m3)
      z(6)=LogFactorial(j3-m3)
      beta=0.5D0*(Sum(z(1:6)))

      Do s=0,smax

         p(1)=j1+j2-j3-s
         p(2)=j1-m1-s
         p(3)=j2+m2-s
         p(4)=j3-j2+m1+s
         p(5)=j3-j1-m2+s

         neg=(p(1)<0).or.(p(2)<0).or.(p(3)<0).or.(p(4)<0).or.(p(5)<0)

         If(.not.neg) Then
            Do k=1,5
               z(k)=-LogFactorial(p(k))
            End Do
            z(6)=-LogFactorial(s)    
            s3j=s3j+(-1)**(j1-j2-m3+s)*exp(sum(z(1:6))+beta+delta)          
         End If
      End Do

      Log3J=s3j

   end function Log3J
!--------------------------------------------------------------------------------------
   pure real(dp) function ClebGord(l1,m1,l2,m2,l3,m3)
      integer,intent(in) :: l1,m1,l2,m2,l3,m3

      ClebGord=sqrt(2D0*l3+1D0)*(-1)**(l1-l2+m3)
      ClebGord=ClebGord*Log3J(l1,m1,l2,m2,l3,-m3)

  end function ClebGord
!--------------------------------------------------------------------------------------
   pure real(dp) function ThreeYlm(l1,m1,l2,m2,l3,m3)
      integer,intent(in)  :: l1,m1,l2,m2,l3,m3

      ThreeYlm=sqrt((2D0*l1+1D0)*(2D0*l2+1D0)*(2D0*l3+1D0)/qpi)
      ThreeYlm=ThreeYlm*Log3J(l1,0,l2,0,l3,0)*Log3J(l1,m1,l2,m2,l3,m3)

   end function ThreeYlm
!--------------------------------------------------------------------------------------
   pure real(dp) function ThreeXlm(l1,m1,l2,m2,l3,m3)
      real(dp),parameter :: c = 1.0_dp/sqrt(2.0_dp)
      integer,intent(in)  :: l1,m1,l2,m2,l3,m3
      complex(dp) :: a(3),b(3),xlm

      !================================
      ! real spherical harmonics Xlm:
      ! Xlm = (-1)**m*sqrt(2)*Im(Ylm) : m<0
      ! Xlm = Ylm : m=0
      ! Xlm = sqrt(2)*Re(Ylm)
      !================================

      a=1D0
      b=0D0
      if(m1<0) then
         a(1)=(-1)**m1*c/iu
         b(1)=-c/iu
      end if
      if(m1>0) then
         a(1)=c
         b(1)=(-1)**m1*c
      end if
      if(m2<0) then
         a(2)=(-1)**m2*c/iu
         b(2)=-c/iu
      end if
      if(m2>0) then
         a(2)=c
         b(2)=(-1)**m2*c
      end if
      if(m3<0) then
         a(3)=(-1)**m3*c/iu
         b(3)=-c/iu
      end if
      if(m3>0) then
         a(3)=c
         b(3)=(-1)**m3*c
      end if

      xlm=a(1)*a(2)*a(3)*ThreeYlm(l1,m1,l2,m2,l3,m3)
      xlm=xlm+b(1)*a(2)*a(3)*ThreeYlm(l1,-m1,l2,m2,l3,m3)
      xlm=xlm+a(1)*b(2)*a(3)*ThreeYlm(l1,m1,l2,-m2,l3,m3)
      xlm=xlm+a(1)*a(2)*b(3)*ThreeYlm(l1,m1,l2,m2,l3,-m3)
      xlm=xlm+b(1)*b(2)*a(3)*ThreeYlm(l1,-m1,l2,-m2,l3,m3)
      xlm=xlm+b(1)*a(2)*b(3)*ThreeYlm(l1,-m1,l2,m2,l3,-m3)
      xlm=xlm+a(1)*b(2)*b(3)*ThreeYlm(l1,m1,l2,-m2,l3,-m3)
      xlm=xlm+b(1)*b(2)*b(3)*ThreeYlm(l1,-m1,l2,-m2,l3,-m3)

      ThreeXlm=dble(xlm)

   end function ThreeXlm
!--------------------------------------------------------------------------------------


!======================================================================================
end module Mangcoeff
module Mevol
!======================================================================================
!! This module defines routines constructing and applying the time-evolution
!! operator.  
   use Mdef,only: dp,iu,one,zero
   use Mlinalg,only: blas_zgemm
   use Mmatrixexp,only: ExpM, ExpM_Thrsv
   implicit none
!--------------------------------------------------------------------------------------
   private
   public :: UnitaryStepFW, UnitaryStepBW, UnitaryStepFBW, GenU_CF2, GenU_CF4
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   pure subroutine UnitaryStepFW(n,Ut,Xt,Xtdt,large_size)
      !! peforms the operation \(X(t+dt) = U(dt)X(t)\)
      integer,intent(in)                        :: n !! matrix dimension
      complex(dp),dimension(:,:),intent(in)     :: Ut !! time evolution operator U(dt)
      complex(dp),dimension(:,:),intent(in)     :: Xt !! matrix X(t)
      complex(dp),dimension(:,:),intent(inout)  :: Xtdt !! matrix X(t+dt)
      logical,intent(in),optional :: large_size !! if .true., ZGEMM will be called for matrix-matrix multiplication
      logical :: use_zgemm

      use_zgemm = .false.
      if(present(large_size)) use_zgemm = large_size

      if(use_zgemm) then
        call blas_zgemm('n','n',n,n,n,one,Ut,n,Xt,n,zero,Xtdt,n)
      else
        Xtdt = matmul(Ut,Xt)
      end if

   end subroutine UnitaryStepFW
!--------------------------------------------------------------------------------------
   pure subroutine UnitaryStepBW(n,Ut,Xt,Xtdt,large_size)
      !! peforms the operation \(X(t+dt) = X(t)U^\dagger(dt)\)
      integer,intent(in)                        :: n !! matrix dimension
      complex(dp),dimension(:,:),intent(in)     :: Ut !! time evolution operator U(dt)
      complex(dp),dimension(:,:),intent(in)     :: Xt !! matrix X(t)
      complex(dp),dimension(:,:),intent(inout)  :: Xtdt !! matrix X(t+dt)
      logical,intent(in),optional :: large_size !! if .true., ZGEMM will be called for matrix-matrix multiplication
      complex(dp),dimension(n,n)::Utcc
      logical :: use_zgemm

      use_zgemm = .false.
      if(present(large_size)) use_zgemm = large_size

      if(use_zgemm) then
        call blas_zgemm('n','c',n,n,n,one,Xt,n,Ut,n,zero,Xtdt,n)
      else
        Utcc = conjg(transpose(Ut))
        Xtdt = matmul(Xt,Utcc)
      end if

   end subroutine UnitaryStepBW
!--------------------------------------------------------------------------------------
   pure subroutine UnitaryStepFBW(n,Ut,Xt,Xtdt,large_size)
      !! peforms the operation X(t+dt) = \(U(dt)X(t)U^\dagger(dt)\)
      integer,intent(in)                        :: n !! matrix dimension
      complex(dp),dimension(:,:),intent(in)     :: Ut !! time evolution operator U(dt)
      complex(dp),dimension(:,:),intent(in)     :: Xt !! matrix X(t)
      complex(dp),dimension(:,:),intent(inout)  :: Xtdt !! matrix X(t+dt)
      logical,intent(in),optional :: large_size !! if .true., ZGEMM will be called for matrix-matrix multiplication
      complex(dp),dimension(n,n)::Utcc,Xtmp
      logical :: use_zgemm

      use_zgemm = .false.
      if(present(large_size)) use_zgemm = large_size

      if(use_zgemm) then
        call blas_zgemm('n','n',n,n,n,one,Ut,n,Xt,n,zero,Xtmp,n)
        call blas_zgemm('n','c',n,n,n,one,Xtmp,n,Ut,n,zero,Xtdt,n)
      else
        Utcc=conjg(transpose(Ut))
        Xtmp = matmul(Ut, Xt)
        Xtdt = matmul(Xtmp, Utcc)
      end if

   end subroutine UnitaryStepFBW
!--------------------------------------------------------------------------------------    
   subroutine GenU_CF2(dt,Hdt2,Udt)
      !! Second-order commutator-free time evolution operator
      !! (mid-point approximation)
      !!  U(t+dt,t) = exp[-i * dt * H(t+dt/2)]
      real(dp),intent(in)       :: dt !! time step
      complex(dp),intent(in)    :: Hdt2(:,:) !! Hamiltonian at t+dt/2
      complex(dp),intent(inout) :: Udt(:,:) !! time evolution operator

      if(size(Hdt2,1)==1) then
         Udt(1,1)=exp(-iu*dt*Hdt2(1,1))
         return
      end if

      Udt(:,:)=-iu*dt*Hdt2
      Udt(:,:)=ExpM_Thrsv(1.0_dp,Udt(:,:))

   end subroutine GenU_CF2
!--------------------------------------------------------------------------------------
   subroutine GenU_CF4(dt,H1,H2,Udt,large_size)
      !! Fourth-order commutator-free time evolution operator from J. Comp. Phys. 230, 5930 (2011):
      !! $$ U(t+dt,t) = \exp[-i dt (a_1 H_1 + a_2*H_2] \exp[-i dt (a_2 H_1 + a_1 H_2)] , $$ where
      !! \(H_1 = H(t + c_1 dt)\), \(H_2 = H(t + c_2 dt)\) and
      !! $$ c_1 = \frac{1}{2} - \frac{\sqrt{3}}{6}, c_2 = \frac{1}{2} + \frac{\sqrt{3}}{6} $$
      real(dp),parameter::sq3=sqrt(3.0_dp)
      real(dp),parameter::a1=(3.0_dp-2.0_dp*sq3)/12.0_dp
      real(dp),parameter::a2=(3.0_dp+2.0_dp*sq3)/12.0_dp
      real(dp),intent(in)         :: dt !! time step
      complex(dp),intent(in)      :: H1(:,:) !! Hamiltonian at t + c1 * dt
      complex(dp),intent(in)      :: H2(:,:) !! Hamiltonian at t + c2 * dt
      complex(dp),intent(inout)   :: Udt(:,:) !! time evolution operator
      logical,intent(in),optional :: large_size !! if .true., ZGEMM will be called for matrix-matrix multiplication
      logical :: use_zgemm
      integer :: n
      complex(dp),dimension(size(H1,1),size(H1,1))::U1,U2

      use_zgemm = .false.
      if(present(large_size)) use_zgemm = large_size

      n = size(H1,1)

      if(n==1) then
         U1=exp(-iu*dt*(a1*H1(1,1)+a2*H2(1,1)))
         U2=exp(-iu*dt*(a2*H1(1,1)+a1*H2(1,1)))
         Udt(1,1)=U1(1,1)*U2(1,1)
         return
      end if

      U1=-iu*dt*(a1*H1+a2*H2)
      U2=-iu*dt*(a2*H1+a1*H2)    
      U1=ExpM_Thrsv(1.0_dp,U1)
      U2=ExpM_Thrsv(1.0_dp,U2)

      if(use_zgemm) then
         call blas_zgemm('N','N',n,n,n,one,U1,n,U2,n,zero,Udt,n)
      else
         Udt = matmul(U1,U2)
      end if

   end subroutine GenU_CF4
!--------------------------------------------------------------------------------------
 
  
!======================================================================================
end module Mevol

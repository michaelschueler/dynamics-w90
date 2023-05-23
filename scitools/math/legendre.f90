module scitools_legendre
!! Contains some tools for constructing Legendre polynomials
!======================================================================================
   use scitools_debug
   use scitools_def,only: dp
   implicit none
!--------------------------------------------------------------------------------------
   private 
   public :: legendre, legendre_roots
!--------------------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------------------
   elemental real(dp) function legendre(n, x)
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: x
      real(dp) :: p0, p1, p2
      integer :: i

      if (n == 0) then
         legendre = 1.0_dp
      else if (n == 1) then
         legendre = x
      else
         p0 = 1.0_dp
         p1 = x
         do i = 2, n
            p2 = ((2.0_dp*i-1.0_dp)*x*p1 - (i-1.0_dp)*p0) / i
            p0 = p1
            p1 = p2
         end do
      legendre = p2
      end if

   end function legendre
!--------------------------------------------------------------------------------------
   function legendre_roots(n) result(roots)
      implicit none
      integer, intent(in) :: n
      real(dp), dimension(n) :: roots
      real(dp), dimension(n,n) :: mat, tau
      integer :: lwork, info
      integer :: i, j
      real(dp),allocatable :: work(:)

      ! Initialize the matrix
      do i = 1, n
         do j = 1, n
            if (i == j) then
               mat(i,j) = (2.0_dp*i-1.0_dp) / i
            else if (j == i+1) then
               mat(i,j) = 0.5_dp*sqrt((i*(i+1.0_dp)))
            else if (j == i-1) then
               mat(i,j) = 0.5_dp*sqrt((i*(i+1.0_dp)))
            else
               mat(i,j) = 0.0_dp
            end if
         end do
      end do

      ! Compute the QR factorization using LAPACK
      ! call dgeqrf(n, n, mat, n, tau, roots, n, info)
      allocate(work(1))
      lwork = -1
      call DGEQRF(n, n, mat, n, tau, work, lwork, info)
      lwork = int(work(1))
      deallocate(work); allocate(work(lwork))
      call DGEQRF(n, n, mat, n, tau, work, lwork, info)
      deallocate(work)

      ! Extract the roots from the diagonal of the upper triangular matrix
      do i = 1, n
         roots(i) = mat(i,i)
      end do

   end function legendre_roots
!--------------------------------------------------------------------------------------

!======================================================================================
end module scitools_legendre
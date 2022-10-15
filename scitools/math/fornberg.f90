module scitools_fornberg
!! Contains routines for generation interpolation and differentation weights 
!! by the Fornberg algorithm.
!======================================================================================
  use scitools_def,only: dp
  implicit none
!--------------------------------------------------------------------------------------
  private
  public :: fornberg_weights
!--------------------------------------------------------------------------------------  
contains
!--------------------------------------------------------------------------------------  
  pure subroutine fornberg_weights(z, x, nd, m, c)
  !! Generates interpolation and finite-differences differentation weights \(c^{(m)}_j\)
  !! at point \(z\), provided the functional values have tabulated at \(y_j=y(x_j), \ j=0,n\).
  !! The \(m\)-th derivative (\(m=0\) corresponds to interpolation) at \(z\) is expressed
  !! as 
  !! $$ y^{(m)}(z) = \sum^n_{j=0}c^{(m)}_j y_j .$$
  !!  Reference:
  !!      Generation of Finite Difference Formulas on Arbitrarily
  !!          Spaced Grids, Bengt Fornberg,
  !!         Mathematics of compuation, 51, 184, 1988, 699-706
    real(dp), intent(in)    :: z !! location where approximations are to be accurate
    integer, intent(in)     :: nd !! dimension of x- and c-arrays in calling 
                                  !! program x(0:nd) and c(0:nd, 0:m), respectively
    integer, intent(in)     :: m !! highest derivative for which weights are sought
    real(dp), intent(in)    :: x(0:) !! grid point locations, found in x(0:n)
    real(dp), intent(out) :: c(0:, 0:) !! weights at grid locations x(0:n) for 
                                       !! derivatives of order 0:m, found in c(0:nd, 0:m)

    real(dp) :: c1, c2, c3, c4, c5
    integer :: i, j, k, mn

    c1 = 1
    c4 = x(0) - z
    c = 0
    c(0, 0) = 1
    do i=1, nd
      mn = min(i, m)
      c2 = 1
      c5 = c4
      c4 = x(i) - z
      do j=0, i-1
        c3 = x(i) - x(j)
        c2 = c2*c3
        if (j == i-1) then
          do k = mn, 1, -1
            c(i, k) = c1*(k*c(i-1, k-1) - c5*c(i-1, k))/c2
          end do
          c(i, 0) = -c1*c5*c(i-1, 0)/c2
        endif
        do k=mn, 1, -1
          c(j, k) = (c4*c(j, k) - k*c(j, k-1))/c3
        end do
        c(j, 0) = c4*c(j, 0)/c3
      end do
      c1 = c2
    end do
  end subroutine fornberg_weights
!-------------------------------------------------------------------------------------
  
!======================================================================================
end module scitools_fornberg

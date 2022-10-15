module Mdebug
!======================================================================================
!! Provides tools for debugging. The debug checks are disabled if compiled in release
!! mode, which is trigged by the precompiler flag -DNDEBUG
!======================================================================================
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use scitools_def,only: dp

#ifdef NDEBUG
   logical,parameter :: debug_mode=.false.
#else
   logical,parameter :: debug_mode=.true.
#endif
!--------------------------------------------------------------------------------------
   private
   public :: debug_mode, assert, assert_shape
!-------------------------------------------------------------------------------------- 
   interface assert_shape
      !! Checks if given matrix has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
      module procedure dassert_shape_1
      module procedure zassert_shape_1
      module procedure dassert_shape_2
      module procedure zassert_shape_2
      module procedure dassert_shape_3
      module procedure zassert_shape_3
      module procedure dassert_shape_4
      module procedure zassert_shape_4
      module procedure dassert_shape_5
      module procedure zassert_shape_5
      module procedure dassert_shape_6
      module procedure zassert_shape_6
   end interface assert_shape
!--------------------------------------------------------------------------------------       
contains
!-------------------------------------------------------------------------------------- 
   subroutine assert(condition,message)
      !! Checks if given condition is fulfilled; otherwise an error is returned
      !! Only active if compiled in debug mode.
      logical,intent(in) :: condition !! logical expression for the condition to be checked
      character(len=*),intent(in),optional :: message !! return message

      if(debug_mode) then

         if(.not.condition) then
            if(present(message)) then
               write(error_unit,'("[ERROR]",1x,A)') 'Assertion '//message//' failed!'
            else
               write(error_unit,'("[ERROR]",1x,A)') 'Assertion failed!'
            end if

            stop

         end if

      end if

   end subroutine assert
!--------------------------------------------------------------------------------------
   subroutine dassert_shape_1(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
      real(dp), intent(in) :: A(:) !! the vector to be checked
      integer, intent(in) :: shap(:) !! the shape the vector should have
      character(len=*) :: routine !! name of routine where this assertion is called
      character(len=*) :: matname !! name of the vector A

      if(debug_mode)then

         if(any(shape(A) /= shap)) then
            write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
            write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
            stop
         end if
      end if

   end subroutine dassert_shape_1
!--------------------------------------------------------------------------------------
  subroutine zassert_shape_1(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    complex(dp), intent(in) :: A(:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode) then
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    endif
  end subroutine zassert_shape_1
!--------------------------------------------------------------------------------------
  subroutine dassert_shape_2(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    real(dp), intent(in) :: A(:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode)then
    
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    end if
       
  end subroutine dassert_shape_2
!--------------------------------------------------------------------------------------
  subroutine zassert_shape_2(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    complex(dp), intent(in) :: A(:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode) then
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    endif
  end subroutine zassert_shape_2
!--------------------------------------------------------------------------------------
  subroutine dassert_shape_3(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    real(dp), intent(in) :: A(:,:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode)then
    
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    end if
       
  end subroutine dassert_shape_3
!--------------------------------------------------------------------------------------
  subroutine zassert_shape_3(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    complex(dp), intent(in) :: A(:,:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode) then
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    endif
  end subroutine zassert_shape_3
!--------------------------------------------------------------------------------------
  subroutine dassert_shape_4(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    real(dp), intent(in) :: A(:,:,:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode)then
    
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    end if
       
  end subroutine dassert_shape_4
!--------------------------------------------------------------------------------------
  subroutine zassert_shape_4(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    complex(dp), intent(in) :: A(:,:,:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode) then
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    endif
  end subroutine zassert_shape_4
!--------------------------------------------------------------------------------------
  subroutine dassert_shape_5(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    real(dp), intent(in) :: A(:,:,:,:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode)then
    
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    end if
       
  end subroutine dassert_shape_5
!--------------------------------------------------------------------------------------
  subroutine zassert_shape_5(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    complex(dp), intent(in) :: A(:,:,:,:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode) then
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    endif
  end subroutine zassert_shape_5
!--------------------------------------------------------------------------------------
  subroutine dassert_shape_6(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    real(dp), intent(in) :: A(:,:,:,:,:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode)then
    
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    end if
       
  end subroutine dassert_shape_6
!--------------------------------------------------------------------------------------
  subroutine zassert_shape_6(A, shap, routine, matname)
      !! Checks if given real vector has a given shape; otherwise an error is returned
      !! Only active if compiled in debug mode.
    complex(dp), intent(in) :: A(:,:,:,:,:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(debug_mode) then
       if(any(shape(A) /= shap)) then
          write(error_unit,'("[ERROR]",1x,A)')  "In routine " // trim(routine) // " matrix " &
            // trim(matname) // " has illegal shape ", shape(A)
          write(error_unit,'("[ERROR]",1x,A)')  "Shape should be ", shap
          stop
       end if
    endif
  end subroutine zassert_shape_6
!--------------------------------------------------------------------------------------
  
!======================================================================================  
end module Mdebug

module Mdebug
!======================================================================================
!! tools for debugging  
!======================================================================================
  use,intrinsic::iso_fortran_env,only:&
       output_unit,error_unit

#ifdef NDEBUG
  logical,parameter :: debug_mode=.false.
#else
  logical,parameter :: debug_mode=.true.
#endif
!--------------------------------------------------------------------------------------
  private
  public &
       debug_mode,&
       assert,&
       assert_shape
!-------------------------------------------------------------------------------------- 
  ! assert shape of matrices:
  interface assert_shape
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
    logical,intent(in) :: condition
    character(len=*),intent(in),optional :: message

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
    ! make sure a given real matrix has a given shape
    real(8), intent(in) :: A(:)
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
       
  end subroutine dassert_shape_1
!--------------------------------------------------------------------------------------
  subroutine zassert_shape_1(A, shap, routine, matname)
    ! make sure a given real matrix has a given shape
    complex(8), intent(in) :: A(:)
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
    ! make sure a given real matrix has a given shape
    real(8), intent(in) :: A(:,:)
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
    ! make sure a given real matrix has a given shape
    complex(8), intent(in) :: A(:,:)
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
    ! make sure a given real matrix has a given shape
    real(8), intent(in) :: A(:,:,:)
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
    ! make sure a given real matrix has a given shape
    complex(8), intent(in) :: A(:,:,:)
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
    ! make sure a given real matrix has a given shape
    real(8), intent(in) :: A(:,:,:,:)
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
    ! make sure a given real matrix has a given shape
    complex(8), intent(in) :: A(:,:,:,:)
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
    ! make sure a given real matrix has a given shape
    real(8), intent(in) :: A(:,:,:,:,:)
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
    ! make sure a given real matrix has a given shape
    complex(8), intent(in) :: A(:,:,:,:,:)
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
    ! make sure a given real matrix has a given shape
    real(8), intent(in) :: A(:,:,:,:,:,:)
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
    ! make sure a given real matrix has a given shape
    complex(8), intent(in) :: A(:,:,:,:,:,:)
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

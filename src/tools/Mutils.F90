module Mutils
!======================================================================================
   !! Various general utilities.
   !! Based on a code by John E. Pask, LLNL.
   use,intrinsic::iso_fortran_env,only: output_unit,error_unit
   use Mdef, only: dp
   implicit none
   include "../formats.h"
!-------------------------------------------------------------------------------------- 
   private
   public &
      upcase, &
      lowcase, &
      whitechar, &
      blank, &
      numstrings, &
      getstring, &
      stop_error, &
      loadtxt, &
      savetxt, &
      newunit, &
      str, &
      load_griddata, &
      save_griddata, &
      checkoption, &
      linspace, &
      print_header, &
      print_tag, &
      print_title, &
      get_file_ext,&
      remove_file_ext, &
      check_file_ext
!-------------------------------------------------------------------------------------- 
   interface str
      module procedure str_int, str_real, str_real_n
   end interface str

   interface load_griddata
      module procedure load_griddata_real, load_griddata_complex
   end interface load_griddata

   interface save_griddata
      module procedure save_griddata_real_1d, save_griddata_complex_1d, &
         save_griddata_real_2d, save_griddata_complex_2d
   end interface save_griddata
!-------------------------------------------------------------------------------------- 
contains
!-------------------------------------------------------------------------------------- 
  function linspace(xmin, xmax, nx) result(x)
   !! analogous to numpys linspace function: returns an equally-spaced array 
   !! from xmin to xmax with nx points
    real(dp),intent(in) :: xmin !! The starting value of the sequence.
    real(dp),intent(in) :: xmax !! The end value of the sequence.
    integer,intent(in)  :: nx !! number of points
    real(dp),allocatable :: x(:)
    integer :: i

    if(.not.allocated(x)) allocate(x(nx))
    do i=1,nx
      x(i) = xmin + (xmax - xmin) * (i-1)/dble(nx-1)
    end do

  end function linspace
!-------------------------------------------------------------------------------------- 
  subroutine print_title(unit,title)
   !! prints the title of the program with date and time
    character(len=*),parameter :: fmt_time='(a,"/",a,"/",a," at ",a,":",a,":",a)'
    integer,intent(in) :: unit
    character(len=*),intent(in) :: title
    character(len=8) :: date
    character(len=10) :: time
    
    write(unit,'("+",68("-"),"+")') 
    write(unit,'(A1)',advance='no') "|"
    call print_tag(unit,title,width=68,advance='no')
    write(unit,'(A1)') "|"
    write(unit,'("+",68("-"),"+")') 
    call date_and_time(date=date,time=time)
    write(unit,'(a)',advance='no') '  Calculation started on '
    write(unit,fmt_time) date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:6)
    write(unit,*)
  end subroutine
!-------------------------------------------------------------------------------------- 
  subroutine print_tag(unit,tag,width,advance)
   !! prints a header of a section of the program
    integer,intent(in) :: unit
    character(len=*),intent(in) :: tag
    integer,intent(in),optional :: width
    character(len=*),intent(in),optional :: advance
    character(len=3) :: advance_
    integer :: w,l,l2
    character(len=100) :: fmt

    w = 70
    if(present(width)) w = width

    advance_ = 'yes'
    if(present(advance)) advance_ = advance

    l = w - len_trim(tag)
    l2 = floor(l/2.0d0)

    fmt = "(" // trim(adjustl(str_int(l2-1))) // "x" // ",A," //trim(adjustl(str_int(l-l2+1))) //"x" // ")"

    write(unit,fmt,advance=trim(advance_)) trim(tag)

  end subroutine print_tag
!-------------------------------------------------------------------------------------- 
  subroutine print_header(unit,tag,symbol,width)
   !! prints a header of a section of the program
    integer,intent(in) :: unit
    character(len=*),intent(in) :: tag
    character(len=1),intent(in) :: symbol
    integer,intent(in),optional :: width
    integer :: w,l,l2
    character(len=100) :: fmt

    w = 70
    if(present(width)) w = width

    l = w - len_trim(tag)
    l2 = floor(l/2.0d0)

    fmt = '(' // trim(adjustl(str_int(l2-1))) // '(' // '"' // symbol // '")' // ',1x,A,1x,' // &
      trim(adjustl(str_int(l-l2-1))) // '("' // symbol // '"))'

    write(unit,fmt) trim(tag)

  end subroutine print_header
!-------------------------------------------------------------------------------------- 
  function upcase(s) result(t)
    !! Returns string 's' in uppercase
    character(*), intent(in) :: s
    character(len(s)) :: t
    integer :: i, diff
    t = s; diff = ichar('A')-ichar('a')
    do i = 1, len(t)
       if (ichar(t(i:i)) >= ichar('a') .and. ichar(t(i:i)) <= ichar('z')) then
          ! if lowercase, make uppercase
          t(i:i) = char(ichar(t(i:i)) + diff)
       end if
    end do
  end function upcase
!-------------------------------------------------------------------------------------- 
  function lowcase(s) result(t)
    !! Returns string 's' in lowercase
    character(*), intent(in) :: s
    character(len(s)) :: t
    integer :: i, diff
    t = s; diff = ichar('A')-ichar('a')
    do i = 1, len(t)
       if (ichar(t(i:i)) >= ichar('A') .and. ichar(t(i:i)) <= ichar('Z')) then
          ! if uppercase, make lowercase
          t(i:i) = char(ichar(t(i:i)) - diff)
       end if
    end do
  end function lowcase
!-------------------------------------------------------------------------------------- 
  logical function checkoption(s,sopt)
  !! checks if string 's' is equivalent to string 'sopt' after converting to lower case
    character(*), intent(in) :: s,sopt

    checkoption = trim(lowcase(s)) == trim(lowcase(sopt))
    
  end function checkoption
!--------------------------------------------------------------------------------------  
  logical function whitechar(char) 
    !! returns .true. if char is space (32) or tab (9), .false. otherwise
    character, intent(in) :: char
    if (iachar(char) == 32 .or. iachar(char) == 9) then
       whitechar = .true.
    else
       whitechar = .false.
    end if
  end function whitechar
!-------------------------------------------------------------------------------------- 
  logical function blank(string)
    !! Returns true if string contains only white characters
    character(*), intent(in) :: string
    integer :: i
    do i = 1, len(string)
       if (.not. whitechar(string(i:i))) exit
    end do
    blank = (i>len(string))
  end function blank
!-------------------------------------------------------------------------------------- 
  integer function numstrings(s) result(n)
    !! Returns number of substrings contained in input string 's' delimited
    !! by white space.
    character(*), intent(in) :: s    !! input string
    character(len(s)+2) :: t         !! temporary string to facilitate analysis
    integer :: i
    t = " " // s // " "
    n = 0
    do i = 1, len(t)-1
       if (whitechar(t(i:i)) .and. .not. whitechar(t(i+1:i+1))) n = n + 1
    end do
  end function numstrings
  !--------------------------------------------------------------------------------------------------!
  subroutine getstring(s,is,ss)
    !! Returns first substring ss in string s, delimited by white space, starting at
    !! index is in s. If ss is found, is is set to (index of last character of ss in
    !! s) + 1; else is is set to 0. If is is out of range on input, routine
    !! terminates with is = -1.
    character(*), intent(in) :: s   !! input string
    integer, intent(inout) :: is    !! on input: starting index for search for ss in
    !! s on output: (index of last character of ss in
    !! s) + 1
    character(*), intent(out) :: ss !! first substring in s, starting from index is
    character(len(s)+1) :: t        !! temporary string to facilitate search
    integer i, i1, i2
    logical prevwhite, curwhite
    if (is <= 0 .or. is > len(s)) then
       ss = ""; is = -1; return
    end if
    t = s // " "
    if (is == 1) then
       prevwhite = .true.
    else
       prevwhite = whitechar(t(is-1:is-1))
    end if
    i1 = 0; i2 = 0
    do i = is, len(t)
       curwhite = whitechar(t(i:i))
       if (prevwhite .and. .not. curwhite) i1 = i   ! beginning of substring
       if (i1>0 .and. curwhite) then                ! end of substring
          i2 = i-1; exit
       end if
       prevwhite=curwhite
    end do
    if (i2 > 0) then
       ss = t(i1:i2); is = i2+1
    else
       ss = ""; is = 0
    end if
  end subroutine getstring
!-------------------------------------------------------------------------------------- 
  integer function newunit(unit) result(n)
    !! Returns lowest i/o unit number not in use (to be used in older compilers).
    !! Starting at 10 to avoid lower numbers which are sometimes reserved.
    !! Note: largest valid unit number may be system-dependent.
    integer, intent(out), optional :: unit

    logical inuse
    integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter :: nmax=999  ! may be system-dependent
    do n = nmin, nmax
       inquire(unit=n, opened=inuse)
       if (.not. inuse) then
          if (present(unit)) unit=n
          return
       end if
    end do
    call stop_error("newunit ERROR: available unit not found.")
  end function newunit
!-------------------------------------------------------------------------------------- 
   subroutine stop_error(msg,root_flag)
   !! Aborts the program with nonzero exit code
   !! The statement "stop msg" will return 0 exit code when compiled using
   !! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
   !! 1 and a print statement to print the message.
   !! If compiled with MPI, the program will terminate all ranks
#ifdef MPI
      use mpi
#endif
      character(len=*) :: msg !! Message to print on error output
      logical,intent(in),optional :: root_flag !! flag for MPI execution: message is only printed on root
      logical :: root_
      integer :: ierr

      root_ = .true.
      if(present(root_flag)) root_ = root_flag

      if(root_) write(error_unit,fmt900) msg

#ifdef MPI
      call MPI_Finalize(ierr)
#endif

      stop 1
   end subroutine stop_error
!-------------------------------------------------------------------------------------- 
  subroutine loadtxt(filename, d, trans)
    !! Loads a 2D array from a text file.
    character(len=*), intent(in) :: filename !! Filename to load the array from
    real(dp), allocatable, intent(out) :: d(:, :) !! The array 'd' will be automatically allocated with the correct dimensions
    logical,intent(in),optional :: trans !! option to transpose the output array
    logical ::  trans_
    character :: c
    integer :: s, ncol, nrow, ios, i
    logical :: lastwhite
    real(dp) :: r

    trans_ = .false.
    if(present(trans)) trans_ = trans

    open(newunit=s, file=filename, status="old")

    ! determine number of columns
    ncol = 0
    lastwhite = .true.
    do
       read(s, '(a)', advance='no', iostat=ios) c
       if (ios /= 0) exit
       if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
       lastwhite = whitechar(c)
    end do

    rewind(s)

    ! determine number or rows
    nrow = 0
    do
       read(s, *, iostat=ios) r
       if (ios /= 0) exit
       nrow = nrow + 1
    end do

    rewind(s)

    if(trans_) then
      allocate(d(ncol, nrow))
      do i = 1, nrow
         read(s, *) d(:, i)
      end do
    else
      allocate(d(nrow, ncol))
      do i = 1, nrow
         read(s, *) d(i, :)
      end do
    end if
    close(s)
  end subroutine loadtxt
!-------------------------------------------------------------------------------------- 
  subroutine savetxt(filename, d, fmt, transp)
    !! Saves a 2D array into a textfile.
    character(len=*), intent(in) :: filename  !! File to save the array to
    real(dp), intent(in) :: d(:, :)           !! The 2D array to save
    character(len=*), intent(in), optional  :: fmt !! string format
    logical, intent(in), optional :: transp !! option to tranpose array d before writing
    integer :: s, i
    logical :: transp_


    transp_ = .false.
    if(present(transp)) transp_ = transp

    open(newunit=s, file=filename, status="replace")

    if(transp_) then
       if(present(fmt)) then
          do i = 1, size(d, 2)
             write(s, fmt) d(:, i)
          end do
       else
          do i = 1, size(d, 2)
             write(s, *) d(:, i)
          end do
       end if
    else
       if(present(fmt)) then
          do i = 1, size(d, 1)
             write(s, fmt) d(i, :)
          end do
       else
          do i = 1, size(d, 1)
             write(s, *) d(i, :)
          end do
       end if
    end if
    close(s)
  end subroutine savetxt
!-------------------------------------------------------------------------------------- 
   subroutine save_griddata_real_1d(filename, x, y)
    !! Saves real grid data x_i, y_i to file
    character(len=*), intent(in) :: filename !! File to save the data to
    real(dp), intent(in) :: x(:) !! x-values
    real(dp), intent(in) :: y(:) !! y-values
    integer :: s, i
    
    open(newunit=s, file=filename, status="replace")
    do i = 1, size(y, 1)
       write(s, *) x(i), y(i)
    end do
    close(s)
    
  end subroutine save_griddata_real_1d
!-------------------------------------------------------------------------------------- 
  subroutine save_griddata_complex_1d(filename, x, y)
    !! Saves complex grid data x_i, y_i to file
    character(len=*), intent(in) :: filename !! File to save the data to
    real(dp), intent(in) :: x(:)
    complex(dp), intent(in) :: y(:)
    integer :: s, i
    
    open(newunit=s, file=filename, status="replace")
    do i = 1, size(y, 1)
       write(s, *) x(i), dble(y(i)), aimag(y(i))
    end do
    close(s)
    
  end subroutine save_griddata_complex_1d
!-------------------------------------------------------------------------------------- 
  subroutine save_griddata_real_2d(filename, x, y, transp)
   !! Saves real grid data x_i, y_i to file where y_i is vector-values
    character(len=*), intent(in) :: filename !! File to save the data to
    real(dp), intent(in) :: x(:) !! x-values 
    real(dp), intent(in) :: y(:,:) !! y-values
    logical, intent(in), optional :: transp !! option to tranpose array y before writing
    logical :: transp_
    integer :: s, i

    transp_ = .false.
    if(present(transp)) transp_ = transp
    
    open(newunit=s, file=filename, status="replace")
    if(transp_) then
       do i = 1, size(y, 2)
          write(s, *) x(i), y(:, i)
       end do
    else
       do i = 1, size(y, 1)
          write(s, *) x(i), y(i, :)
       end do
    end if
    close(s)
    
  end subroutine save_griddata_real_2d
!-------------------------------------------------------------------------------------- 
  subroutine save_griddata_complex_2d(filename, x, y, transp)
   !! Saves complex grid data x_i, y_i to file where y_i is vector-values
    character(len=*), intent(in) :: filename !! File to save the data to
    real(dp), intent(in) :: x(:) !! x-values 
    complex(dp), intent(in) :: y(:,:)  !! y-values
    logical, intent(in), optional :: transp !! option to tranpose array y before writing
    logical :: transp_
    integer :: s, i

    transp_ = .false.
    if(present(transp)) transp_ = transp
    
    open(newunit=s, file=filename, status="replace")
    if(transp_) then
       do i = 1, size(y, 2)
          write(s, *) x(i), dble(y(:,i)), aimag(y(:,i))
       end do
    else
       do i = 1, size(y, 1)
          write(s, *) x(i), dble(y(i, :)), aimag(y(i,:))
       end do
    end if
    close(s)
    
  end subroutine save_griddata_complex_2d
!-------------------------------------------------------------------------------------- 
  subroutine load_griddata_real(filename, x, y)
   !! loads real grid data written by save_griddata_real_2d
    character(len=*), intent(in) :: filename  !! Filename to load the array from
    real(dp), allocatable, intent(out) :: x(:) !! x-values; will be automatically allocated with the correct dimensions
    real(dp), allocatable, intent(out) :: y(:, :) !! y-values; will be automatically allocated with the correct dimensions
    
    character :: c
    integer :: s, ncol, nrow, ios, i
    logical :: lastwhite
    real(dp) :: r

    open(newunit=s, file=filename, status="old")

    ! determine number of columns
    ncol = 0
    lastwhite = .true.
    do
       read(s, '(a)', advance='no', iostat=ios) c
       if (ios /= 0) exit
       if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
       lastwhite = whitechar(c)
    end do

    rewind(s)

    ! determine number or rows
    nrow = 0
    do
       read(s, *, iostat=ios) r
       if (ios /= 0) exit
       nrow = nrow + 1
    end do

    rewind(s)

    if(.not.allocated(x)) allocate(x(nrow))
    if(.not.allocated(y)) allocate(y(nrow, ncol))

    do i = 1, nrow
       read(s, *) x(i), y(i, :)
    end do
    close(s)
    
  end subroutine load_griddata_real
!-------------------------------------------------------------------------------------- 
  subroutine load_griddata_complex(filename, x, y)
   !! loads complex grid data written by save_griddata_complex_2d
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: x(:)
    complex(dp), allocatable, intent(out) :: y(:,:)
    
    character :: c
    integer :: s, ncol, nrow, ios, i
    logical :: lastwhite
    real(dp) :: r
    real(dp), allocatable :: colvec(:)
    
    open(newunit=s, file=filename, status="old")

    ! determine number of columns
    ncol = 0
    lastwhite = .true.
    do
       read(s, '(a)', advance='no', iostat=ios) c
       if (ios /= 0) exit
       if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
       lastwhite = whitechar(c)
    end do

    rewind(s)

    ! determine number or rows
    nrow = 0
    do
       read(s, *, iostat=ios) r
       if (ios /= 0) exit
       nrow = nrow + 1
    end do

    rewind(s)

    if(.not.allocated(x)) allocate(x(nrow))
    if(.not.allocated(y)) allocate(y(nrow, (ncol-1)/2))
    allocate(colvec(ncol-1))
    do i = 1, nrow
       read(s, *) x(i), colvec(1:ncol-1)
       y(i,1:(ncol-1)/2)=cmplx(colvec(1:(ncol-1)/2),colvec((ncol-1)/2+1:ncol-1),dp)
    end do
    close(s)
    deallocate(colvec)
    
  end subroutine load_griddata_complex
!-------------------------------------------------------------------------------------- 
  pure integer function str_int_len(i) result(sz)
    !! Returns the length of the string representation of 'i'
    !! If 's' is too short (MAX_STR too small), Fortan will abort with:
    !! "Fortran runtime error: End of record"
    integer, intent(in) :: i
    integer, parameter :: MAX_STR = 100
    character(MAX_STR) :: s

    write(s, '(i0)') i
    sz = len_trim(s)
  end function str_int_len
!-------------------------------------------------------------------------------------- 
  pure function str_int(i) result(s)
    !! Converts integer "i" to string
    integer, intent(in) :: i
    character(len=str_int_len(i)) :: s
    write(s, '(i0)') i
  end function str_int
!-------------------------------------------------------------------------------------- 
  pure integer function str_real_len(r, fmt) result(sz)
    !! Returns the length of the string representation of 'i'
    !! If 's' is too short (MAX_STR too small), Fortan will abort with:
    !! "Fortran runtime error: End of record"
    real(dp), intent(in) :: r
    character(len=*), intent(in) :: fmt
    integer, parameter :: MAX_STR = 100
    character(MAX_STR) :: s

    write(s, fmt) r
    sz = len_trim(s)
  end function str_real_len
!-------------------------------------------------------------------------------------- 
  pure function str_real(r) result(s)
    !! Converts the real number "r" to string with 7 decimal digits.
    real(dp), intent(in) :: r
    character(len=*), parameter :: fmt="(f0.6)"
    character(len=str_real_len(r, fmt)) :: s
    write(s, fmt) r
  end function str_real
!--------------------------------------------------------------------------------------
  pure function str_real_n(r, n) result(s)
    !! Converts the real number "r" to string with 'n' decimal digits.
    real(dp), intent(in) :: r
    integer, intent(in) :: n
    character(len=str_real_len(r, "(f0." // str_int(n) // ")")) :: s
    write(s, "(f0." // str_int(n) // ")") r
  end function str_real_n
!--------------------------------------------------------------------------------------
  function get_file_ext(fname) result(ext)
   !! returns the file extension of a file name
    character(len=*),intent(in) :: fname !! the file name
    character(len=32) :: ext !! file extention
    integer :: ppos

    ppos = scan(trim(fname),".", BACK= .true.)
    if( ppos > 0 ) then
      ext = trim(fname(ppos+1:))
    else
      ext = ""
    end if

  end function get_file_ext
!--------------------------------------------------------------------------------------
  function remove_file_ext(fname) result(output_file)
    !! returns file name with the file extension removed
    character(len=*),intent(in) :: fname !! the file name
    character(len=256) :: output_file !! file name without file extension
    integer :: ppos

    ppos = scan(trim(fname),".", BACK= .true.)
    if( ppos > 0 ) then
      output_file = trim(fname(1:ppos-1))
    else
      output_file = trim(fname)
    end if

  end function remove_file_ext
!--------------------------------------------------------------------------------------
  logical function check_file_ext(fname,ext)
  !! checks if given file name name has specific file extension 
    character(len=*),intent(in) :: fname !! the file name
    character(len=*),intent(in) :: ext !! file extention

    check_file_ext = checkoption(get_file_ext(fname), trim(ext))

  end function check_file_ext
!--------------------------------------------------------------------------------------

!======================================================================================
end module Mutils

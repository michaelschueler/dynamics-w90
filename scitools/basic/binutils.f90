module scitools_binutils
!======================================================================================
!! Provides routines for writing multi-dimensional grid data to binary files. 
!! For distinguishing scalar, vector or matrix-valued quantities 
!! the integer flag DTYPE states the type of data is contained. The reading routines
!! read this flag first and then refer to specializations for the different data types:
!!
!!  * DTYPE =  0  ...  grid only
!!  * DTYPE =  1  ...  real scalar
!!  * DTYPE =  2  ...  complex scalar
!!  * DTYPE =  3  ...  real vector
!!  * DTYPE =  4  ...  complex vector
!!  * DTYPE =  5  ...  real matrix
!!  * DTYPE =  6  ...  complex matrix      
!======================================================================================  
  use scitools_def,only: dp
  use scitools_utils,only: newunit
  implicit none
!--------------------------------------------------------------------------------------   
  private
  public &
       WriteData1D_to_BinaryFile,&
       WriteData2D_to_BinaryFile,&
       ReadData1D_from_BinaryFile,&
       ReadData2D_from_BinaryFile
!-------------------------------------------------------------------------------------- 
  interface WriteData1D_to_BinaryFile
  !! Interface for writing data \(x_i, y_i\) to binary file, where \(x_i\) is a real grid and
  !! \(y_i\) can be a real or complex scalar / vector / square matrix
     module procedure &
          WriteData1D_to_BinaryFile_dtype0,&
          WriteData1D_to_BinaryFile_dtype1,&
          WriteData1D_to_BinaryFile_dtype2,&
          WriteData1D_to_BinaryFile_dtype3,&
          WriteData1D_to_BinaryFile_dtype4,&
          WriteData1D_to_BinaryFile_dtype5,&
          WriteData1D_to_BinaryFile_dtype6
  end interface WriteData1D_to_BinaryFile

  interface WriteData2D_to_BinaryFile
  !! Interface for writing data \(x_i, y_j, z_{ij}\) to binary file, where \(x_i, y_j\) is a real grid and
  !! \(z_{ij}\) can be a real or complex scalar / vector / square matrix
     module procedure &
          WriteData2D_to_BinaryFile_dtype0,&
          WriteData2D_to_BinaryFile_dtype1,&
          WriteData2D_to_BinaryFile_dtype2,&
          WriteData2D_to_BinaryFile_dtype3,&
          WriteData2D_to_BinaryFile_dtype4,&
          WriteData2D_to_BinaryFile_dtype5,&
          WriteData2D_to_BinaryFile_dtype6
  end interface WriteData2D_to_BinaryFile

  interface ReadData1D_from_BinaryFile
  !! Interface for reading data \(x_i, y_i\) from binary file, where \(x_i\) is a real grid and
  !! \(y_i\) can be a real or complex scalar / vector / square matrix 
     module procedure &
          ReadData1D_from_BinaryFile_dtype0,&
          ReadData1D_from_BinaryFile_dtype1,&
          ReadData1D_from_BinaryFile_dtype2,&
          ReadData1D_from_BinaryFile_dtype3,&
          ReadData1D_from_BinaryFile_dtype4,&
          ReadData1D_from_BinaryFile_dtype5,&
          ReadData1D_from_BinaryFile_dtype6
  end interface ReadData1D_from_BinaryFile

  interface ReadData2D_from_BinaryFile
  !! Interface for reading data \(x_i, y_j, z_{ij}\) from binary file, where \(x_i, y_j\) is a real grid and
  !! \(z_{ij}\) can be a real or complex scalar / vector / square matrix
     module procedure &
          ReadData2D_from_BinaryFile_dtype0,&
          ReadData2D_from_BinaryFile_dtype1,&
          ReadData2D_from_BinaryFile_dtype2,&
          ReadData2D_from_BinaryFile_dtype3,&
          ReadData2D_from_BinaryFile_dtype4,&
          ReadData2D_from_BinaryFile_dtype5,&
          ReadData2D_from_BinaryFile_dtype6
  end interface ReadData2D_from_BinaryFile
!-------------------------------------------------------------------------------------- 
contains
!--------------------------------------------------------------------------------------
!                        ++++   1D routines ++++ 
!--------------------------------------------------------------------------------------  
  subroutine WriteData1D_to_BinaryFile_dtype0(Flname,x)
    integer,parameter::dtype=0
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:) !! sample points \(x_i\)
    integer::Nx
    integer::Flunit

    Nx=size(x) 

    
    open(newunit(Flunit),form='UNFORMATTED',file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx
    write(Flunit) x
    
    close(FlUnit)

  end subroutine WriteData1D_to_BinaryFile_dtype0
!--------------------------------------------------------------------------------------
  subroutine WriteData1D_to_BinaryFile_dtype1(Flname,x,data)
    integer,parameter::dtype=1
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:) !! sample points \(x_i\)
    real(dp),intent(in)::data(:) !! sample values \(y_i\)
    integer::Nx
    integer::Flunit

    Nx=size(x)
    
    open(newunit(Flunit),form='UNFORMATTED',file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx
    write(Flunit) x
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData1D_to_BinaryFile_dtype1
!--------------------------------------------------------------------------------------
  subroutine WriteData1D_to_BinaryFile_dtype2(Flname,x,data)
    integer,parameter::dtype=2
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:) !! sample points \(x_i\)
    complex(dp),intent(in)::data(:) !! sample values \(y_i\)
    integer::Nx
    integer::Flunit

    Nx=size(x)
    
    open(newunit(Flunit),form='UNFORMATTED',file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx
    write(Flunit) x
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData1D_to_BinaryFile_dtype2
!--------------------------------------------------------------------------------------
  subroutine WriteData1D_to_BinaryFile_dtype3(Flname,x,data)
    integer,parameter::dtype=3
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:) !! sample points \(x_i\)
    real(dp),intent(in)::data(:,:) !! sample values \(y_i\)
    integer::Nx,size1
    integer::Flunit

    size1=size(data,1)
    Nx=size(x)
    
    open(newunit(Flunit),form='UNFORMATTED',file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,size1
    write(Flunit) x
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData1D_to_BinaryFile_dtype3
!--------------------------------------------------------------------------------------
  subroutine WriteData1D_to_BinaryFile_dtype4(Flname,x,data)
    integer,parameter::dtype=4
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:) !! sample points \(x_i\)
    complex(dp),intent(in)::data(:,:) !! sample values \(y_i\)
    integer::Nx,size1
    integer::Flunit

    size1=size(data,1)
    Nx=size(x)
    
    open(newunit(Flunit),form='UNFORMATTED',file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,size1
    write(Flunit) x
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData1D_to_BinaryFile_dtype4
!--------------------------------------------------------------------------------------  
  subroutine WriteData1D_to_BinaryFile_dtype5(Flname,x,data)
    integer,parameter::dtype=5
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:) !! sample points \(x_i\)
    real(dp),intent(in)::data(:,:,:) !! sample values \(y_i\)
    integer::Nx,size1,size2
    integer::Flunit

    size1=size(data,1) ; size2=size(data,2)
    Nx=size(x)
    
    open(newunit(Flunit),form='UNFORMATTED',file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,size1,size2
    write(Flunit) x
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData1D_to_BinaryFile_dtype5
!--------------------------------------------------------------------------------------
  subroutine WriteData1D_to_BinaryFile_dtype6(Flname,x,data)
    integer,parameter::dtype=6
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:) !! sample points \(x_i\)
    complex(dp),intent(in)::data(:,:,:) !! sample values \(y_i\)
    integer::Nx,size1,size2
    integer::Flunit

    size1=size(data,1) ; size2=size(data,2)
    Nx=size(x) 
    
    open(newunit(Flunit),form='UNFORMATTED',file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,size1,size2
    write(Flunit) x
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData1D_to_BinaryFile_dtype6
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
  subroutine ReadData1D_from_BinaryFile_dtype0(Flname,x,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:) !! sample points \(x_i\)
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,dtype_dummy
    integer::Flunit
    logical::isthere

    iflag=0
    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
   
    open(newunit(Flunit),form='UNFORMATTED',file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.0) then
       iflag=1
       close(FlUnit)
       return
    end if
    
    read(Flunit) Nx
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    
    close(FlUnit)

  end subroutine ReadData1D_from_BinaryFile_dtype0
!--------------------------------------------------------------------------------------
  subroutine ReadData1D_from_BinaryFile_dtype1(Flname,x,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:) !! sample points \(x_i\)
    real(dp),allocatable,intent(out)::data(:) !! sample values \(y_i\)
    integer,intent(out)::iflag  !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,dtype_dummy
    integer::Flunit
    logical::isthere

    iflag=0
    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.1) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(data)) allocate(data(Nx))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData1D_from_BinaryFile_dtype1
!--------------------------------------------------------------------------------------
  subroutine ReadData1D_from_BinaryFile_dtype2(Flname,x,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:) !! sample points \(x_i\)
    complex(dp),allocatable,intent(out)::data(:) !! sample values \(y_i\)
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,dtype_dummy
    integer::Flunit
    logical::isthere

    iflag=0
    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.2) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(data)) allocate(data(Nx))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData1D_from_BinaryFile_dtype2
!--------------------------------------------------------------------------------------
  subroutine ReadData1D_from_BinaryFile_dtype3(Flname,x,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:) !! sample points \(x_i\)
    real(dp),allocatable,intent(out)::data(:,:) !! sample values \(y_i\)
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,size1,dtype_dummy
    integer::Flunit
    logical::isthere

    iflag=0
    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.3) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx,size1
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(data)) allocate(data(size1,Nx))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData1D_from_BinaryFile_dtype3
!--------------------------------------------------------------------------------------
  subroutine ReadData1D_from_BinaryFile_dtype4(Flname,x,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:) !! sample points \(x_i\)
    complex(dp),allocatable,intent(out)::data(:,:)  !! sample values \(y_i\)
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,size1,dtype_dummy
    integer::Flunit
    logical::isthere

    iflag=0
    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.4) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx,size1
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(data)) allocate(data(size1,Nx))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData1D_from_BinaryFile_dtype4
!--------------------------------------------------------------------------------------  
  subroutine ReadData1D_from_BinaryFile_dtype5(Flname,x,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:) !! sample points \(x_i\)
    real(dp),allocatable,intent(out)::data(:,:,:) !! sample values \(y_i\) 
    integer,intent(out)::iflag  !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,size1,size2,dtype_dummy
    integer::Flunit
    logical::isthere

    iflag=0
    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.5) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx,size1,size2
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(data)) allocate(data(size1,size2,Nx))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData1D_from_BinaryFile_dtype5
!--------------------------------------------------------------------------------------
  subroutine ReadData1D_from_BinaryFile_dtype6(Flname,x,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:) !! sample points \(x_i\)
    complex(dp),allocatable,intent(out)::data(:,:,:) !! sample values \(y_i\) 
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,size1,size2,dtype_dummy
    integer::Flunit
    logical::isthere

    iflag=0
    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.6) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx,size1,size2
    if(.not.allocated(x)) allocate(x(Nx))   
    read(Flunit) x
    if(.not.allocated(data)) allocate(data(size1,size2,Nx))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData1D_from_BinaryFile_dtype6
!--------------------------------------------------------------------------------------
  

!--------------------------------------------------------------------------------------
!                        ++++   2D routines ++++ 
!--------------------------------------------------------------------------------------  
  subroutine WriteData2D_to_BinaryFile_dtype0(Flname,x,y)
    integer,parameter::dtype=0
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:),y(:) !! sample points \(x_i\), \(y_j\)
    integer::Nx,Ny
    integer::Flunit

    Nx=size(x) ; Ny = size(y)
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,Ny
    write(Flunit) x
    write(Flunit) y
    
    close(FlUnit)

  end subroutine WriteData2D_to_BinaryFile_dtype0
!--------------------------------------------------------------------------------------
  subroutine WriteData2D_to_BinaryFile_dtype1(Flname,x,y,data)
    integer,parameter::dtype=1
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:),y(:) !! sample points \(x_i\), \(y_j\)
    real(dp),intent(in)::data(:,:) !! sample values \(z_{ij}\)
    integer::Nx,Ny
    integer::Flunit

    Nx=size(x) ; Ny = size(y)
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,Ny
    write(Flunit) x
    write(Flunit) y
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData2D_to_BinaryFile_dtype1
!--------------------------------------------------------------------------------------
  subroutine WriteData2D_to_BinaryFile_dtype2(Flname,x,y,data)
    integer,parameter::dtype=2
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:),y(:) !! sample points \(x_i\), \(y_j\)
    complex(dp),intent(in)::data(:,:) !! sample values \(z_{ij}\)
    integer::Nx,Ny
    integer::Flunit

    Nx=size(x) ; Ny = size(y)
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,Ny
    write(Flunit) x
    write(Flunit) y
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData2D_to_BinaryFile_dtype2
!--------------------------------------------------------------------------------------
  subroutine WriteData2D_to_BinaryFile_dtype3(Flname,x,y,data)
    integer,parameter::dtype=3
    character(len=*),intent(in)::Flname !! file name for output
    real(dp),intent(in)::x(:),y(:) !! sample points \(x_i\), \(y_j\)
    real(dp),intent(in)::data(:,:,:) !! sample values \(z_{ij}\)
    integer::Nx,Ny,size1
    integer::Flunit

    size1=size(data,1)
    Nx=size(x) ; Ny = size(y)
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,Ny,size1
    write(Flunit) x
    write(Flunit) y
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData2D_to_BinaryFile_dtype3
!--------------------------------------------------------------------------------------
  subroutine WriteData2D_to_BinaryFile_dtype4(Flname,x,y,data)
    integer,parameter::dtype=4
    character(len=*),intent(in)::Flname  !! file name for output
    real(dp),intent(in)::x(:),y(:) !! sample points \(x_i\), \(y_j\)
    complex(dp),intent(in)::data(:,:,:) !! sample values \(z_{ij}\)
    integer::Nx,Ny,size1
    integer::Flunit

    size1=size(data,1)
    Nx=size(x) ; Ny = size(y)
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,Ny,size1
    write(Flunit) x
    write(Flunit) y
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData2D_to_BinaryFile_dtype4
!--------------------------------------------------------------------------------------  
  subroutine WriteData2D_to_BinaryFile_dtype5(Flname,x,y,data)
    integer,parameter::dtype=5
    character(len=*),intent(in)::Flname  !! file name for output
    real(dp),intent(in)::x(:),y(:) !! sample points \(x_i\), \(y_j\)
    real(dp),intent(in)::data(:,:,:,:) !! sample values \(z_{ij}\)
    integer::Nx,Ny,size1,size2
    integer::Flunit

    size1=size(data,1) ; size2=size(data,2)
    Nx=size(x) ; Ny = size(y)
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,Ny,size1,size2
    write(Flunit) x
    write(Flunit) y
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData2D_to_BinaryFile_dtype5
!--------------------------------------------------------------------------------------
  subroutine WriteData2D_to_BinaryFile_dtype6(Flname,x,y,data)
    integer,parameter::dtype=6
    character(len=*),intent(in)::Flname  !! file name for output
    real(dp),intent(in)::x(:),y(:) !! sample points \(x_i\), \(y_j\)
    complex(dp),intent(in)::data(:,:,:,:) !! sample values \(z_{ij}\)
    integer::Nx,Ny,size1,size2
    integer::Flunit

    size1=size(data,1) ; size2=size(data,2)
    Nx=size(x) ; Ny = size(y)
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    write(Flunit) dtype
    write(Flunit) Nx,Ny,size1,size2
    write(Flunit) x
    write(Flunit) y
    write(Flunit) data

    close(FlUnit)

  end subroutine WriteData2D_to_BinaryFile_dtype6
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine ReadData2D_from_BinaryFile_dtype0(Flname,x,y,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:),y(:) !! sample points \(x_i\) , \(y_j\)
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,Ny,dtype_dummy
    integer::Flunit
    logical::isthere
    
    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
   
    open(newunit(Flunit),form='UNFORMATTED',file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.0) then
       iflag=1
       close(FlUnit)
       return
    end if
    
    read(Flunit) Nx,Ny
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(y)) allocate(y(Ny))
    read(Flunit) y
    
    close(FlUnit)

  end subroutine ReadData2D_from_BinaryFile_dtype0
!--------------------------------------------------------------------------------------
  subroutine ReadData2D_from_BinaryFile_dtype1(Flname,x,y,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:),y(:) !! sample points \(x_i\) , \(y_j\)
    real(dp),allocatable,intent(out)::data(:,:) !! sample values \(z_{ij}\) 
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,Ny,dtype_dummy
    integer::Flunit
    logical::isthere

    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.1) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx,Ny
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(y)) allocate(y(Ny))
    read(Flunit) y
    if(.not.allocated(data)) allocate(data(Nx,Ny))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData2D_from_BinaryFile_dtype1
!--------------------------------------------------------------------------------------
  subroutine ReadData2D_from_BinaryFile_dtype2(Flname,x,y,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:),y(:) !! sample points \(x_i\) , \(y_j\)
    complex(dp),allocatable,intent(out)::data(:,:) !! sample values \(z_{ij}\) 
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,Ny,dtype_dummy
    integer::Flunit
    logical::isthere

    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.2) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx,Ny
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(y)) allocate(y(Ny))
    read(Flunit) y
    if(.not.allocated(data)) allocate(data(Nx,Ny))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData2D_from_BinaryFile_dtype2
!--------------------------------------------------------------------------------------
  subroutine ReadData2D_from_BinaryFile_dtype3(Flname,x,y,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:),y(:) !! sample points \(x_i\) , \(y_j\)
    real(dp),allocatable,intent(out)::data(:,:,:) !! sample values \(z_{ij}\) 
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,Ny,size1,dtype_dummy
    integer::Flunit
    logical::isthere

    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.3) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx,Ny,size1
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(y)) allocate(y(Ny))
    read(Flunit) y
    if(.not.allocated(data)) allocate(data(size1,Nx,Ny))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData2D_from_BinaryFile_dtype3
!--------------------------------------------------------------------------------------
  subroutine ReadData2D_from_BinaryFile_dtype4(Flname,x,y,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:),y(:) !! sample points \(x_i\) , \(y_j\)
    complex(dp),allocatable,intent(out)::data(:,:,:) !! sample values \(z_{ij}\) 
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,Ny,size1,dtype_dummy
    integer::Flunit
    logical::isthere

    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.4) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx,Ny,size1
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(y)) allocate(y(Ny))
    read(Flunit) y
    if(.not.allocated(data)) allocate(data(size1,Nx,Ny))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData2D_from_BinaryFile_dtype4
!--------------------------------------------------------------------------------------  
  subroutine ReadData2D_from_BinaryFile_dtype5(Flname,x,y,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:),y(:) !! sample points \(x_i\) , \(y_j\)
    real(dp),allocatable,intent(out)::data(:,:,:,:) !! sample values \(z_{ij}\) 
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,Ny,size1,size2,dtype_dummy
    integer::Flunit
    logical::isthere

    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.5) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx,Ny,size1,size2
    if(.not.allocated(x)) allocate(x(Nx))
    read(Flunit) x
    if(.not.allocated(y)) allocate(y(Ny))
    read(Flunit) y
    if(.not.allocated(data)) allocate(data(size1,size2,Nx,Ny))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData2D_from_BinaryFile_dtype5
!--------------------------------------------------------------------------------------
  subroutine ReadData2D_from_BinaryFile_dtype6(Flname,x,y,data,iflag)
    character(len=*),intent(in)::Flname !! file name for reading data
    real(dp),allocatable,intent(out)::x(:),y(:) !! sample points \(x_i\) , \(y_j\)
    complex(dp),allocatable,intent(out)::data(:,:,:,:) !! sample values \(z_{ij}\) 
    integer,intent(out)::iflag !! status flag. `iflag = -1`: file does not exist, `iflag = 1`: wrong data type
    integer::Nx,Ny,size1,size2,dtype_dummy
    integer::Flunit
    logical::isthere

    inquire(file=trim(Flname), exist=isthere)
    if(.not.isthere) then
       iflag=-1
       return
    end if
    
    open(newunit(Flunit),form='UNFORMATTED', file=trim(Flname))
    read(Flunit) dtype_dummy
    if(dtype_dummy.ne.6) then
       iflag=1
       close(FlUnit)
       return
    end if
    read(Flunit) Nx,Ny,size1,size2
    if(.not.allocated(x)) allocate(x(Nx))   
    read(Flunit) x
    if(.not.allocated(y)) allocate(y(Ny))
    read(Flunit) y
    if(.not.allocated(data)) allocate(data(size1,size2,Nx,Ny))
    read(Flunit) data

    close(FlUnit)

  end subroutine ReadData2D_from_BinaryFile_dtype6
!--------------------------------------------------------------------------------------
  
!======================================================================================    
end module scitools_binutils

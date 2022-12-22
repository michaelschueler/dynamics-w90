cmake_minimum_required(VERSION 3.0)
project("DYNAMICS-W90" Fortran)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


SET(CMAKE_Fortran_MODULE_DIRECTORY
   ${CMAKE_BINARY_DIR}/mod CACHE PATH "Single Directory for all fortran modules."
)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/exe)

option(wannier_lib "Building Wannier library" OFF)
if(wannier_lib)
    message(STATUS "Building Wannier library")
endif(wannier_lib)


# ~~ Commenly shared libs ~~
set(COMMON_LIBS "")

# ~~ OpenMP + SIMD ~~
find_package(OpenMP REQUIRED)
option(omp "Build with OpenMP support" ON)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
option(simd "Build with OpenMP SIMD support" OFF)
if (simd)
   message(STATUS "Building with SIMD support")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DOMPSIMD")
endif (simd)

# Require MPI for this project:
if(mpi)
   message(STATUS "Building with MPI")
   find_package(MPI REQUIRED)
   set(CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
   set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_LINK_FLAGS})
   include_directories(${MPI_INCLUDE_PATH})
   list(APPEND COMMON_LIBS ${MPI_LIBRARIES})
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMPI")
endif(mpi)

# ~~ Add SciTools ~~
if(scitools_prefix)
    find_library(SCITOOLS_LIB
        NAMES scitools
        HINTS "${scitools_prefix}/lib"
    )
   include_directories(${scitools_prefix}/include/scitools/)
   link_directories(${scitools_prefix}/lib)
   message(STATUS "SciTools library: ${SCITOOLS_LIB}")
endif(scitools_prefix)

# ~~ Add HDF5 ~~
if (hdf5)
  message(STATUS "Building with HDF5")
  find_library(HDF5_LIB hdf5)
  find_library(HDF5FOR_LIB hdf5_fortran)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DWITHHDF5")
  list(APPEND COMMON_LIBS ${HDF5_LIB})
  list(APPEND COMMON_LIBS ${HDF5FOR_LIB})
endif (hdf5)

# ~~ Add BLAS/LAPACK ~~
if(lapackblas_libraries)
   list(APPEND COMMON_LIBS ${lapackblas_libraries})
   message(STATUS LAPACK/BLAS: ${lapackblas_libraries})
else(lapackblas_libraries)
   message("Warning: LAPACK/BLAS not specified!")
endif(lapackblas_libraries)

# ~~ Add modern quadrature ~~
if (modernquad)
  message(STATUS "Building with modern quadrature")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMODERNQUAD")
endif (modernquad)

# ~~ Test programs ~~
option(compile_test "Compiling test programs" OFF)
if (compile_test)
  message(STATUS "Building test programs")
endif (compile_test)

# ~~ CMake Configuration ~~
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ~~ Libraries ~~
link_directories(${CMAKE_LIBRARY_PATH})

# ~~ Paths and Subdirs ~~
include_directories(${CMAKE_INCLUDE_PATH})

if(scitools_prefix)

else(scitools_prefix)
    message(STATUS "Path to SciTools not specified. Building internal copy.")
    add_subdirectory(scitools)
endif(scitools_prefix)

add_subdirectory(src)


# ~~ Include Guard ~~
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing, that CMake has configured at least once.")

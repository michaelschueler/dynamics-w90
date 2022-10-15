include_directories(
  SYSTEM
  )

set(OBJECT_FILES
		binutils.f90
		debug.F90
		def.f90
		linalg.f90
		time.f90
		utils.F90
		sorting.f90
)

if(hdf5)
	list(APPEND OBJECT_FILES hdf5_utils.f90)
endif(hdf5)

if(mpi)
	list(APPEND OBJECT_FILES array1d_dist.F90)
endif(mpi)


# ~~ Installation targets ~~
add_library(basic OBJECT ${OBJECT_FILES})
target_link_libraries(basic ${COMMON_LIBS})
# install(TARGETS basic ARCHIVE DESTINATION lib)

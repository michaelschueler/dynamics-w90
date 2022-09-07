include_directories(
  SYSTEM
  )

set(OBJECT_FILES
		Mdebug.F90
		Mdef.f90
		Mlinalg.f90
		Mtime.f90
		Mutils.f90
		Msorting.f90
)

if(hdf5)
	list(APPEND OBJECT_FILES Mhdf5_utils.f90)
endif(hdf5)

if(mpi)
	list(APPEND OBJECT_FILES Marray1d_dist.F90)
endif(mpi)


# ~~ Installation targets ~~
add_library(tools STATIC ${OBJECT_FILES})
target_link_libraries(tools ${COMMON_LIBS})
# install(TARGETS tools ARCHIVE DESTINATION lib)

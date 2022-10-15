
set(OBJECT_FILES
	basic/binutils.f90
	basic/debug.F90
	basic/def.f90
	basic/linalg.f90
	basic/sorting.f90
	basic/time.f90
	basic/utils.F90
	math/bsplines.f90
	math/fornberg.f90
	math/filon.F90
	math/gausslegendre.F90
	math/gregory.f90
	math/lebedev_weights.F90
	math/lebedev_quad.F90
	math/matrixexp.F90
	math/quadpack.f90
	math/quadrature.F90
	math/root.f90
	math/rungekutta.f90
	math/special.F90
	math/vector_bsplines.F90
	phys/evol.F90
	phys/laserpulse.F90
)

if (hdf5)
   list(APPEND OBJECT_FILES
	   basic/hdf5_utils.f90
   )
endif (hdf5)

if (mpi)
   list(APPEND OBJECT_FILES
	   basic/array1d_dist.F90
   )
endif (mpi)

if (enable_static)
	add_library(scitools STATIC ${OBJECT_FILES})
	install(TARGETS scitools ARCHIVE DESTINATION lib)
else (enable_static)
	add_library(scitools SHARED ${OBJECT_FILES})
	install(TARGETS scitools LIBRARY DESTINATION lib)
endif (enable_static)

target_link_libraries(scitools ${COMMON_LIBS})
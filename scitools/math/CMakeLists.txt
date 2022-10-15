include_directories(
  SYSTEM
 )

set(OBJECT_FILES
	bsplines.f90
	fornberg.f90
	filon.F90
	gausslegendre.F90
	gregory.f90
	lebedev_quad.F90
	lebedev_weights.F90
	matrixexp.F90
	quadpack.f90
	quadrature.F90
	root.f90
	rungekutta.f90
	special.F90
	vector_bsplines.F90
)


# ~~ Installation targets ~~
add_library(math OBJECT ${OBJECT_FILES})
target_link_libraries(math basic ${COMMON_LIBS})
add_dependencies(math basic)
# install(TARGETS math ARCHIVE DESTINATION lib)


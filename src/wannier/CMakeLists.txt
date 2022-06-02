include_directories(
  SYSTEM
 )

set(OBJECT_FILES
	Mlatt.f90
	Mham_w90.F90
	Mlatt_utils.f90
	Mlatt_kpts.F90
	Mwann_compress.f90
	Mwann_dyn.F90
	Mwannier_calc.F90
)

# ~~ Installation targets ~~
add_library(wannier STATIC ${OBJECT_FILES})
target_link_libraries(wannier math tools ${COMMON_LIBS})
add_dependencies(wannier math tools)
install(TARGETS wannier ARCHIVE DESTINATION lib)


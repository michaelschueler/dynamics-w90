include_directories(
  SYSTEM
 )

set(OBJECT_FILES
	Mham_w90.f90
	Mlatt_utils.f90
	Mwannier_orbitals.F90
	Mwannier_soc.F90
)

# ~~ Installation targets ~~
if(enable_shared)
	add_library(wannier SHARED ${OBJECT_FILES})
	target_link_libraries(wannier math tools ${COMMON_LIBS})
	add_dependencies(wannier math tools)
	install(TARGETS wannier LIBRARY DESTINATION lib)
else(enable_shared)
	add_library(wannier STATIC ${OBJECT_FILES})
	target_link_libraries(wannier math tools ${COMMON_LIBS})
	add_dependencies(wannier math tools)
	install(TARGETS wannier ARCHIVE DESTINATION lib)
endif(enable_shared)


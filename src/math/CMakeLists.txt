include_directories(
  SYSTEM
 )

set(OBJECT_FILES
	Mbsplines.f90
	Mfornberg.f90
	Mstencil.f90
	Mfourierintegral.F90
	Mlegendre.F90
	Mpolyfit.f90
	Mquadpack.f90
	Mroot.f90
  Mrungekutta.f90
  Msparse.f90
  Mintegration.f90
  Mspecial.F90
  Mgausslegendre.F90
  Mquadrature.F90
)


# ~~ Installation targets ~~
if(enable_shared)
	add_library(math SHARED ${OBJECT_FILES})
	target_link_libraries(math tools ${COMMON_LIBS})
	add_dependencies(math tools)
	install(TARGETS math LIBRARY DESTINATION lib)
else(enable_shared)
	add_library(math STATIC ${OBJECT_FILES})
	target_link_libraries(math tools ${COMMON_LIBS})
	add_dependencies(math tools)
	install(TARGETS math ARCHIVE DESTINATION lib)
endif(enable_shared)


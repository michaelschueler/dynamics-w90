include_directories(
  SYSTEM
  )

set(OBJECT_FILES
	COUL90.f
	Mangcoeff.F90
	Matomic.F90
	Mmatrix_elements.F90
	Mradialwf.F90
	Mscattwf.F90
	Mwignerd.f90
	Mpes_intensity.F90
)

add_library(pes STATIC ${OBJECT_FILES})
target_link_libraries(pes tools math wannier ${COMMON_LIBS})
add_dependencies(pes tools math wannier )
# install(TARGETS pes ARCHIVE DESTINATION lib)

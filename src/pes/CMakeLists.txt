include_directories(
  SYSTEM
  )

set(OBJECT_FILES
	Mangcoeff.F90
	Matomic.F90
	Mmatrix_elements.F90
	Mradialwf.F90
	Mscattwf.F90
	Mwignerd.f90
)

add_library(pes STATIC ${OBJECT_FILES})
target_link_libraries(pes tools math ${COMMON_LIBS})
add_dependencies(pes tools math )
# install(TARGETS pes ARCHIVE DESTINATION lib)

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
	Mradialintegral.F90
	Mpes_intensity.F90
)

add_library(pes OBJECT ${OBJECT_FILES})
target_link_libraries(pes basic wannier scitools ${COMMON_LIBS})
if(scitools_prefix)
    add_dependencies(pes basic wannier )
else(scitools_prefix)
    add_dependencies(pes basic wannier scitools)
endif(scitools_prefix)

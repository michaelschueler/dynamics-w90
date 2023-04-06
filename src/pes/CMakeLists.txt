include_directories(
  SYSTEM
  )

set(OBJECT_FILES
	COUL90.f
	pes_angcoeff.F90
	pes_atomic.F90
	pes_matel.F90
	pes_radialwf.F90
	pes_scattwf.F90
	pes_wignerd.f90
	pes_radialintegral.F90
	pes_main.F90
)

add_library(pes STATIC ${OBJECT_FILES})
target_link_libraries(pes basic wannier scitools ${COMMON_LIBS})
if(scitools_prefix)
    add_dependencies(pes basic wannier )
else(scitools_prefix)
    add_dependencies(pes basic wannier scitools)
endif(scitools_prefix)

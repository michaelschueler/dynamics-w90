include_directories(
  SYSTEM
  )

set(OBJECT_FILES
		Mdebug.F90
		Mconstants.f90
		Mbinutils.f90
		Mdef.f90
		Mgrid1d.f90
		Mhdf5_utils.f90
		Mlinalg.f90
		Mmesh.f90
		Mtime.f90
		Munits.f90
		Mutils.f90
		Msorting.f90
)


# ~~ Installation targets ~~
if(enable_shared)
	add_library(tools SHARED ${OBJECT_FILES})
	target_link_libraries(tools ${COMMON_LIBS})
	install(TARGETS tools LIBRARY DESTINATION lib)
else(enable_shared)
	add_library(tools STATIC ${OBJECT_FILES})
	target_link_libraries(tools ${COMMON_LIBS})
	install(TARGETS tools ARCHIVE DESTINATION lib)
endif(enable_shared)

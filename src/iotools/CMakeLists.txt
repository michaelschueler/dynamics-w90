include_directories(
  SYSTEM
  )

set(OBJECT_FILES
    io_params.F90
		io_obs.F90
		io_hamiltonian.F90
    io_orbitals.F90
)

add_library(iotools OBJECT ${OBJECT_FILES})
target_link_libraries(iotools basic wannier scitools ${COMMON_LIBS})
if(scitools_prefix)
  add_dependencies(iotools basic wannier )
else(scitools_prefix)
  add_dependencies(iotools basic wannier scitools)
endif(scitools_prefix)

# install(TARGETS iotools ARCHIVE DESTINATION lib)

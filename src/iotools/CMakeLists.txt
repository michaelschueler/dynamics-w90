include_directories(
  SYSTEM
  )

set(OBJECT_FILES
    Mio_params.F90
		Mio_obs.F90
		Mio_hamiltonian.F90
)

add_library(iotools STATIC ${OBJECT_FILES})
target_link_libraries(iotools tools math wannier ${COMMON_LIBS})
add_dependencies(iotools tools math wannier )
install(TARGETS iotools ARCHIVE DESTINATION lib)

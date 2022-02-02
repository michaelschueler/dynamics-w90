include_directories(
  SYSTEM

)

# A list of executables to build
set( project_EXECUTABLES
     latt_wann.x
)

# mpi excecutables
if(mpi)
	set(mpi_EXECUTABLES
        latt_wann_mpi.x
	)
endif(mpi)

set( EXE_latt_wann.x_SOURCES Mlatt_wann.F90 latt_wann.F90)

foreach( exe ${project_EXECUTABLES} )
    # declare the target...
    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
    add_dependencies(${exe} tools math evol wannier)
    target_link_libraries(${exe} tools tools math evol wannier ${COMMON_LIBS})
endforeach()

# ---- MPI programs ----
set( EXE_latt_wann_mpi.x_SOURCES Mlatt_wann_mpi.F90 latt_wann_mpi.F90)

foreach( exe ${mpi_EXECUTABLES} )
    # declare the target...
    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
    add_dependencies(${exe} tools math evol wannier)
    target_link_libraries(${exe} tools tools math evol wannier ${COMMON_LIBS})
endforeach()

include_directories(
  SYSTEM

)

# A list of executables to build
set( project_EXECUTABLES
     wann_evol.x
     wann_calc.x
     wann_prune.x
     wann_soc.x
     arpes.x
)

# mpi excecutables
if(mpi)
	set(mpi_EXECUTABLES
        wann_evol_mpi.x
        arpes_mpi.x
	)
endif(mpi)

set( EXE_wann_evol.x_SOURCES Mwann_evol.F90 wann_evol.F90)
set( EXE_wann_calc.x_SOURCES Mwannier_calc.F90 wann_calc.F90)
set( EXE_wann_prune.x_SOURCES wann_prune.F90)
set( EXE_wann_soc.x_SOURCES wann_soc.F90)
set( EXE_arpes.x_SOURCES Marpes_calc.F90 arpes.F90)

foreach( exe ${project_EXECUTABLES} )
    # declare the target...
    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
    add_dependencies(${exe} tools math iotools evol wannier pes)
    target_link_libraries(${exe} tools math iotools evol wannier pes ${COMMON_LIBS})
endforeach()

# ---- MPI programs ----
set( EXE_wann_evol_mpi.x_SOURCES Mwann_evol_mpi.F90 wann_evol_mpi.F90)
set( EXE_arpes_mpi.x_SOURCES Marpes_calc_mpi.F90 arpes_mpi.F90)

foreach( exe ${mpi_EXECUTABLES} )
    # declare the target...
    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
    add_dependencies(${exe} tools math iotools evol wannier pes)
    target_link_libraries(${exe} tools math iotools evol wannier pes ${COMMON_LIBS})
endforeach()

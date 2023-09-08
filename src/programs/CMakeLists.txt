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

if(mpi)
    set(mpi_EXECUTABLES
        wann_evol_mpi.x
        arpes_mpi.x
    )
endif(mpi)

if(fftw_mpi)
    list(APPEND mpi_EXECUTABLES
        test_fftw_mpi.x
    )    
    set( EXE_test_fftw_mpi.x_SOURCES test_fftw_mpi.F90)
endif(fftw_mpi)


set( EXE_wann_evol.x_SOURCES Mwann_evol.F90 wann_evol.F90)
set( EXE_wann_calc.x_SOURCES Mwannier_calc.F90 wann_calc.F90)
set( EXE_wann_prune.x_SOURCES wann_prune.F90)
set( EXE_wann_soc.x_SOURCES wann_soc.F90)
set( EXE_arpes.x_SOURCES Marpes_calc.F90 arpes.F90)



foreach( exe ${project_EXECUTABLES} )
    # declare the target...
    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
    if(scitools_prefix)
        add_dependencies(${exe} basic iotools wannier pes)
    else(scitools_prefix)
        add_dependencies(${exe} basic iotools wannier pes scitools)
    endif(scitools_prefix)
    target_link_libraries(${exe} basic iotools wannier pes scitools ${COMMON_LIBS})
    
    if(CMAKE_INSTALL_PREFIX)
        INSTALL( TARGETS ${exe} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)   
    else(CMAKE_INSTALL_PREFIX)
        INSTALL( TARGETS ${exe} DESTINATION ${PROJECT_SOURCE_DIR}/bin)
    endif(CMAKE_INSTALL_PREFIX)
    
endforeach()

# ---- MPI programs ----
set( EXE_wann_evol_mpi.x_SOURCES Mwann_evol_mpi.F90 wann_evol_mpi.F90)
set( EXE_arpes_mpi.x_SOURCES Marpes_calc_mpi.F90 arpes_mpi.F90)

foreach( exe ${mpi_EXECUTABLES} )
    # declare the target...
    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
    if(scitools_prefix)
        add_dependencies(${exe} basic iotools wannier pes)
    else(scitools_prefix)
        add_dependencies(${exe} basic iotools wannier pes scitools)
    endif(scitools_prefix)
    target_link_libraries(${exe} basic iotools wannier pes scitools ${COMMON_LIBS})
    
    if(CMAKE_INSTALL_PREFIX)
        INSTALL( TARGETS ${exe} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)   
    else(CMAKE_INSTALL_PREFIX)
        INSTALL( TARGETS ${exe} DESTINATION ${PROJECT_SOURCE_DIR}/bin)
    endif(CMAKE_INSTALL_PREFIX)
    
endforeach()






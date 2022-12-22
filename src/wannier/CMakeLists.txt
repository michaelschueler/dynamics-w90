include_directories(
  SYSTEM
 )

set(OBJECT_FILES
    wan_read_xyz.F90
    wan_hamiltonian.F90
    wan_latt_utils.f90
    wan_latt_kpts.F90
    wan_compress.F90
    wan_dynamics.F90
    wan_slab.F90
    wan_soc.F90
    wan_orbitals.F90
)

# ~~ Installation targets ~~
add_library(wannier STATIC ${OBJECT_FILES})
target_link_libraries(wannier basic scitools ${COMMON_LIBS})

if(scitools_prefix)
    add_dependencies(wannier basic)
else(scitools_prefix)
    add_dependencies(wannier basic scitools)
endif(scitools_prefix)
# install(TARGETS wannier ARCHIVE DESTINATION lib)

if(wannier_lib)
    set(LIB_FILES
        wan_read_xyz.F90
        wan_hamiltonian.F90
        wan_latt_utils.f90    
    )
    add_library(hamw90 SHARED ${LIB_FILES})
    target_link_libraries(hamw90 basic scitools ${COMMON_LIBS})

    if(scitools_prefix)
        add_dependencies(hamw90 basic)
    else(scitools_prefix)
        add_dependencies(hamw90 basic scitools)
    endif(scitools_prefix)

    install(TARGETS hamw90 LIBRARY DESTINATION lib)

    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/wan_read_xyz.mod DESTINATION include/hamw90)
    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/wan_hamiltonian.mod DESTINATION include/hamw90)
    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/wan_latt_utils.mod DESTINATION include/hamw90)
endif(wannier_lib)
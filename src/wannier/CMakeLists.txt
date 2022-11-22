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
add_library(wannier OBJECT ${OBJECT_FILES})
target_link_libraries(wannier basic scitools ${COMMON_LIBS})

if(scitools_prefix)
    add_dependencies(wannier basic)
else(scitools_prefix)
    add_dependencies(wannier basic scitools)
endif(scitools_prefix)
# install(TARGETS wannier ARCHIVE DESTINATION lib)


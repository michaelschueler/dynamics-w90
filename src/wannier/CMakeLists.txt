include_directories(
  SYSTEM
 )

set(OBJECT_FILES
    Mlatt.f90
    Mread_xyz.F90
    Mham_w90.F90
    Mlatt_utils.f90
    Mlatt_kpts.F90
    Mwann_compress.F90
    Mwann_dyn.F90
    Mwann_slab.F90
    Mwann_soc.F90
    Mwannier_orbitals.F90
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


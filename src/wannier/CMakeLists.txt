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
add_library(wannier STATIC ${OBJECT_FILES})
target_link_libraries(wannier math tools evol ${COMMON_LIBS})
add_dependencies(wannier math tools evol)
# install(TARGETS wannier ARCHIVE DESTINATION lib)

if(wannier_lib)
  set(LIB_FILES
    Mlatt_utils.f90
    Mread_xyz.F90
    Mham_w90.F90
  )

  add_library(ham_w90 SHARED ${LIB_FILES})
  install(TARGETS ham_w90 LIBRARY DESTINATION lib)

  target_link_libraries(ham_w90 tools math)
  add_dependencies(ham_w90 tools math)

  install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/Mham_w90.mod DESTINATION include/ham_w90)

endif(wannier_lib)

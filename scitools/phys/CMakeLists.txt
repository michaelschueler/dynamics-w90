include_directories(
  SYSTEM
 )

set(OBJECT_FILES
	evol.F90
	laserpulse.F90
)


# ~~ Installation targets ~~
add_library(phys STATIC ${OBJECT_FILES})
target_link_libraries(phys basic math ${COMMON_LIBS})
add_dependencies(phys basic math )
# install(TARGETS math ARCHIVE DESTINATION lib)


include_directories(
  SYSTEM
 )

set(OBJECT_FILES
	Mbsplines.f90
	Mintegration.f90
	Mquadpack.f90
	Mroot.f90
	Mrungekutta.f90
)


# ~~ Installation targets ~~
add_library(math STATIC ${OBJECT_FILES})
target_link_libraries(math tools ${COMMON_LIBS})
add_dependencies(math tools)
install(TARGETS math ARCHIVE DESTINATION lib)


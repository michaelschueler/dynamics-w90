include_directories(
  SYSTEM
  )

set(OBJECT_FILES
		Mmatrixexp.F90
		Mevol.F90
		Mlaserpulse.F90
)

add_library(evol STATIC ${OBJECT_FILES})
target_link_libraries(evol tools math ${COMMON_LIBS})
add_dependencies(evol tools math )
install(TARGETS evol ARCHIVE DESTINATION lib)

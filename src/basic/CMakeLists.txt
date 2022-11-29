include_directories(
  SYSTEM
  )

set(OBJECT_FILES
		Mdebug.F90)


# ~~ Installation targets ~~
add_library(basic STATIC ${OBJECT_FILES})
target_link_libraries(basic scitools ${COMMON_LIBS})
if(scitools_prefix)

else(scitools_prefix)
  add_dependencies(basic scitools)
endif(scitools_prefix)


# install(TARGETS basic ARCHIVE DESTINATION lib)

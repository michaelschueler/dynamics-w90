include_directories(
  SYSTEM
  )

set(OBJECT_FILES
		#Mcoeffs.f90
		Mintegration.f90
		Mmatrixexp.F90
		Mevol.f90
		Mmfevol.F90
		Mtdarpes_gkba.F90
		Mtdarpes_kbe.F90
		Megkba.f90
		Mspectral_pade.f90
		Mcntr_array.F90
		Mmatsubara_legendre.F90
		Mpack_green.F90
		Mfourier_timeslice.F90
)

if(mpi)
	list(APPEND OBJECT_FILES Mtdarpes_kbe_mpi.F90)
endif(mpi)

# ~~ Installation targets ~~
if(enable_shared)
	add_library(evol SHARED ${OBJECT_FILES})
	if(mpi)
		target_link_libraries(evol tools steady mpitools kbkit ${COMMON_LIBS})
		add_dependencies(evol tools steady mpitools)
	else(mpi)
		target_link_libraries(evol tools steady kbkit ${COMMON_LIBS})
		add_dependencies(evol tools steady)
	endif(mpi)
	install(TARGETS evol
	        LIBRARY DESTINATION lib)
else(enable_shared)
	add_library(evol STATIC ${OBJECT_FILES})
	if(mpi)
		target_link_libraries(evol tools steady mpitools kbkit ${COMMON_LIBS})
		add_dependencies(evol tools steady mpitools)
	else(mpi)
		target_link_libraries(evol tools steady kbkit ${COMMON_LIBS})
		add_dependencies(evol tools steady)
	endif(mpi)
	install(TARGETS evol ARCHIVE DESTINATION lib)
endif(enable_shared)


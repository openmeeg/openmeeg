if (USE_MKL)
	find_package(MKL)
	if (MKL_FOUND)
		include_directories(${MKL_INCLUDE_DIR})
		set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
        get_filename_component(MKL_DLL_DIR "${MKL_LIBRARIES}" DIRECTORY) 
        #message(${LAPACK_LIBRARIES}) # for debug
		if(UNIX AND NOT APPLE) # MKL on linux requires to link with the pthread library
			set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} pthread)
		endif()
	else()
		message(FATAL_ERROR "MKL not found. Please set environment variable MKLDIR")
	endif()
endif()

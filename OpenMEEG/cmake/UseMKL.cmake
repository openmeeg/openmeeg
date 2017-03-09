
if (USE_MKL)
    find_package(MKL ${REQUIRED}) # e.g REQUIRED can be REQUIRED or QUIET
    if (MKL_FOUND)
        include_directories(${MKL_INCLUDE_DIR})
        set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
        get_filename_component(MKL_DLL_DIR "${MKL_LIBRARIES}" DIRECTORY)
        if (UNIX AND NOT APPLE) # MKL on linux requires to link with the pthread library
            set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} pthread)
        endif()
    endif()
endif()

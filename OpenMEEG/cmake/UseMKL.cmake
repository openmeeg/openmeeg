if (USE_MKL)
    # user defined options for MKL
    option(MKL_USE_parallel "Use MKL parallel" True)
    option(MKL_USE_sdl "Single Dynamic Library or static/dynamic" False)
    set(MKL_USE_interface "lp64" CACHE STRING "for Intel(R)64 compatible arch: ilp64/lp64 or for ia32 arch: cdecl/stdcall")

    find_package(MKL ${FIND_MODE}) # e.g FIND_MODE can be REQUIRED or QUIET
    if (MKL_FOUND)
        include_directories(${MKL_INCLUDE_DIR})
        set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
        set(LAPACK_DEFINITIONS ${MKL_DEFINITIONS})
        list(APPEND OpenMEEG_DEPENDENCIES MKL)
        set(CMAKE_MSVCIDE_RUN_PATH ${MKL_LIBRARY_DIR} ${CMAKE_MSVCIDE_RUN_PATH})
    elseif(NOT FIND_MODE STREQUAL "REQUIRED") # case we were in Auto Mode
        unset(MKL_USE_parallel CACHE)
        unset(MKL_USE_sdl CACHE)
        unset(MKL_USE_interface CACHE)
        unset(MKL_ROOT_DIR CACHE)
    endif()
endif()

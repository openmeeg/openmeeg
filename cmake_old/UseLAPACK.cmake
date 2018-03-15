
if (USE_LAPACK)
    # first look for the compiled clapack
    find_package(LAPACK QUIET CONFIG NAMES clapack PATHS ${lapack_DIR})
    if (LAPACK_FOUND)
        set(LAPACK_INCLUDE_DIRS "${CLAPACK_INCLUDE_DIRS}")
        set(LAPACK_LIBRARY_DIRS "${CLAPACK_LIBRARY_DIRS}")
        set(LAPACK_LIBRARIES    "${CLAPACK_LIBRARIES}")
    else()
        # if not found, try to use the system lapack
        find_package(LAPACK QUIET MODULE)
    endif()
    if (LAPACK_FOUND)
        include_directories(${LAPACK_INCLUDE_DIRS})
        set(CMAKE_MSVCIDE_RUN_PATH ${LAPACK_LIBRARY_DIRS} ${CMAKE_MSVCIDE_RUN_PATH})
    endif()
endif()

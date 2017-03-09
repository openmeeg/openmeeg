
if (USE_LAPACK)
    # first look for the compiled clapack
    find_package(LAPACK QUIET CONFIG NAMES clapack PATHS ${lapack_DIR})
    if (LAPACK_FOUND)
        set(LAPACK_INCLUDE_DIRS "${CLAPACK_INCLUDE_DIRS}")
        set(LAPACK_LIBRARY_DIRS "${CLAPACK_LIBRARY_DIRS}")
        set(LAPACK_LIBRARIES    "${CLAPACK_LIBRARIES}")
        include_directories(${CLAPACK_LIBRARY_DIRS})
    else()
        # if not found, try to use the system lapack
        find_package(LAPACK QUIET MODULE)
    endif()
endif()

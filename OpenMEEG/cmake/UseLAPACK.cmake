
if (USE_LAPACK)
    # first look for the compiled clapack
    set(LAPACK_DIR ${lapack_DIR})
    find_package(LAPACK QUIET NAMES clapack PATHS ${LAPACK_DIR} CONFIG)
    if (LAPACK_FOUND)
        set(LAPACK_INCLUDE_DIRS "${CLAPACK_INCLUDE_DIRS}")
        set(LAPACK_LIBRARY_DIRS "${CLAPACK_LIBRARY_DIRS}")
        set(LAPACK_LIBRARIES    "${CLAPACK_LIBRARIES}")
        include_directories(${CLAPACK_LIBRARY_DIRS})
    else()
        # if not found, try to use the system lapack
        find_package(LAPACK ${REQUIRED} MODULE)
        if (NOT LAPACK_FOUND)
            message(FATAL_ERROR "Clapack was not built, please re-run cmake 
                super-project with \"-DBLASLAPACK_IMPLEMENTATION=LAPACK\" to 
                build/install it.")
        endif()
    endif()
endif()

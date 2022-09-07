set(BLA_SIZEOF_INTEGER 4)
set(BLA_IMPLEMENTATION "OpenBLAS" CACHE STRING "BLAS/LAPACK implementation")
set_property(CACHE BLA_IMPLEMENTATION PROPERTY STRINGS "OpenBLAS" "mkl" "mkl-findblas")

if (BLA_IMPLEMENTATION STREQUAL "mkl-findblas")

    # Just use the findblas interface
    message(STATUS "Using BLA_IMPLEMENTATION=mkl-findblas")
    set(USE_MKL ON)
    set(HAVE_BLAS ON)
    set(HAVE_LAPACK ON)

elseif (BLA_IMPLEMENTATION STREQUAL "mkl")
    message(STATUS "Using BLA_IMPLEMENTATION=mkl")
    set(MKL_PARALLELISM "parallel" CACHE STRING "MKL parallelism")
    set_property(CACHE MKL_PARALLELISM PROPERTY STRINGS "parallel" "sequential" "sdl")

    set(MKL_PARALLEL_SUFFIX)
    if (${MKL_PARALLELISM} STREQUAL "sequential")
        set(MKL_PARALLEL_SUFFIX "_seq")
        set(MKL_THREADING sequential)
    elseif (${MKL_PARALLELISM} STREQUAL "sdl")
        if (BLA_STATIC)
            message(FATAL_ERROR "mkl sdl mode incompatible with static linking")
        endif()
        set(MKL_PARALLEL_SUFFIX "_dyn")
        set(MKL_LINK sdl)
    endif()

    if (BLA_STATIC)
        set(MKL_LINK static)
    endif()

    set(MKL_ARCH "intel64" CACHE STRING "MKL archirecture name")
    set_property(CACHE MKL_ARCH PROPERTY STRINGS "ia32" "intel64")

    set(MKL_ARCH_BITS 32)
    unset(MKL_INT_MODEL)
    if (${MKL_ARCH} STREQUAL "intel64")
        set(MKL_INTERFACE "lp64" CACHE STRING "MKL interface")
        set_property(CACHE MKL_INTERFACE PROPERTY STRINGS "lp64" "ilp64")
        set(MKL_ARCH_BITS 64)
        set(MKL_INT_MODEL lp)
        if (${MKL_INTERFACE} STREQUAL "ilp64")
            set(MKL_INTERFACE_DEF MKL_ILP64)
            set(MKL_INT_MODEL ilp)
            set(BLA_SIZEOF_INTEGER 8)
        endif()
    else()
        unset(MKL_INTERFACE)
        unset(MKL_INTERFACE_DEF)
        set(MKL_ARCH_BITS 32)
    endif()

    # MKL stuff: This is needed because FindBLAS does not define includes variable.
    # It also allows to avoid the need to define global shell variables (with a small hack).
    # So this simplifies the use of MKL.

    find_package(MKL REQUIRED)

    # This is a hack to help FindBLAS finding MKL without needing source the mkl script
    # setting global variables.

    set(ENV{MKLROOT} ${MKL_ROOT})
    if (${OMP_LIBRARY})
        get_filename_component(OMP_LIBRARY_DIR ${OMP_LIBRARY} DIRECTORY)
        set(ENV{LD_LIBRARY_PATH} ${OMP_LIBRARY_DIR})
        set(_libs ${OMP_LIBRARY_DIR} ${MKL_LIBRARIES})
        set(ENV{LIBRARY_PATH} "${_libs}")
    endif()

    # For some reason ilp version of MKL does not work. TODO.
    # So for the time being, we force lp mode.

    set(BLA_VENDOR Intel10_${MKL_ARCH_BITS}${MKL_INT_MODEL}${MKL_PARALLEL_SUFFIX})
    set(USE_MKL ON)
    set(HAVE_BLAS ON)
    set(HAVE_LAPACK ON)
    set(BLA_INCLUDE_DIR ${MKL_INCLUDE})

elseif (BLA_IMPLEMENTATION STREQUAL "OpenBLAS")

    message(STATUS "Using BLA_IMPLEMENTATION=OpenBLAS")
    set(USE_OPENBLAS ON)
    set(HAVE_BLAS ON)
    set(HAVE_LAPACK ON)
    set(BLA_VENDOR ${BLA_IMPLEMENTATION})

else()

    message(STATUS "Using no BLAS implementation")

endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Add targets for compatibility with older cmake versions < 3.18

if (NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    if (BLAS_LIBRARIES)
        set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
    endif()
    if (BLAS_LINKER_FLAGS)
        set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_OPTIONS ${BLAS_LINKER_FLAGS})
    endif()
endif()

if (NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)

    # Filter out redundant BLAS info and replace with the BLAS target

    set(_lapack_libs ${LAPACK_LIBRARIES})
    set(_lapack_flags ${LAPACK_LINKER_FLAGS})
    if (TARGET BLAS::BLAS)
        if (_lapack_libs AND BLAS_LIBRARIES)
            foreach(_blas_lib IN LISTS BLAS_LIBRARIES)
                list(REMOVE_ITEM _lapack_libs "${_blas_lib}")
            endforeach()
        endif()
        if (_lapack_flags AND BLAS_LINKER_FLAGS)
            foreach(_blas_flag IN LISTS BLAS_LINKER_FLAGS)
                list(REMOVE_ITEM _lapack_flags "${_blas_flag}")
            endforeach()
        endif()
        list(APPEND _lapack_libs BLAS::BLAS)
    endif()
    if (_lapack_libs)
        set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${_lapack_libs}")
    endif()
    if (_lapack_flags)
        set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_OPTIONS "${_lapack_flags}")
    endif()
endif()

# OpenBLAS may or may not include lapacke.
# Check which version is used.

set(CMAKE_REQUIRED_LIBRARIES LAPACK::LAPACK BLAS::BLAS)
check_function_exists(LAPACKE_dlange LAPACKE_WORKS)
mark_as_advanced(LAPACKE_WORKS)
if (NOT LAPACKE_WORKS)
    find_library(LAPACKE lapacke REQUIRED)
    list(PREPEND _lapack_libs ${LAPACKE})
    set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${_lapack_libs}")
endif()

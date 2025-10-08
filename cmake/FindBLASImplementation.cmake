include(CheckFunctionExists)
include(CheckSymbolExists)

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
    if (USE_SCIPY_OPENBLAS)
        message(STATUS "  Using SciPy OpenBLAS variant with scipy_ prefix")
        # In theory this should work, but it doesn't for some reason:
        #
        # find_package(PkgConfig REQUIRED)
        # pkg_search_module(SP_OB REQUIRED scipy_openblas)
        # message(STATUS ".. Found OpenBLAS version ${SP_OB_VERSION} at ${SP_OB_LIBRARY_DIRS} and includes ${SP_OB_INCLUDE_DIRS}")
        # set(CMAKE_REQUIRED_INCLUDES ${SP_OB_INCLUDE_DIRS})
        # set(BLA_PREFER_PKGCONFIG ON)
        # set(BLA_PKGCONFIG_BLAS scipy_openblas)
        # set(BLA_PKGCONFIG_LAPACK scipy_openblas)
        #
        # So we take care of it manually:
        find_package(OpenBLAS REQUIRED CONFIG)
        message(STATUS "  Libraries:    ${OpenBLAS_LIBRARIES}")
        message(STATUS "  Include dirs: ${OpenBLAS_INCLUDE_DIRS}")
        message(STATUS "  Linker flags: ${OpenBLAS_LDFLAGS}")
        get_filename_component(OpenBLAS_DLL ${OpenBLAS_LIBRARIES} NAME)
        set(OpenBLAS_DLL "${OpenBLAS_DLL}.dll")
        mark_as_advanced(OpenBLAS_INCLUDE_DIRS OpenBLAS_LIBRARIES OpenBLAS_LDFLAGS OpenBLAS_DLL)
        set(CMAKE_REQUIRED_INCLUDES ${OpenBLAS_INCLUDE_DIRS})
        foreach (TARG BLAS LAPACK)
            add_library(${TARG}::${TARG} SHARED IMPORTED)
            set_target_properties(${TARG}::${TARG} PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${OpenBLAS_INCLUDE_DIRS}"
                LINK_FLAGS "${OpenBLAS_LDFLAGS}"
            )
            # If on Windows we need the implib and the actual DLL
            if (WIN32)
                set_target_properties(${TARG}::${TARG} PROPERTIES
                    IMPORTED_IMPLIB "${OpenBLAS_LIBRARIES}"
                    IMPORTED_LOCATION "${OpenBLAS_DLL}"
                )
            else()
                set_target_properties(${TARG}::${TARG} PROPERTIES
                    IMPORTED_LOCATION "${OpenBLAS_LIBRARIES}"
                )
            endif()
        endforeach()
    endif()

else()

    message(STATUS "Using no BLAS implementation")

endif()

if (NOT USE_SCIPY_OPENBLAS)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()

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
if (USE_SCIPY_OPENBLAS)
    set(CHECK_DLANGE scipy_LAPACKE_dlange)
else()
    set(CHECK_DLANGE LAPACKE_dlange)
endif()
mark_as_advanced(CHECK_DLANGE)
check_symbol_exists(${CHECK_DLANGE} "lapacke.h" HAVE_DLANGE)
mark_as_advanced(HAVE_DLANGE)
if (NOT HAVE_DLANGE)
    message(STATUS "${CHECK_DLANGE} not found!")
endif()
check_function_exists(${CHECK_DLANGE} LAPACKE_WORKS)
mark_as_advanced(LAPACKE_WORKS)
# TODO: unclear why this test doesn't work for scipy_openblas, so we skip it for now
if (NOT LAPACKE_WORKS AND NOT SCIPY_USE_OPENBLAS)
    find_library(LAPACKE lapacke REQUIRED)
    list(PREPEND _lapack_libs ${LAPACKE})
    set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${_lapack_libs}")
endif()

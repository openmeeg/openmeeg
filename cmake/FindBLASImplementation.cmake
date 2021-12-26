set(BLA_DEFINITIONS)
set(BLA_VENDOR "OpenBLAS" CACHE STRING "BLAS/LAPACK implementation")
set(BLAS_IMPLEMENTATION "OpenBLAS" CACHE STRING "BLAS/LAPACK implementation")
set_property(CACHE BLAS_IMPLEMENTATION PROPERTY STRINGS "OpenBlas" "mkl")

if (BLAS_IMPLEMENTATION STREQUAL "MKL")
    set(MKL_PARALLELISM "parallel" CACHE STRING "MKL parallelism")
    set_property(CACHE MKL_PARALLELISM PROPERTY STRINGS "parallel" "sequential" "sdl")

    set(MKL_PARALLEL_SUFFIX)
    if (${MKL_PARALLELISM} STREQUAL "seq")
        set(MKL_PARALLEL_SUFFIX "_seq")
        set(MKL_THREADING sequential)
    elif (${MKL_PARALLELISM} STREQUAL "sdl")
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
    set(MKL_ARCH_BITS 64)
    if (${MKL_ARCH} STREQUAL "ia32")
        set(MKL_ARCH_BITS 32)
    endif()

    # MKL stuff: This is needed because FindBLAS does not define includes variable.
    # It also allows to avoid the need to define global shell variables (with a small hack).
    # So this simplifies the use of MKL.

    set(MKL_INTERFACE_FULL intel_lp${MKL_ARCH_BITS})

    find_package(MKL REQUIRED)

    # This is a hack to help FindBLAS fiding MKL without needing source the mkl script
    # setting global variables.

    set(ENV{MKLROOT} ${MKL_ROOT})
    get_filename_component(OMP_LIBRARY_DIR ${OMP_LIBRARY} DIRECTORY)
    set(ENV{LD_LIBRARY_PATH} ${OMP_LIBRARY_DIR})

    # For some reason ilp version of MKL does not work. TODO.
    # So for the time being, we force lp mode.

    set(BLA_VENDOR Intel10_${MKL_ARCH_BITS}lp${MKL_PARALLEL_SUFFIX})
    set(BLA_DEFINITIONS USE_MKL HAVE_BLAS HAVE_LAPACK)
    set(BLA_INCLUDE_DIR ${MKL_INCLUDE})

elseif (BLAS_IMPLEMENTATION STREQUAL "OpenBLAS")

    set(BLA_DEFINITIONS USE_OPENBLAS HAVE_BLAS HAVE_LAPACK)
    set(BLA_VENDOR ${BLAS_IMPLEMENTATION})

endif()

set(BLA_SIZEOF_INTEGER 4)
find_package(BLAS)

if (NOT BLAS_FOUND)
    message(FATAL_ERROR "No BLAS library was found. Please define BLA_VENDOR appropriately.\
            For Intel MKL you can try to source the setvars.sh file.")
endif()

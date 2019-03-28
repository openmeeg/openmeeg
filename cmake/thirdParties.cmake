###############
# FindMKL stuff
###############

set(BLA_DEFINITIONS)
set(BLA_VENDOR "OpenBLAS" CACHE STRING "BLAS/LAPACK implementation")


if (BLA_VENDOR MATCHES Intel)
    if ("$ENV{MKLROOT}" STREQUAL "")
        message(FATAL_ERROR "MKLROOT is not set. Please source the Intel MKL mklvars.sh file.")
    endif()

    # user defined options for MKL
    option(MKL_USE_parallel "Use MKL parallel" True)
    option(MKL_USE_sdl "Single Dynamic Library or static/dynamic" False)
    set(MKL_USE_interface "lp64" CACHE STRING "for Intel(R)64 compatible arch: ilp64/lp64 or for ia32 arch: cdecl/stdcall")

    if (BLA_VENDOR MATCHES "_seq")
        set(MKL_USE_parallel OFF)
    else()
        set(MKL_USE_parallel ON)
    endif()

    find_package(MKL REQUIRED)

    if (MKL_FOUND)
        set(BLA_INCLUDE_DIR ${MKL_INCLUDE_DIR})
        set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
        set(BLA_DEFINITIONS USE_MKL HAVE_BLAS HAVE_LAPACK)
    endif()
elseif(BLA_VENDOR MATCHES OpenBLAS) # XXX: OpenBLAS should be set up using find_package(BLAS)
    find_package(OpenBLAS REQUIRED)
    set(BLA_INCLUDE_DIR ${OpenBLAS_INCLUDE_DIR})
    set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})
    set(BLA_DEFINITIONS USE_OPENBLAS HAVE_BLAS HAVE_LAPACK)
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    if (BLA_VENDOR MATCHES ATLAS)
        find_library(CBLAS_LIB NAMES cblas)
        set(BLAS_LIBRARIES "${BLAS_LIBRARIES};${CBLAS_LIB}")
        find_library(LAPACK_ATLAS_LIB NAMES lapack_atlas)
        set(LAPACK_LIBRARIES "${LAPACK_LIBRARIES};${LAPACK_ATLAS_LIB}")
        set(BLA_DEFINITIONS USE_ATLAS HAVE_BLAS HAVE_LAPACK)
    elseif(BLA_VENDOR MATCHES Apple)
        set(BLA_DEFINITIONS USE_VECLIB HAVE_BLAS HAVE_LAPACK)
    elseif(BLA_VENDOR MATCHES "")
        message(FATAL_ERROR "BLA_VENDOR Unset. Please Specify a library")
    endif()
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()

find_package(Threads)

find_package(OpenMP)
if(BLA_STATIC)
    set(MATIO_USE_STATIC_LIBRARIES TRUE) # XXX This should be an option
endif()

find_package(matio REQUIRED)

################
# VTK stuff
###############

if (USE_VTK)
    # what components do we want:
    set(VTK_FIND_COMPONENTS vtkIOXML vtkIOLegacy)
    mark_as_advanced(VTK_FIND_COMPONENTS)

    find_package(VTK REQUIRED COMPONENTS ${VTK_FIND_COMPONENTS})
    if (VTK_FOUND)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            add_compile_options(-Wno-inconsistent-missing-override)
        endif()
        # set(CMAKE_MSVCIDE_RUN_PATH ${VTK_RUNTIME_LIBRARY_DIRS} ${CMAKE_MSVCIDE_RUN_PATH}) # specially for windows
    endif()
endif()

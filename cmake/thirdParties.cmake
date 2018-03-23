###############
# FindMKL stuff
###############

set(BLA_DEFINITIONS)


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

################
# CGAL stuff
###############

if (USE_CGAL)
    # what components do we want:
    set(CGAL_FIND_COMPONENTS Core)
    mark_as_advanced(CGAL_FIND_COMPONENTS)

    # find_package(CGAL REQUIRED COMPONENTS ImageIO)
    # set(CGAL_LIBRARIES CGAL::CGAL_ImageIO)
    # add_definitions(-DCGAL_ImageIO)

    # # find_package(CGAL REQUIRED COMPONENTS Core OPTIONAL_COMPONENTS ImageIO) <- cannot since CGAL do not support OPTIONAL_COMPONENTS
    find_package(CGAL REQUIRED COMPONENTS ${CGAL_FIND_COMPONENTS})
    if (CGAL_FOUND)
        set(CGAL_LIBRARIES CGAL::CGAL_Core)
    else()
        message(FATAL_ERROR "Please set CGAL_DIR")
    endif()

    # find_package(CGAL QUIET COMPONENTS ImageIO)
    # # find_package(CGAL REQUIRED COMPONENTS Core OPTIONAL_COMPONENTS ImageIO) <- cannot since CGAL do not support OPTIONAL_COMPONENTS

    # if (CGAL_ImageIO_FOUND) # for handling images (.inr format only for the moment!)
    #     set(CGAL_LIBRARIES CGAL::CGAL_ImageIO)
    #     add_definitions(-DCGAL_ImageIO)
    # else()
    #     find_package(CGAL REQUIRED COMPONENTS ${CGAL_FIND_COMPONENTS})
    #     if (CGAL_FOUND)
    #         set(CGAL_LIBRARIES CGAL::CGAL_Core)
    #     else()
    #         message(FATAL_ERROR "Please set CGAL_DIR")
    #     endif()
    # endif()
    # set(CGAL_CXX_FLAGS ${CGAL_CXX_FLAGS_INIT} ${CGAL_SHARED_LINKER_FLAGS_INIT} ${CGAL_CXX_FLAGS_RELEASE_INIT} )
    # separate_arguments(CGAL_CXX_FLAGS) # needed to remove quotes/spaces problems
    # list(APPEND OpenMEEG_OTHER_LIBRARY_DIRS ${CGAL_LIBRARY_DIRS})
    # list(APPEND OpenMEEG_OTHER_INCLUDE_DIRS ${CGAL_INCLUDE_DIRS})
    # list(APPEND OpenMEEG_DEPENDENCIES CGAL)
    # if (CGAL_3RD_PARTY_LIBRARIES)
    #     # old CGAL (trusty 4.2.5.ubuntu)
    #     set(CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_Core_LIBRARY} ${CGAL_ImageIO_LIBRARY} ${MPFR_LIBRARIES} ${GMP_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_ImageIO_3RD_PARTY_LIBRARIES})
    #     set(CGAL_CXX_FLAGS ${CGAL_CXX_FLAGS} ${CGAL_ImageIO_3RD_PARTY_DEFINITIONS})
    # endif()
endif()

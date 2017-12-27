set(Atlas    USE_ATLAS)
set(OpenBLAS USE_OPENBLAS)
set(MKL      USE_MKL)
set(LAPACK   USE_LAPACK)
set(vecLib   USE_VECLIB)
set(Auto     USE_AUTO)

# the default case is Auto unless user tells something else
set(BLASLAPACK_IMPLEMENTATION_DEFAULT Auto)
if (BLASLAPACK_IMPLEMENTATION)
    set(BLASLAPACK_IMPLEMENTATION_DEFAULT ${BLASLAPACK_IMPLEMENTATION})
endif()

# the list of possibilites depending on the OS
set(LIST_IMPL "Auto" "MKL" "OpenBLAS" "LAPACK")
if (APPLE)
    set(LIST_IMPL ${LIST_IMPL} "vecLib")
elseif(UNIX)
    set(LIST_IMPL ${LIST_IMPL} "Atlas")
endif()

set(BLASLAPACK_IMPLEMENTATION_DOCSTRING "Choose the proper Blas/Lapack implementation: ${LIST_IMPL}")
set(BLASLAPACK_IMPLEMENTATION "${BLASLAPACK_IMPLEMENTATION_DEFAULT}" CACHE STRING "${BLASLAPACK_IMPLEMENTATION_DOCSTRING}" FORCE)

# Set the possible values of build type for cmake-gui
set_property(CACHE BLASLAPACK_IMPLEMENTATION PROPERTY STRINGS ${LIST_IMPL})

# Ensure that only one lapack implementation is selected by clearing all variable before setting the one chosen.
foreach (i Auto Atlas OpenBLAS MKL LAPACK vecLib)
    unset(${${i}})
endforeach()
set(${${BLASLAPACK_IMPLEMENTATION}} ON)

# unset unused variables
foreach (i Auto Atlas OpenBLAS MKL LAPACK vecLib)
    unset(${i})
endforeach()

### now actually include the files

# include the wanted BlasLapack
set(FIND_MODE "REQUIRED")
include(Use${BLASLAPACK_IMPLEMENTATION})

if (${CMAKE_PROJECT_NAME} STREQUAL "OpenMEEG" OR LAPACK_LIBRARIES)

    message(STATUS "Lapack package: ${LAPACK}: ${LAPACK_LIBRARIES}")

    # Detect Fortran to C interface.
    if (NOT USE_MKL AND NOT WIN32 AND NOT USE_OPENBLAS)
        include(FortranCInterface)
        FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_" FROM_LIBRARY blas[daxpy] HINTS ${lapack_DIR}/lib)
    endif()

    # TODO windows goes in this ....
    # Last chance check...
    if (NOT LAPACK_LIBRARIES)
        set(LAPACK_DIR ${LAPACK_DIR_SAVE})
        set(LAPACK_libs_dir ${LAPACK_DIR}/${INSTALL_LIB_DIR})
        message("Searching LAPACK in ${LAPACK_libs_dir}")
        set(CMAKE_FIND_DEBUG_MODE 1)
        find_library(lapack NAMES lapack lapackd HINTS ${LAPACK_libs_dir})
        find_library(blas   NAMES blas   blasd   HINTS ${LAPACK_libs_dir})
        find_library(f2c    NAMES f2c f2cd libf2c libf2cd HINTS ${LAPACK_libs_dir})
        if (NOT (lapack AND blas AND f2c))
            message(SEND_ERROR "clapack is needed")
        endif()
        message("Lapack package: ${LAPACK}: ${LAPACK_LIBRARIES}")
        get_filename_component(LAPACK_DLL_DIR "${lapack}" DIRECTORY)
        set(LAPACK_LIBRARIES ${lapack} ${blas} ${f2c})
        if (NOT BUILD_SHARED_LIBS)
            file(GLOB GCC_fileS "/usr/lib/gcc/*/*")
            find_file(GFORTRAN_LIB libgfortran.a ${GCC_fileS})
            set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${GFORTRAN_LIB})
        endif()
    endif()

    set(HAVE_LAPACK TRUE)
    set(HAVE_BLAS TRUE)
else()
    # if not we can build our LAPACK
    # if still no LAPACK and we are in the superproject then build it
    unset(USE_AUTO CACHE)
    set(USE_LAPACK True CACHE BOOL "Force Lapack build" FORCE)
    set(BLASLAPACK_IMPLEMENTATION "LAPACK" CACHE STRING "Forced to LAPACK as no other BLASLAPACK_IMPLEMENTATION found." FORCE)
    message("Forcing BLASLAPACK_IMPLEMENTATION to LAPACK to build it as no other BLASLAPACK_IMPLEMENTATION found.")
endif()

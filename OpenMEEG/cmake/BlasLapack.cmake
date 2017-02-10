
if (NOT BUILD_SHARED_LIBS)
    if (WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib;.dll")
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
        if (APPLE)
            set(CMAKE_FIND_LIBRARY_SUFFIXES "${CMAKE_FIND_LIBRARY_SUFFIXES};.dylib")
        endif()
    endif()
endif()

# include the wanted BlasLapack
set(REQUIRED "REQUIRED")
include(Use${BLASLAPACK_IMPLEMENTATION})
message("Lapack package: ${LAPACK}: ${LAPACK_LIBRARIES}")

#   Detect Fortran to C interface.
if (NOT USE_MKL AND NOT WIN32)
    set(FC_INTERFACE)
    include(FortranCInterface)
    FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_" FROM_LIBRARY blas[daxpy])
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

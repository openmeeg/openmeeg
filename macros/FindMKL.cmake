# - Try to find the Intel Math Kernel Library
# Once done this will define
#
#  MKL_FOUND - system has MKL
#  MKL_ROOT_DIR - path to the MKL base directory
#  MKL_INCLUDE_DIR - the MKL include directory
#  MKL_LIBRARIES - MKL libraries
#
# There are few sets of libraries:
# Array indexes modes:
# LP - 32 bit indexes of arrays
# ILP - 64 bit indexes of arrays
# Threading:
# SEQUENTIAL - no threading
# INTEL - Intel threading library
# GNU - GNU threading library
# MPI support
# NOMPI - no MPI support
# INTEL - Intel MPI library
# OPEN - Open MPI library
# SGI - SGI MPT Library

# linux
if(UNIX AND NOT APPLE)
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        set(MKL_ARCH_DIR "em64t")
    else()
        set(MKL_ARCH_DIR "32")
    endif()
endif()

# macos
if(APPLE)
    set(MKL_ARCH_DIR "em64t")
endif()

IF(FORCE_BUILD_32BITS)
    set(MKL_ARCH_DIR "32")
ENDIF()

set (MKL_THREAD_VARIANTS SEQUENTIAL GNUTHREAD INTELTHREAD)
set (MKL_MODE_VARIANTS ILP LP)
set (MKL_MPI_VARIANTS NOMPI INTELMPI OPENMPI SGIMPT)


find_path(MKL_ROOT_DIR
    include/mkl_cblas.h
    PATHS
    $ENV{MKLDIR}
    /opt/intel/mkl/*/
    /opt/intel/cmkl/*/
    /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal
)

# MESSAGE("MKL_ROOT_DIR : ${MKL_ROOT_DIR}") # for debug

find_path(MKL_INCLUDE_DIR
  mkl_cblas.h
  PATHS
    ${MKL_ROOT_DIR}/include
    ${INCLUDE_INSTALL_DIR}
)

find_path(MKL_FFTW_INCLUDE_DIR
  fftw3.h
  PATH_SUFFIXES fftw
  PATHS
    ${MKL_ROOT_DIR}/include
    ${INCLUDE_INSTALL_DIR}
  NO_DEFAULT_PATH
)

find_library(MKL_CORE_LIBRARY
  mkl_core
  PATHS
    ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
)

# Threading libraries
find_library(MKL_SEQUENTIAL_LIBRARY
  mkl_sequential
  PATHS
    ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
)

find_library(MKL_INTELTHREAD_LIBRARY
  mkl_intel_thread
  PATHS
    ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
)

find_library(MKL_GNUTHREAD_LIBRARY
  mkl_gnu_thread
  PATHS
    ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
)

# Intel Libraries
IF("${MKL_ARCH_DIR}" STREQUAL "32")
    find_library(MKL_LP_LIBRARY
      mkl_intel
      PATHS
        ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
    )

    find_library(MKL_ILP_LIBRARY
      mkl_intel
      PATHS
        ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
    )
else()
    find_library(MKL_LP_LIBRARY
      mkl_intel_lp64
      PATHS
        ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
    )

    find_library(MKL_ILP_LIBRARY
      mkl_intel_ilp64
      PATHS
        ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
    )
ENDIF()

# Lapack
find_library(MKL_LAPACK_LIBRARY
  mkl_lapack
  PATHS
    ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}
)

# iomp5
find_library(MKL_IOMP5_LIBRARY
  iomp5
  PATHS
    ${MKL_ROOT_DIR}/../lib/intel64
)

foreach (MODEVAR ${MKL_MODE_VARIANTS})
    foreach (THREADVAR ${MKL_THREAD_VARIANTS})
        if (MKL_CORE_LIBRARY AND MKL_${MODEVAR}_LIBRARY AND MKL_${THREADVAR}_LIBRARY)
            set(MKL_${MODEVAR}_${THREADVAR}_LIBRARIES
                ${MKL_${MODEVAR}_LIBRARY} ${MKL_${THREADVAR}_LIBRARY} ${MKL_CORE_LIBRARY}
                ${MKL_LAPACK_LIBRARY} ${MKL_IOMP5_LIBRARY})
            # message("${MODEVAR} ${THREADVAR} ${MKL_${MODEVAR}_${THREADVAR}_LIBRARIES}") # for debug
        endif()
    endforeach()
endforeach()

set(MKL_LIBRARIES ${MKL_LP_SEQUENTIAL_LIBRARIES})
LINK_DIRECTORIES(${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}) # hack

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARIES)

mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARIES
    MKL_CORE_LIBRARY MKL_LP_LIBRARY MKL_ILP_LIBRARY
    MKL_SEQUENTIAL_LIBRARY MKL_INTELTHREAD_LIBRARY MKL_GNUTHREAD_LIBRARY
)

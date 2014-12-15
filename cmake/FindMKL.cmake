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

set(MKL_ARCH_DIR "ia32")
if (${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    set(MKL_ARCH_DIR "intel64")
endif()

if (FORCE_BUILD_32BITS)
    set(MKL_ARCH_DIR "ia32")
endif()

set(MKL_THREAD_VARIANTS SEQUENTIAL GNUTHREAD INTELTHREAD)
set(MKL_MODE_VARIANTS ILP LP)
set(MKL_MPI_VARIANTS NOMPI INTELMPI OPENMPI SGIMPT)

set(CMAKE_FIND_DEBUG_MODE 1)

set(MKL_POSSIBLE_LOCATIONS
    $ENV{MKLDIR}
    /opt/intel/mkl
    /opt/intel/cmkl
    /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal
    "C:/Program Files (x86)/Intel/ComposerXE-2011/mkl"
    "C:/Program Files (x86)/Intel/Composer XE 2013/mkl"
    "C:/Program Files/Intel/MKL/*/"
    "C:/Program Files/Intel/ComposerXE-2011/mkl"
    "C:/Program Files/Intel/Composer XE 2013/mkl"
    "C:/Program Files (x86)/Intel/Composer XE 2015/mkl/"
    "C:/Program Files/Intel/Composer XE 2015/mkl/"
)

foreach (i ${MKL_POSSIBLE_LOCATIONS})
    if (EXISTS ${i}/include/mkl_cblas.h)
        set(MKL_ROOT_DIR ${i})
        break()
    endif()
endforeach()

#   Does This work at all ?
find_path(MKL_ROOT_DIR NAMES include/mkl_cblas.h PATHS ${MKL_POSSIBLE_LOCATIONS})

IF (NOT MKL_ROOT_DIR)
	MESSAGE(WARNING "Could not find MKL: disabling it")
	set(USE_MKL FALSE)
endif()

if (USE_MKL)
    find_path(MKL_INCLUDE_DIR mkl_cblas.h PATHS ${MKL_ROOT_DIR}/include ${INCLUDE_INSTALL_DIR})

    find_path(MKL_FFTW_INCLUDE_DIR fftw3.h PATH_SUFFIXES fftw PATHS ${MKL_ROOT_DIR}/include ${INCLUDE_INSTALL_DIR} NO_DEFAULT_PATH)

    if (WIN32)
        set(MKL_LIB_SEARCHPATH $ENV{ICC_LIB_DIR} $ENV{MKL_LIB_DIR} "${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}" "${MKL_ROOT_DIR}/../compiler" "${MKL_ROOT_DIR}/../compiler/lib/${MKL_ARCH_DIR}")
        
        if (MKL_INCLUDE_DIR MATCHES "10.")
            set(MKL_LIBS mkl_solver mkl_core mkl_intel_c mkl_intel_s mkl_intel_thread libguide mkl_lapack95 mkl_blas95)
            if (CMAKE_CL_64)
                set(MKL_LIBS mkl_solver_lp64 mkl_core mkl_intel_lp64 mkl_intel_thread libguide mkl_lapack95_lp64 mkl_blas95_lp64)
            endif()
        elseif(MKL_INCLUDE_DIR MATCHES "2013") # version 11 ...
            set(MKL_LIBS mkl_core mkl_intel_c mkl_intel_s mkl_intel_thread libiomp5md mkl_lapack95 mkl_blas95)
            if (CMAKE_CL_64)
                set(MKL_LIBS mkl_core mkl_intel_lp64 mkl_intel_thread libiomp5md mkl_lapack95_lp64 mkl_blas95_lp64)
            endif()
        elseif(MKL_INCLUDE_DIR MATCHES "2015")
            if(CMAKE_CL_64)
                SET(MKL_LIBS mkl_intel_lp64 mkl_core mkl_intel_thread mkl_lapack95_lp64 mkl_blas95_lp64 )
            else()
                SET(MKL_LIBS mkl_intel_c mkl_core mkl_intel_thread mkl_lapack95 mkl_blas95 )
            endif()
        else() # old MKL 9
            set(MKL_LIBS mkl_solver mkl_c libguide mkl_lapack mkl_ia32)
        endif()

        if (MKL_INCLUDE_DIR MATCHES "10.3")
            set(MKL_LIBS ${MKL_LIBS} libiomp5md)
        endif()
        
        foreach (LIB ${MKL_LIBS})
            find_library(${LIB}_PATH ${LIB} PATHS ${MKL_LIB_SEARCHPATH} ENV LIBRARY_PATH)
            if (${LIB}_PATH)
                set(MKL_LIBRARIES ${MKL_LIBRARIES} ${${LIB}_PATH})
            else()
                message(FATAL_ERROR "Could not find ${LIB}: disabling MKL")
                BREAK()
            endif()
        endforeach()
        set(MKL_FOUND ON)

    else() # UNIX and macOS

        set(MKL_LIBRARY_LOCATIONS ${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR} ${MKL_ROOT_DIR}/lib)

        find_library(MKL_CORE_LIBRARY mkl_core PATHS ${MKL_LIBRARY_LOCATIONS})

        # Threading libraries

        find_library(MKL_RT_LIBRARY mkl_rt PATHS ${MKL_LIBRARY_LOCATIONS})
        find_library(MKL_SEQUENTIAL_LIBRARY mkl_sequential PATHS ${MKL_LIBRARY_LOCATIONS})
        find_library(MKL_INTELTHREAD_LIBRARY mkl_intel_thread PATHS ${MKL_LIBRARY_LOCATIONS})
        find_library(MKL_GNUTHREAD_LIBRARY mkl_gnu_thread PATHS ${MKL_LIBRARY_LOCATIONS})

        # Intel Libraries

        if (NOT "${MKL_ARCH_DIR}" STREQUAL "ia32")
            set(INTEL_LP_SUFFIX  "_lp64")
            set(INTEL_ILP_SUFFIX "_ilp64")
        endif()

        find_library(MKL_LP_LIBRARY mkl_intel%{INTEL_LP_SUFFIX} PATHS ${MKL_LIBRARY_LOCATIONS})
        find_library(MKL_ILP_LIBRARY mkl_intel${INTEL_ILP_SUFFIX} PATHS ${MKL_LIBRARY_LOCATIONS})

        # Lapack

        find_library(MKL_LAPACK_LIBRARY mkl_lapack PATHS ${MKL_LIBRARY_LOCATIONS})

        if (NOT MKL_LAPACK_LIBRARY)
            find_library(MKL_LAPACK_LIBRARY mkl_lapack95_lp64 PATHS ${MKL_LIBRARY_LOCATIONS})
        endif()

        # iomp5

        if (UNIX AND NOT APPLE)
            find_library(MKL_IOMP5_LIBRARY iomp5 PATHS ${MKL_ROOT_DIR}/../lib/${MKL_ARCH_DIR})
        endif()

        foreach (MODEVAR ${MKL_MODE_VARIANTS})
            foreach (THREADVAR ${MKL_THREAD_VARIANTS})
                if (MKL_CORE_LIBRARY AND MKL_${MODEVAR}_LIBRARY AND MKL_${THREADVAR}_LIBRARY)
                    set(MKL_${MODEVAR}_${THREADVAR}_LIBRARIES
                        ${MKL_${MODEVAR}_LIBRARY} ${MKL_${THREADVAR}_LIBRARY} ${MKL_CORE_LIBRARY}
                        ${MKL_LAPACK_LIBRARY} ${MKL_IOMP5_LIBRARY})
                    message("${MODEVAR} ${THREADVAR} ${MKL_${MODEVAR}_${THREADVAR}_LIBRARIES}") # for debug
                endif()
            endforeach()
        endforeach()

        set(MKL_LIBRARIES ${MKL_RT_LIBRARY})
        mark_as_advanced(MKL_CORE_LIBRARY MKL_LP_LIBRARY MKL_ILP_LIBRARY
            MKL_SEQUENTIAL_LIBRARY MKL_INTELTHREAD_LIBRARY MKL_GNUTHREAD_LIBRARY)
    endif()
            
    #link_directories(${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}) # hack

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARIES)

    mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARIES)
endif()

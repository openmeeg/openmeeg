IF (WIN32)       # WINDOWS
#        OPTION(USE_ATLAS "Build the project using ATLAS" OFF)
        OPTION(USE_MKL "Build the project with MKL" ON)
        MARK_AS_ADVANCED(USE_MKL)
#        MARK_AS_ADVANCED(USE_ATLAS)
        IF (NOT BUILD_SHARED)
            set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib;.dll")
        ENDIF()
ELSE()
    OPTION(USE_ATLAS "Build the project using ATLAS" ON)
    OPTION(USE_MKL "Build the project with MKL" OFF)
    IF (APPLE)   # MACOSX
        IF (NOT BUILD_SHARED)
            set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so;.dylib")
        ENDIF()
    ELSE() # LINUX
        IF (NOT BUILD_SHARED)
            set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
        ENDIF()
    ENDIF()
ENDIF()

IF (USE_MKL)
    IF (WIN32)
        FILE(GLOB MKL_PATH "C:/Program Files/Intel/MKL/*")
        SET(MKL_LIB_SEARCHPATH "${MKL_PATH}/ia32/lib")
        FIND_PATH(MKL_INCLUDE_PATH mkl.h "${MKL_PATH}/include")
        IF (MKL_INCLUDE_PATH MATCHES "10.")
            SET(MKL_LIBS mkl_solver mkl_core mkl_intel_c mkl_intel_s mkl_intel_thread libguide mkl_lapack95 mkl_blas95)
        ELSE()
            SET(MKL_LIBS mkl_solver mkl_c libguide mkl_lapack mkl_ia32)
        ENDIF()
        INCLUDE_DIRECTORIES(${MKL_INCLUDE_PATH})

        FOREACH (LIB ${MKL_LIBS})
            FIND_LIBRARY(${LIB}_PATH ${LIB} PATHS ${MKL_LIB_SEARCHPATH} ENV LIBRARY_PATH)

            IF(${LIB}_PATH)
                SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${${LIB}_PATH})
                #MESSAGE("${LIB} found in ${${LIB}_PATH}")
                MARK_AS_ADVANCED(${LIB}_PATH)
            ELSE()
                MESSAGE("Could not find ${LIB}")
            ENDIF()
        ENDFOREACH()
    ELSE()
        FIND_PACKAGE(MKL)
        IF (MKL_FOUND)
            INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
            SET(LAPACK_LIBRARIES ${MKL_LIBRARIES})
            # MESSAGE(${LAPACK_LIBRARIES}) # for debug

            IF(UNIX AND NOT APPLE) # MKL on linux requires to link with the pthread library
                SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} pthread)
            ENDIF()
        ELSE()
            MESSAGE(FATAL_ERROR "MKL not found. Please set environment variable MKLDIR")
        ENDIF()
    ENDIF()
ELSE()

    #   ATLAS OR LAPACK/BLAS
    IF (UNIX AND NOT APPLE)
        IF (USE_ATLAS)
            SET(ATLAS_LIB_SEARCHPATH
                /usr/lib64/
                /usr/lib64/atlas
                /usr/lib64/atlas/sse2
                /usr/lib/atlas/sse2
                /usr/lib/sse2
                /usr/lib64/atlas/sse3
                /usr/lib/atlas/sse3
                /usr/lib/sse3
                /usr/lib/
                /usr/lib/atlas
                )
            SET(ATLAS_LIBS atlas cblas f77blas clapack lapack blas)

            FIND_PATH(ATLAS_INCLUDE_PATH cblas.h /usr/include/ /usr/include/atlas)
            MARK_AS_ADVANCED(ATLAS_INCLUDE_PATH)
            INCLUDE_DIRECTORIES(${ATLAS_INCLUDE_PATH})
            FOREACH (LIB ${ATLAS_LIBS})
                SET(LIBNAMES ${LIB})
                IF (${LIB} STREQUAL "clapack")
                    SET(LIBNAMES ${LIB} lapack_atlas)
                ENDIF()
                FIND_LIBRARY(${LIB}_PATH
                    NAMES ${LIBNAMES}
                    PATHS ${ATLAS_LIB_SEARCHPATH}
                    NO_DEFAULT_PATH
                    NO_CMAKE_ENVIRONMENT_PATH
                    NO_CMAKE_PATH
                    NO_SYSTEM_ENVIRONMENT_PATH
                    NO_CMAKE_SYSTEM_PATH)
                IF(${LIB}_PATH)
                    SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${${LIB}_PATH})
                    MARK_AS_ADVANCED(${LIB}_PATH)
                ELSE()
                    MESSAGE(WARNING "Could not find ${LIB}")
                ENDIF()
            ENDFOREACH()
        ELSE()
            IF (lapack_PATH AND blas_PATH)
                SET(LAPACKBLAS_LIB_SEARCHPATH
                    /usr/lib64/
                    /usr/lib/
                    )
                #FOREACH (LIB cblas f77blas clapack lapack blas)
                FOREACH (LIB lapack blas)
                    FIND_LIBRARY(${LIB}_PATH
                        NAMES ${LIB}
                        PATHS ${LAPACKBLAS_LIB_SEARCHPATH}
                        NO_DEFAULT_PATH
                        NO_CMAKE_ENVIRONMENT_PATH
                        NO_CMAKE_PATH
                        NO_SYSTEM_ENVIRONMENT_PATH
                        NO_CMAKE_SYSTEM_PATH)
                    IF(${LIB}_PATH)
                        SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${${LIB}_PATH})
                        MARK_AS_ADVANCED(${LIB}_PATH)
                    ELSE()
                        MESSAGE(WARNING "Could not find ${LIB}")
                    ENDIF()
                ENDFOREACH()
            ENDIF()
        ENDIF()
        SET(HAVE_LAPACK TRUE)
        SET(HAVE_BLAS TRUE)

        IF (NOT BUILD_SHARED)
            FILE(GLOB GCC_FILES "/usr/lib/gcc/*/*")
            FIND_FILE(GFORTRAN_LIB libgfortran.a ${GCC_FILES})
            SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${GFORTRAN_LIB})
        ENDIF()
    ELSE() # Assume APPLE
        SET(LAPACK_LIBRARIES "-framework vecLib")
        SET(HAVE_LAPACK TRUE)
        SET(HAVE_BLAS TRUE)
        INCLUDE_DIRECTORIES(/System/Library/Frameworks/vecLib.framework/Headers)
    ENDIF()
ENDIF()

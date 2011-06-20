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
        FIND_PATH(MKL_ROOT_DIR
            include/mkl.h
            PATHS
            $ENV{MKL_DIR}
            "C:/Program Files/Intel/MKL/*/"
            "C:/Program Files/Intel/ComposerXE-2011/mkl"
            "C:/Program Files (x86)/Intel/ComposerXE-2011/mkl"
        )

        IF (MKL_ROOT_DIR)
            if (WIN32)
                if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
                    set(MKL_ARCH_DIR "intel64")
                else()
                    set(MKL_ARCH_DIR "ia32")
                endif()
            endif()

            SET(MKL_LIB_SEARCHPATH $ENV{ICC_LIB_DIR} $ENV{MKL_LIB_DIR} "${MKL_ROOT_DIR}/lib/${MKL_ARCH_DIR}" "${MKL_ROOT_DIR}/../compiler")

            FIND_PATH(MKL_INCLUDE_PATH mkl.h ${MKL_ROOT_DIR}/include)

            IF (MKL_INCLUDE_PATH MATCHES "10." OR MKL_INCLUDE_PATH MATCHES "11.")
                IF(CMAKE_CL_64)
                    SET(MKL_LIBS mkl_solver_lp64 mkl_core mkl_intel_lp64 mkl_intel_thread mkl_lapack95_lp64 mkl_blas95_lp64)
                ELSE()
                    SET(MKL_LIBS mkl_solver mkl_core mkl_intel_c mkl_intel_s mkl_intel_thread libguide mkl_lapack95 mkl_blas95)
                ENDIF()
            ELSE() # old MKL 9
                SET(MKL_LIBS mkl_solver mkl_c libguide mkl_lapack mkl_ia32)
            ENDIF()

            IF (MKL_INCLUDE_PATH MATCHES "10.3")
                SET(MKL_LIBS ${MKL_LIBS} libiomp5md)
            ENDIF()

            INCLUDE_DIRECTORIES(${MKL_INCLUDE_PATH})

            FOREACH (LIB ${MKL_LIBS})
                FIND_LIBRARY(${LIB}_PATH ${LIB} PATHS ${MKL_LIB_SEARCHPATH} ENV LIBRARY_PATH)

                IF(${LIB}_PATH)
                    SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${${LIB}_PATH})
                    MARK_AS_ADVANCED(${LIB}_PATH)
                ELSE()
                    MESSAGE(FATAL_ERROR "Could not find ${LIB}: disabling MKL")
                    BREAK()
                ENDIF()
            ENDFOREACH()
        ELSE()
            MESSAGE(WARNING "Could not find MKL: disabling it")
            SET(USE_MKL FALSE)
        ENDIF()
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
ENDIF()

IF (NOT USE_MKL)
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
    ELSE() # Assume APPLE or local lapack/blas (treated in contrib)
        IF (APPLE)
            SET(LAPACK_LIBRARIES "-framework vecLib")
            INCLUDE_DIRECTORIES(/System/Library/Frameworks/vecLib.framework/Headers)
            SET(HAVE_LAPACK TRUE)
            SET(HAVE_BLAS TRUE)
        ENDIF()
    ENDIF()
ENDIF()

IF (NOT HAVE_LAPACK)
    SET(HAVE_LAPACK TRUE)
    SET(HAVE_BLAS TRUE)
    SET(NEED_CLAPACK TRUE)
    SET(LAPACK_LIBRARIES lapack blas f2c)
ENDIF()

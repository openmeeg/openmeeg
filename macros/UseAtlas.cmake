IF (WIN32)       # WINDOWS
        OPTION(USE_ATLAS "Build the project using ATLAS" OFF)
        OPTION(USE_MKL "Build the project with MKL" ON)
        MARK_AS_ADVANCED(USE_MKL)
        MARK_AS_ADVANCED(USE_ATLAS)
        IF (NOT BUILD_SHARED)
            set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib;.dll")
        ENDIF()
ELSE()
    IF (APPLE)   # MACOSX
        OPTION(USE_ATLAS "Build the project using ATLAS" OFF)
        OPTION(USE_MKL "Build the project with MKL" ON)
        IF (NOT BUILD_SHARED)
            set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so;.dylib")
        ENDIF()
    ELSE() # LINUX
        OPTION(USE_ATLAS "Build the project using ATLAS" ON)
        OPTION(USE_MKL "Build the project with MKL" OFF)
        IF (NOT BUILD_SHARED)
            set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
        ENDIF()
    ENDIF()
ENDIF()

IF(USE_ATLAS)
    IF (WIN32)
        MESSAGE(FATAL_ERROR "Atlas not supported on Windows. Please use MKL")
    ENDIF()

    IF (UNIX AND NOT APPLE)
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

        IF (NOT BUILD_SHARED)
            FIND_FILE(GFORTRAN_LIB libgfortran.a
                /usr/lib/gcc/i486-linux-gnu/4.2.4/
                /usr/lib/gcc/i386-redhat-linux/4.3.2/
                /usr/lib/gcc/x86_64-manbo-linux-gnu/4.3.2/
                /usr/lib/gcc/x86_64-redhat-linux/4.4.3/
            )
    #	    SET(GFORTRAN_LIB
    #    /usr/lib/gcc/i486-linux-gnu/4.2.4/libgfortran.a)
    #    /usr/lib/gcc/i386-redhat-linux/4.3.2/libgfortran.a)
    #    /usr/lib/gcc/x86_64-manbo-linux-gnu/4.3.2/libgfortran.a)
            SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${GFORTRAN_LIB})
        ENDIF()
# 
        FIND_PATH(ATLAS_INCLUDE_PATH cblas.h /usr/include/ /usr/include/atlas)
        MARK_AS_ADVANCED(ATLAS_INCLUDE_PATH)
        INCLUDE_DIRECTORIES(${ATLAS_INCLUDE_PATH})
    ELSE() # Assume APPLE
        SET(LAPACK_LIBRARIES "-framework vecLib")
        INCLUDE_DIRECTORIES(/System/Library/Frameworks/vecLib.framework/Headers)
    ENDIF()
ENDIF()

IF (USE_MKL)

    IF ( WIN32 )
        FILE(GLOB MKL_PATH "C:/Program Files/Intel/MKL/*")
    ENDIF( WIN32 )
    FIND_PATH(MKL_INCLUDE_PATH mkl.h "${MKL_PATH}/include")

    IF ( MKL_INCLUDE_PATH )
        #MESSAGE("mkl.h found in ${MKL_INCLUDE_PATH}")
        INCLUDE_DIRECTORIES(${MKL_INCLUDE_PATH})
    ELSE()
        MESSAGE("Can not find mkl.h")
    ENDIF()

    IF (UNIX AND NOT APPLE)
        SET(MKL_LIB_SEARCHPATH # add here some paths to look for mkl libs
            ""
        )
    ENDIF()

    IF (APPLE)
        SET(MKL_LIB_SEARCHPATH # add here some paths to look for mkl libs
            /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal
            /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/32
            /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/64
            #/Library/Frameworks/Intel_MKL.framework/Libraries/32
            #/Library/Frameworks/Intel_MKL.framework/Libraries/universal
            /opt/intel/Compiler/11.0/056/lib
        )
    ENDIF()

    IF (WIN32)
        SET(MKL_LIB_SEARCHPATH "${MKL_PATH}/ia32/lib")
    ENDIF()

    IF (WIN32)
        IF (MKL_INCLUDE_PATH MATCHES "10.")
            SET(MKL_LIBS mkl_solver mkl_core mkl_intel_c mkl_intel_s mkl_intel_thread libguide mkl_lapack95 mkl_blas95)
        ELSE (MKL_INCLUDE_PATH MATCHES "10.")
            SET(MKL_LIBS mkl_solver mkl_c libguide mkl_lapack mkl_ia32)
        ENDIF (MKL_INCLUDE_PATH MATCHES "10.")
    ELSE()
        SET(MKL_LIBS mkl_intel mkl_intel_thread mkl_core iomp5 pthread)
        #SET(MKL_LIBS mkl_intel mkl_intel_thread mkl_core iomp5md pthread)
        #SET(MKL_LIBS mkl_intel mkl_core mkl_lapack)
        #SET(MKL_LIBS mkl_intel_lp64 mkl_core mkl_lapack)
        #SET(MKL_LIBS mkl guide mkl_lapack) % for old MKL
    ENDIF()

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

    IF(UNIX AND NOT APPLE) # MKL on linux requires to link with the pthread library
        SET(LAPACK_LIBRARIES "${LAPACK_LIBRARIES} pthread")
    ENDIF()
ENDIF()

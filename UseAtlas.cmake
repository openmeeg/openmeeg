IF (WIN32)       # WINDOWS
        OPTION(USE_ATLAS "Build the project using ATLAS" OFF)
        OPTION(USE_MKL "Build the project with MKL" ON)
ELSE (WIN32)
    IF (APPLE)   # MACOSX
        OPTION(USE_ATLAS "Build the project using ATLAS" OFF)
        OPTION(USE_MKL "Build the project with MKL" ON)
    ELSE (APPLE) # LINUX
        OPTION(USE_ATLAS "Build the project using ATLAS" ON)
        OPTION(USE_MKL "Build the project with MKL" OFF)
    ENDIF (APPLE)
ENDIF (WIN32)

IF(USE_ATLAS)

    SET(ATLAS_LIB_SEARCHPATH
        /usr/lib64/atlas
        /usr/lib/sse2
        /usr/lib/
    )

    SET(ATLAS_OTHER_LIBS lapack cblas)

    # Find lib atlas and assume ${ATLAS_OTHER_LIBS} are in the same directory
    FIND_LIBRARY(ATLAS_LIB
                        NAMES atlas
                        PATHS ${ATLAS_LIB_SEARCHPATH}
                        NO_DEFAULT_PATH
                        NO_CMAKE_ENVIRONMENT_PATH
                        NO_CMAKE_PATH
                        NO_SYSTEM_ENVIRONMENT_PATH
                        NO_CMAKE_SYSTEM_PATH)

    SET(OPENMEEG_OTHER_LIBRARIES
        ${OPENMEEG_OTHER_LIBRARIES} ${ATLAS_LIB} ${ATLAS_OTHER_LIBS})
    #MARK_AS_ADVANCED(${ATLAS_LIB})

    FIND_PATH(ATLAS_INCLUDE_PATH atlas/cblas.h
                /usr/include/
    )
    INCLUDE_DIRECTORIES(${ATLAS_INCLUDE_PATH})

ENDIF(USE_ATLAS)

IF ( USE_MKL )

    FIND_PATH(MKL_INCLUDE_PATH mkl.h
                "C:/Program Files/Intel/MKL/8.1.1/include"
                ../../mkl/include
                ~/intel/mkl/8.1/include
                ~/Intel/MKL/8.1/include
                /Library/Frameworks/Intel_MKL.framework/Headers

    )
    INCLUDE_DIRECTORIES(${MKL_INCLUDE_PATH})

    SET(MKL_LIB_SEARCHPATH
        "C:/Program Files/Intel/MKL/8.1.1/ia32/lib"
        ../../mkl/ia32/lib
        ~/intel/mkl/8.1/lib/32
        ~/Intel/MKL/8.1/lib/32
        /Library/Frameworks/Intel_MKL.framework/Libraries/32
    )

    IF ( WIN32 )
        SET(MKL_LIBS mkl_solver mkl_c guide)
    ELSE ( WIN32 )
        SET(MKL_LIBS mkl mkl_lapack mkl_ia32 guide)
    ENDIF ( WIN32 )

    FOREACH ( LIB ${MKL_LIBS})
        FIND_LIBRARY(${LIB}_PATH ${LIB} ${MKL_LIB_SEARCHPATH})

        IF(${LIB}_PATH)
            SET(OPENMEEG_OTHER_LIBRARIES
                ${OPENMEEG_OTHER_LIBRARIES} ${${LIB}_PATH})
            #MESSAGE("${LIB} found in ${${LIB}_PATH}")
            MARK_AS_ADVANCED(${LIB}_PATH)
        ELSE(${LIB}_PATH)
            MESSAGE("Could not find ${LIB}")
        ENDIF(${LIB}_PATH)

    ENDFOREACH ( LIB )

ENDIF ( USE_MKL )

IF (WIN32)       # WINDOWS
        OPTION(USE_ATLAS "Build the project using ATLAS" OFF)
        OPTION(USE_MKL "Build the project with MKL" ON)
        MARK_AS_ADVANCED(USE_MKL)
        MARK_AS_ADVANCED(USE_ATLAS)
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
        /usr/lib/atlas
    )

    SET(ATLAS_OTHER_LIBS lapack_atlas lapack cblas)

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
                "C:/Program Files/Intel/MKL/9.1.027/include"
                "C:/Program Files/Intel/MKL/8.1.1/include"
    )
    IF ( MKL_INCLUDE_PATH )
        #MESSAGE("mkl.h found in ${MKL_INCLUDE_PATH}")
        INCLUDE_DIRECTORIES(${MKL_INCLUDE_PATH})
    ELSE ( MKL_INCLUDE_PATH )
        MESSAGE("Can not find mkl.h")
    ENDIF ( MKL_INCLUDE_PATH )

    IF ( LINUX )
        SET(MKL_LIB_SEARCHPATH # add here some paths to look for mkl libs
            ""
        )
    ENDIF ( LINUX )

    IF ( APPLE )
        SET(MKL_LIB_SEARCHPATH # add here some paths to look for mkl libs
            /Library/Frameworks/Intel_MKL.framework/Libraries/universal
        )
    ENDIF ( APPLE )

    IF ( WIN32 )
        SET(MKL_LIB_SEARCHPATH
            "C:/Program Files/Intel/MKL/9.1.027/ia32/lib"
            "C:/Program Files/Intel/MKL/8.1.1/ia32/lib"
        )
    ENDIF ( WIN32 )

    IF ( WIN32 )
        SET(MKL_LIBS mkl_solver mkl_c libguide mkl_lapack mkl_ia32)
    ELSE ( WIN32 )
        SET(MKL_LIBS mkl guide mkl_lapack )
    ENDIF ( WIN32 )

    FOREACH ( LIB ${MKL_LIBS} )
        FIND_LIBRARY(${LIB}_PATH ${LIB}
            PATHS "${MKL_LIB_SEARCHPATH}"
            ENV LIBRARY_PATH
            )

        IF(${LIB}_PATH)
            SET(OPENMEEG_OTHER_LIBRARIES
                ${OPENMEEG_OTHER_LIBRARIES} ${${LIB}_PATH})
            #MESSAGE("${LIB} found in ${${LIB}_PATH}")
            MARK_AS_ADVANCED(${LIB}_PATH)
        ELSE(${LIB}_PATH)
            MESSAGE("Could not find ${LIB}")
        ENDIF(${LIB}_PATH)

    ENDFOREACH ( LIB )

    IF( LINUX ) # MKL on linux requires to link with the pthread library
        SET(OPENMEEG_OTHER_LIBRARIES "${OPENMEEG_OTHER_LIBRARIES} pthread")
    ENDIF( LINUX )

ENDIF ( USE_MKL )

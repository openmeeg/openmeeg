#######################################################################
# Setting compilation options
#######################################################################

IF(${CMAKE_C_COMPILER} MATCHES "gcc")
    SET(USE_GCC YES)
ENDIF()

IF(${CMAKE_C_COMPILER} MATCHES "icc")
    SET(USE_ICC YES)
ENDIF()

IF(UNIX)
    OPTION(FORCE_BUILD_32BITS "Force 32 bits compilation" OFF)
    IF(FORCE_BUILD_32BITS)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
    ENDIF()

    MARK_AS_ADVANCED(FORCE_BUILD_32BITS)

    IF(NOT BUILD_SHARED_LIBS)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    ENDIF()
ENDIF()


IF (USE_GCC)
    SET(CXX_WARNING_OPTIONS "-Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings")
    SET(CC_WARNING_OPTIONS "-Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wwrite-strings")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CXX_WARNING_OPTIONS}")
    SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${CC_WARNING_OPTIONS}")
    #IF ( APPLE )
    #    IF ( NOT XCODE ) # Test if not xcode
    #        IF ( NOT ENABLE_PYTHON AND NOT USE_VTK )
    #            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
    #            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
    #        ENDIF ( NOT ENABLE_PYTHON AND NOT USE_VTK )
    #    ENDIF ( NOT XCODE  )
    #ENDIF ( APPLE )
ENDIF()

IF(APPLE)
    OPTION(BUILD_UNIVERSAL "Build Universal Binaries" OFF)
    IF(BUILD_UNIVERSAL)
        SET(GCC_UNIVERSAL_FLAGS "-arch i386 -arch x86_64")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_UNIVERSAL_FLAGS}")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_UNIVERSAL_FLAGS}")
    ENDIF()
ENDIF()


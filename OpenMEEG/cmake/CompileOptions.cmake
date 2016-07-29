#######################################################################
# Setting compilation options
#######################################################################

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(USE_GCC YES)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(USE_ICC YES)
endif()

if (UNIX)
    option(FORCE_BUILD_32BITS "Force 32 bits compilation" OFF)
    if (FORCE_BUILD_32BITS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
    endif()

    mark_as_advanced(FORCE_BUILD_32BITS)

    if (NOT BUILD_SHARED_LIBS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()


if (USE_GCC)
    set(CXX_WARNING_OPTIONS "-Wall -W -Wno-unknown-pragmas -Wshadow -Wunused-variable -Wunused-parameter -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings")
    set(CC_WARNING_OPTIONS "-Wall -W -Wno-unknown-pragmas -Wshadow -Wunused-variable -Wunused-parameter -Wunused -Wno-system-headers -Wwrite-strings")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CXX_WARNING_OPTIONS}")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${CC_WARNING_OPTIONS}")
endif()

if (APPLE)
    option(BUILD_UNIVERSAL "Build Universal Binaries" OFF)
    if (BUILD_UNIVERSAL)
        set(GCC_UNIVERSAL_FLAGS "-arch i386 -arch x86_64")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_UNIVERSAL_FLAGS}")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_UNIVERSAL_FLAGS}")
    endif()
endif()


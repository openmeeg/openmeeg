cmake_minimum_required (VERSION 3.14)

if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.31")
    cmake_policy(SET CMP0174 NEW)  # empty values different from not defined
endif()

project (OpenMEEG VERSION 2.5.16 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)  # use c++23 so that #warning remains portable

set(OpenMEEG_HEADER_INSTALLDIR include) # XXX I don't know why this should be here

# Modify where to find the modules
list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils
    ${CMAKE_MODULE_PATH}
)

# OpenMEEG Configuration and needed things

if (NOT WIN32)
    option(USE_MEMORY_CHECKER "Use memory checker" OFF)
    if (USE_MEMORY_CHECKER)
        set(CMAKE_CXX_FLAGS -fsanitize=address ${CMAKE_CXX_FLAGS})
    else()
        list(REMOVE_ITEM CMAKE_CXX_FLAGS "-fsanitize=address")
    endif()
endif()

include(GenerateExportHeader)
include(CompilerOptions)
include(options)
include(thirdParties)
include(install)  # install setup (do it before adding targets)
include(Uninstall)  # install setup (do it before adding targets)

message("-- Configuring OpenMEEGConfigure.h ...")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/OpenMEEGConfigure.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/OpenMEEGConfigure.h)

# OpenMEEG stuff

add_subdirectory(OpenMEEGMaths)
add_subdirectory(OpenMEEG)
if (ENABLE_APPS)
    add_subdirectory(apps)
endif()
add_subdirectory(wrapping)
if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
add_subdirectory(documentation)

# OpenMEEG Export

include(cmake/Packaging.cmake)

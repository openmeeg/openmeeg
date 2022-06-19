cmake_minimum_required (VERSION 3.12)
project (OpenMEEG VERSION 2.4.9999 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)  # use c++17

set(OpenMEEG_HEADER_INSTALLDIR include) # XXX I don't know why this should be here

# Modify where to find the modules
list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils
    ${CMAKE_MODULE_PATH}
)

# OpenMEEG Configuration and needed things

include(GenerateExportHeader)
include(CompilerOptions)
include(options)
include(thirdParties)
include(install)  # install setup (do it before adding targets)
include(Uninstall)  # install setup (do it before adding targets)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/OpenMEEGConfigure.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/OpenMEEGConfigure.h)

# OpenMEEG stuff

add_subdirectory(OpenMEEGMaths)
add_subdirectory(OpenMEEG)
add_subdirectory(apps)
add_subdirectory(documentation)
add_subdirectory(wrapping)
add_subdirectory(tests)

# OpenMEEG Export

include(cmake/Packaging.cmake)

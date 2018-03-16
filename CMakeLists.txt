cmake_minimum_required (VERSION 3.2)
project (OpenMEEG VERSION 2.4.0 LANGUAGES C CXX)

set(OpenMEEG_HEADER_INSTALLDIR include) # XXX I don't know why this should be here

# Modify where to find the modules
list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FortranCInterface
    ${CMAKE_MODULE_PATH}
)

################
# OpenMEEG Configuration and needed things
################
include(GenerateExportHeader)
include(cmake/CompilerOptions.cmake)
include(cmake/options.cmake)
include(cmake/thirdParties.cmake)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/OpenMEEGConfigure.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/OpenMEEGConfigure.h)

################
# OpenMEEG stuff
################
add_subdirectory(OpenMEEGMaths)
add_subdirectory(OpenMEEG)
add_subdirectory(apps)

################
# OpenMEEG Export
################
include(cmake/install.cmake)
include(cmake/Packaging.cmake)

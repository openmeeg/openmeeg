cmake_minimum_required (VERSION 3.9)
project (OpenMEEG VERSION 2.4.9999 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)  # use c++11

# Modify where to find the modules
list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FortranCInterface
)

################
# OpenMEEG Configuration and needed things
################
include(GenerateExportHeader)
include(CompilerOptions)
include(options)
include(thirdParties)
include(install)  # install setup (do it before adding targets)
include(Uninstall)  # install setup (do it before adding targets)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/OpenMEEGConfigure.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/OpenMEEGConfigure.h)

################
# OpenMEEG stuff
################
add_subdirectory(OpenMEEGMaths)
add_subdirectory(OpenMEEG)
add_subdirectory(apps)
add_subdirectory(documentation)
add_subdirectory(wrapping)

################
# OpenMEEG Export
################
include(Packaging)

# XXX: Hot fix to be able to ship win binaries when OpenMP is active
#
# see:https://braintrekking.wordpress.com/2013/04/27/dll-hell-how-to-include-microsoft-redistributable-runtime-libraries-in-your-cmakecpack-project/
# TODO:
#   - Build a test that fails if this is not active
#   - Find a way to find vcompXXX.dll in the system
#   - Place this in the proper place (probably the entire install.cmake/thirdparty.cmake needs refactor)
if (WIN32 AND OPENMP_FOUND AND VCOMP_WORKAROUND)
    set(VCOMP_PATH "C:/Program Files (x86)/Microsoft Visual Studio 14.0/VC/redist/x64/Microsoft.VC140.OpenMP/vcomp140.dll")
    INSTALL(FILES "${VCOMP_PATH}" DESTINATION bin COMPONENT Libraries)
endif()

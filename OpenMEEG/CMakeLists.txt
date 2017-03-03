################################################################
#
#    Compilation file for OpenMEEG library
#
################################################################

# CMakeLists files in this project can refer to the root source
# directory of the project as ${OpenMEEG_SOURCE_DIR} and
# to the root binary directory of the project as ${OpenMEEG_BINARY_DIR}.

#------------------------------------------------------------
# Project Name
#------------------------------------------------------------

set(PACKAGE_VERSION_MAJOR "2")
set(PACKAGE_VERSION_MINOR "3")
set(PACKAGE_VERSION_PATCH "dev")

#   The second path is temporary untill these are integrated in cmake...

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_SOURCE_DIR}/cmake/FortranCInterface)

#   The next two lines are a kludge to make rpm generation work. If you change this,
#   please check that rpm generation still works.

set(CPACK_PACKAGE_VERSION_MAJOR ${PACKAGE_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PACKAGE_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PACKAGE_VERSION_PATCH})

SET(VERSION_STRING "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

project(OpenMEEG C CXX)

cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

# By default, do not warn when built on machines using only VS Express:

if (NOT DEFINED CMAKE_install_SYSTEM_RUNTIME_LIBS_NO_WARNINGS)
    set(CMAKE_install_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON)
endif()

include(NewExecutable)
include(GenerateConfigFile)
include(FindLibrary)
include(Win32Compat)

include(CompileOptions)

set(OPENMEEG_HEADER_INSTALLDIR ${INSTALL_INCLUDE_DIR}/openmeeg)

# ==============
# = Endianess =
# ==============

include(TestBigEndian)
test_big_endian(WORDS_BIGENDIAN)

#######################################################################
# External Libraries
#######################################################################

include(UseOpenMP)
include(UseVTK)
include(UseGifti)
include(UseCGAL)
include(ProgressBar)
include(BlasLapackOption)
include(BlasLapack)

Dependencies(OpenMEEG matio)

#   We may want to make these warning fatal errors in a while.

if (NOT matio_FOUND)
    message(FATAL_ERROR "No matio found.")
endif()
if (matio_VERSION VERSION_LESS "1.5.0")
    message(WARNING "Found matio version: ${MATIO_VERSION}")
    message(FATAL_ERROR "OpenMEEG will not compile with matio versions less than 1.5.0.")
endif()

#---------------------------------------------------------------
# Include OpenMEEG Directories
#---------------------------------------------------------------

# Include the file with all source tree directories for OpenMEEG

set(OPENMEEG_INCLUDE_DIRECTORIES
    ${OpenMEEG_BINARY_DIR}
    ${OpenMEEG_SOURCE_DIR}
    ${OpenMEEG_SOURCE_DIR}/src
    ${OpenMEEG_SOURCE_DIR}/libs/OpenMEEGMaths/include
    ${OpenMEEG_SOURCE_DIR}/libs/OpenMEEG/include
    ${OpenMEEG_BINARY_DIR}/libs/OpenMEEG/include
    ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES}
    ${matio_INCLUDE_DIRS}
)

include_directories(${OPENMEEG_INCLUDE_DIRECTORIES})

#---------------------------------------------------------------
# Test the C++ compiler features
#---------------------------------------------------------------

include(CheckCXXFeatures)
CHECK_CXX_STANDARD_LIBRARY()

if (NOT HAVE_ISNORMAL_IN_NAMESPACE_STD)
    include(CheckSymbolExists)
    check_symbol_exists(isnormal math.h HAVE_ISNORMAL_IN_MATH_H)
endif()

#-----------------------------------------------
# tests
#-----------------------------------------------

include(TestOption)
include(Dart)

#---------------------------------------------------------------
# Setting sources and libraries
#---------------------------------------------------------------

add_subdirectory(libs)
add_subdirectory(src)

#######################################################################
# Configure files with settings for use by the build.
#######################################################################

configure_file(${OpenMEEG_SOURCE_DIR}/OpenMEEGConfigure.h.in
               ${OpenMEEG_BINARY_DIR}/OpenMEEGConfigure.h)

if (UNIX)
    include(AddDemoTarget)
    ADD_DEMO_TARGET(Head1 forinv)
    ADD_DEMO_TARGET(Head2 forward)
    ADD_DEMO_TARGET(Head3 forward)
    ADD_DEMO_TARGET(nerve1 gain)
endif()

#-----------------------------------------------
# tools
#-----------------------------------------------

option(BUILD_TOOLS "Build tools" ON)

if (BUILD_TOOLS)
    subdirs(tools)
endif()

#-----------------------------------------------
# documentation
#-----------------------------------------------

option(BUILD_DOCUMENTATION "Build Doxygen documentation" OFF)
option(BUILD_TUTORIALS "Build Tutorials" OFF)

add_subdirectory(documentation)

#-----------------------------------------------
# Wrapping
#-----------------------------------------------

option(ENABLE_PYTHON "Enable Python Wrapping" OFF)

if (ENABLE_PYTHON AND NOT BUILD_SHARED_LIBS)
    message(SEND_ERROR "Python wrapping does not work with static build, consider changing the BUILD_SHARED_LIBS option.")
endif ()

if (APPLE AND ENABLE_PYTHON)
    if (NOT FORCE_BUILD_32BITS) # Builtin Python2.5 on Mac OS X in 32 bits...
        message("You may have a 32bits Python installed.
                 In case of compilation problem try to force 32 bits compilation.")
    endif()
endif()

add_subdirectory(Wrapping)

#   Tests

if (BUILD_TESTING)
    include(CreateTestEnvironment)
    add_subdirectory(tests)
endif()

include(Uninstall)
include(KWStyle)
include(Packaging)

#######################################################################
# Help other projects use OpenMEEG
# DEFINED FOR OpenMEEGConfig.cmake.in file
#######################################################################

#-----------------------------------------------
# Creating files for external projects
#-----------------------------------------------

# Include the libraries and their directories for OpenMEEG

set(OPENMEEG_LIBRARY_DIRECTORIES
    ${OpenMEEG_BINARY_DIR}/src
    ${OpenMEEG_BINARY_DIR}/libs/OpenMEEGMaths/src
    ${OpenMEEG_BINARY_DIR}/libs/OpenMEEG/src
    ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES}
)

set(OPENMEEG_SOURCE_DIRECTORIES ${OpenMEEG_SOURCE_DIR})
set(OPENMEEG_LIBRARIES OpenMEEG OpenMEEGMaths ${OPENMEEG_LIBRARIES})

# This section uses the OPENMEEG_INCLUDE_DIRECTORIES variable.
# Copy the UseOpenMEEG.cmake file to the binary tree for backward compatability.

configure_file(
    ${OpenMEEG_SOURCE_DIR}/cmake/UseOpenMEEG.cmake.in
    ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake COPYONLY IMMEDIATE)

# The "use" file.

set(OPENMEEG_USE_FILE ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake)

# Configure OpenMEEGConfig.cmake for the build tree.

configure_file(
    ${OpenMEEG_SOURCE_DIR}/cmake/OpenMEEGConfig.cmake.in
    ${OpenMEEG_BINARY_DIR}/OpenMEEGConfig.cmake @ONLY IMMEDIATE)

# Generate another configure for installation

set(OPENMEEG_INCLUDE_DIRECTORIES
    ${OPENMEEG_HEADER_INSTALLDIR}
    ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES}
)
set(OPENMEEG_LIBRARY_DIRECTORIES
    ${CMAKE_install_PREFIX}/lib
    ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES}
)

# The "use" file.

set(OPENMEEG_USE_FILE ${CMAKE_install_PREFIX}/lib/UseOpenMEEG.cmake)

# Configure OpenMEEGConfig.cmake for the install tree.

configure_file(
    ${OpenMEEG_SOURCE_DIR}/cmake/OpenMEEGConfig.cmake.in
    ${OpenMEEG_BINARY_DIR}/cmake/OpenMEEGConfig.cmake @ONLY IMMEDIATE)

install(FILES ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake
    DESTINATION ${INSTALL_DATA_DIR}/OpenMEEG/cmake)

install(FILES ${OpenMEEG_BINARY_DIR}/cmake/OpenMEEGConfig.cmake
    DESTINATION ${INSTALL_DATA_DIR}/OpenMEEG/cmake)

# Install README

install(FILES ${PROJECT_SOURCE_DIR}/LICENSE.txt  ${PROJECT_SOURCE_DIR}/../README.rst
        DESTINATION ${INSTALL_DATA_DIR}/doc/OpenMEEG)

# install generated headers OpenMEEGConfigure and FC

install(FILES ${OpenMEEG_BINARY_DIR}/OpenMEEGConfigure.h
        DESTINATION ${OPENMEEG_HEADER_INSTALLDIR} COMPONENT Development)

if (USE_ATLAS)
    install(FILES ${OpenMEEG_BINARY_DIR}/FC.h
            DESTINATION ${OPENMEEG_HEADER_INSTALLDIR} COMPONENT Development)
endif()

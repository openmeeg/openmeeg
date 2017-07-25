################################################################
#
#    Compilation file for OpenMEEG library
#
################################################################

# CMakeLists files in this project can refer to the root source
# directory of the project as ${OpenMEEG_SOURCE_DIR} and
# to the root binary directory of the project as ${OpenMEEG_BINARY_DIR}.

cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
if (CMAKE_MAJOR_VERSION EQUAL 3)
    cmake_policy(SET CMP0048 NEW)
    cmake_policy(SET CMP0037 NEW)
    set(PROJECT_VERSION VERSION 2.3.99)
endif()
# The second path is temporary untill these are integrated in cmake...
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_SOURCE_DIR}/cmake/FortranCInterface)

#------------------------------------------------------------
# Project Name
#------------------------------------------------------------

project(OpenMEEG LANGUAGES C CXX ${PROJECT_VERSION})

# The next two lines are a kludge to make rpm generation work. If you change this,
# please check that rpm generation still works.

set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})

SET(VERSION_STRING "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

# By default, do not warn when built on machines using only VS Express:

if (NOT DEFINED CMAKE_install_SYSTEM_RUNTIME_LIBS_NO_WARNINGS)
    set(CMAKE_install_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON)
endif()

include(NewExecutable)
include(GenerateConfigFile)
include(FindLibrary)
include(Win32Compat)
include(CompileOptions)

set(OpenMEEG_HEADER_INSTALLDIR ${INSTALL_INCLUDE_DIR}/openmeeg)
set(OpenMEEG_DEPENDENCIES matio)

# ==============
# = Endianess =
# ==============

include(TestBigEndian)
test_big_endian(WORDS_BIGENDIAN)

#######################################################################
# External Libraries
#######################################################################

# Options
include(OpenMEEGOptions)

Dependencies(OpenMEEG ${OpenMEEG_DEPENDENCIES})

if (NOT matio_FOUND)
    message(FATAL_ERROR "No matio found.")
endif()
# check matio version. TODO remove this check in 2018 (matio 1.5.0 was out in 2012!)
if (matio_VERSION VERSION_LESS "1.5.0")
    message(WARNING "Found matio version: ${MATIO_VERSION}")
    message(FATAL_ERROR "OpenMEEG will not compile with matio versions less than 1.5.0.")
endif()

#---------------------------------------------------------------
# Include OpenMEEG Directories
#---------------------------------------------------------------

# Include the file with all source tree directories for OpenMEEG

set(OpenMEEG_INCLUDE_DIRECTORIES
    ${OpenMEEG_BINARY_DIR}
    ${OpenMEEG_SOURCE_DIR}
    ${OpenMEEG_SOURCE_DIR}/src
    ${OpenMEEG_SOURCE_DIR}/libs/OpenMEEGMaths/include
    ${OpenMEEG_SOURCE_DIR}/libs/OpenMEEG/include
    ${OpenMEEG_BINARY_DIR}/libs/OpenMEEG/include
    ${OpenMEEG_BINARY_DIR}/libs/OpenMEEGMaths/src
    ${OpenMEEG_BINARY_DIR}/libs/OpenMEEG/src
    ${OpenMEEG_OTHER_INCLUDE_DIRS}
    ${matio_INCLUDE_DIRS}
)

include_directories(${OpenMEEG_INCLUDE_DIRECTORIES})

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

if (BUILD_TOOLS)
    subdirs(tools)
endif()

#-----------------------------------------------
# documentation
#-----------------------------------------------

add_subdirectory(documentation)

#-----------------------------------------------
# Wrapping
#-----------------------------------------------
if (ENABLE_PYTHON AND NOT BUILD_SHARED_LIBS)
    message(SEND_ERROR "Python wrapping does not work with static build, "
        "consider changing the BUILD_SHARED_LIBS option.")
endif()

if (APPLE AND ENABLE_PYTHON)
    if (NOT FORCE_BUILD_32BITS) # Builtin Python2.5 on Mac OS X in 32 bits...
        message("You may have a 32bits Python installed."
                "In case of compilation problem try to force 32 bits compilation.")
    endif()
endif()

add_subdirectory(wrapping)

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
set(OpenMEEG_LIBRARIES OpenMEEG OpenMEEGMaths)

# Add all targets to the build-tree export set
export(TARGETS OpenMEEGMaths OpenMEEG 
       FILE "${PROJECT_BINARY_DIR}/OpenMEEGTargets.cmake")
 
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE OpenMEEG)

# Configure OpenMEEGConfig.cmake for the install tree.
GenerateConfigFile(OpenMEEG LIBRARIES ${OpenMEEG_LIBRARIES})

# Configure UseOpenMEEG.cmake for simple support.
configure_file(${OpenMEEG_SOURCE_DIR}/cmake/UseOpenMEEG.cmake.in
               ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake COPYONLY IMMEDIATE)

# Install the export set for use with the install-tree
install(EXPORT OpenMEEGTargets DESTINATION "${INSTALL_CMAKE_DIR}"
        COMPONENT DevelopmentCOMPONENT)

# Install README
install(FILES ${PROJECT_SOURCE_DIR}/../LICENSE.txt
        ${PROJECT_SOURCE_DIR}/../README.rst
        DESTINATION ${INSTALL_DATA_DIR}/doc/OpenMEEG)

# install generated headers OpenMEEGConfigure and FC
install(FILES ${OpenMEEG_BINARY_DIR}/OpenMEEGConfigure.h
        DESTINATION ${OpenMEEG_HEADER_INSTALLDIR}
        COMPONENT Development)

if (NOT USE_MKL AND NOT WIN32 AND NOT USE_OPENBLAS)
    install(FILES ${OpenMEEG_BINARY_DIR}/FC.h
            DESTINATION ${OpenMEEG_HEADER_INSTALLDIR}
            COMPONENT Development)
endif()

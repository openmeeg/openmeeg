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

SET(PACKAGE_VERSION_MAJOR "1")
SET(PACKAGE_VERSION_MINOR "1")
SET(PACKAGE_VERSION_PATCH "svn")

#   The next two lines are a kludge to make rpm generation work. If you change this,
#   please check that rpm generation still works.

SET(CPACK_PACKAGE_VERSION_MAJOR ${PACKAGE_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${PACKAGE_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${PACKAGE_VERSION_PATCH})

IF(APPLE OR WIN32)
    PROJECT(OpenMEEG C CXX)
ELSE()
    PROJECT(OpenMEEG C CXX Fortran)
ENDIF()

CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

INCLUDE(macros/NewExecutable.cmake)
INCLUDE(macros/FindLibrary.cmake)
INCLUDE(macros/UnitTest.cmake)
INCLUDE(macros/Revision.cmake)
INCLUDE(macros/Win32Compat.cmake)

GENERATE_VERSION_INFO()

INCLUDE(macros/CompileOptions.cmake)

# ==============
# = Endianess =
# ==============

INCLUDE(${CMAKE_ROOT}/Modules/TestBigEndian.cmake)
TEST_BIG_ENDIAN(WORDS_BIGENDIAN)

#######################################################################
# External Libraries
#######################################################################

INCLUDE(macros/VtkOption.cmake)
INCLUDE(macros/MatioOption.cmake)
INCLUDE(macros/UseAtlas.cmake)
INCLUDE(macros/UseOpenMP.cmake)

#######################################################################
# Configure files with settings for use by the build.
#######################################################################

CONFIGURE_FILE(${OpenMEEG_SOURCE_DIR}/OpenMEEGConfigure.h.in
               ${OpenMEEG_BINARY_DIR}/OpenMEEGConfigure.h)

#---------------------------------------------------------------
# Include OpenMEEG Directories
#---------------------------------------------------------------

# Include the file with all source tree directories for OpenMEEG

SET (OPENMEEG_INCLUDE_DIRECTORIES
    ${OpenMEEG_BINARY_DIR}
    ${OpenMEEG_SOURCE_DIR}
    ${OpenMEEG_SOURCE_DIR}/src
    ${OpenMEEG_SOURCE_DIR}/libs/MatLib/MatLib
)

INCLUDE_DIRECTORIES(${OPENMEEG_INCLUDE_DIRECTORIES})

#---------------------------------------------------------------
# Setting sources and libraries
#---------------------------------------------------------------

SUBDIRS(libs src)

IF(UNIX)
    INCLUDE(macros/AddDemoTarget.cmake)
    ADD_DEMO_TARGET(Head1 forinv)
    ADD_DEMO_TARGET(Head2 forward)
    ADD_DEMO_TARGET(Head3 forward)
    ADD_DEMO_TARGET(nerve1 gain)
ENDIF()

#-----------------------------------------------
# tools
#-----------------------------------------------

OPTION(BUILD_TOOLS "Build tools" ON)

IF (BUILD_TOOLS)
    SUBDIRS(tools)
ENDIF (BUILD_TOOLS)

#-----------------------------------------------
# tests
#-----------------------------------------------

INCLUDE(macros/TestOption.cmake)

#-----------------------------------------------
# documentation
#-----------------------------------------------

OPTION(BUILD_DOCUMENTATION "Build Doxygen documentation" OFF)
OPTION(BUILD_TUTORIALS "Build Tutorials" OFF)

SUBDIRS(${OpenMEEG_SOURCE_DIR}/documentation)

#-----------------------------------------------
# Wrapping
#-----------------------------------------------

OPTION(PYTHON_WRAP "Enable Python Wrapping" OFF)

IF (PYTHON_WRAP AND NOT BUILD_SHARED)
    MESSAGE(SEND_ERROR "Python wrapping does not work with static build, consider changing the BUILD_SHARED option.")
ENDIF ()

IF (PYTHON_WRAP AND APPLE AND USE_ICC) # Builtin Python2.5 on Mac OS X in 32 bits...
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
ENDIF ()

SUBDIRS(Wrapping)

#-----------------------------------------------
# Support for make uninstall
#-----------------------------------------------

CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/macros/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
    "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

#-----------------------------------------------
# Code Syntax Checking with KWStyle
#-----------------------------------------------

OPTION(CHECK_SYNTAX "Checking code syntax using KWStyle" OFF)
MARK_AS_ADVANCED(CHECK_SYNTAX)

IF ( CHECK_SYNTAX )
    ADD_CUSTOM_TARGET(check_syntax
        COMMAND KWStyle -xml ${OpenMEEG_SOURCE_DIR}/OpenMEEGConfig.kws.xml -html ${OpenMEEG_BINARY_DIR}/KWStyleCheck -D KWStyleFilesToCheck.txt -v
        WORKING_DIRECTORY ${OpenMEEG_SOURCE_DIR})
ENDIF ( CHECK_SYNTAX )

#-----------------------------------------------
# Progress bar
#-----------------------------------------------

OPTION(USE_PROGRESSBAR "Show ascii progress bar when assembling matrices" OFF)

IF ( USE_PROGRESSBAR )
    ADD_DEFINITIONS(-DUSE_PROGRESSBAR)
ENDIF ( USE_PROGRESSBAR )

#-----------------------------------------------
# packaging
#-----------------------------------------------

INCLUDE(macros/Packaging.cmake)

#######################################################################
# Help other projects use OpenMEEG
# DEFINED FOR OpenMEEGConfig.cmake.in file
#######################################################################

#-----------------------------------------------
# Creating files for external projects
#-----------------------------------------------

# This section uses the OPENMEEG_INCLUDE_DIRECTORIES variable.
# Copy the UseOpenMEEG.cmake file to the binary tree for backward compatability.
CONFIGURE_FILE(
    ${OpenMEEG_SOURCE_DIR}/macros/UseOpenMEEG.cmake.in
    ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake COPYONLY IMMEDIATE)

# The "use" file.
SET(OPENMEEG_USE_FILE ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake)

# Configure OpenMEEGConfig.cmake for the build tree.
CONFIGURE_FILE(
    ${OpenMEEG_SOURCE_DIR}/macros/OpenMEEGConfig.cmake.in
    ${OpenMEEG_BINARY_DIR}/OpenMEEGConfig.cmake @ONLY IMMEDIATE)

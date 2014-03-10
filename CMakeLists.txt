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

SET(PACKAGE_VERSION_MAJOR "2")
SET(PACKAGE_VERSION_MINOR "3")
SET(PACKAGE_VERSION_PATCH "dev")

#   The second path is temporary untill these are integrated in cmake...

# SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/macros;${CMAKE_SOURCE_DIR}/macros/FortranCInterface;")
SET(CMAKE_MODULE_PATH
"${CMAKE_SOURCE_DIR}/macros;${CMAKE_SOURCE_DIR}/macros/FortranCInterface")

#   The next two lines are a kludge to make rpm generation work. If you change this,
#   please check that rpm generation still works.

SET(CPACK_PACKAGE_VERSION_MAJOR ${PACKAGE_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${PACKAGE_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${PACKAGE_VERSION_PATCH})

PROJECT(OpenMEEG C CXX)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

# By default, do not warn when built on machines using only VS Express:

IF (NOT DEFINED CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS)
    SET(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON)
ENDIF()

INCLUDE(NewExecutable)
INCLUDE(LibCorrections)
INCLUDE(FindLibrary)
INCLUDE(UnitTest)
INCLUDE(Revision)
INCLUDE(Win32Compat)

GENERATE_VERSION_INFO()

INCLUDE(CompileOptions)

# ==============
# = Endianess =
# ==============

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(WORDS_BIGENDIAN)

#######################################################################
# External Libraries
#######################################################################

INCLUDE(VtkOption)
INCLUDE(GiftiOption)
INCLUDE(UseMatio)
INCLUDE(UseAtlas)
INCLUDE(UseOpenMP)

#   Detect Fortran to C interface.
IF (USE_ATLAS)
    INCLUDE(FortranCInterface)
    FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_" FROM_LIBRARY blas[daxpy])
ENDIF()

#---------------------------------------------------------------
# Include OpenMEEG Directories
#---------------------------------------------------------------

# Include the file with all source tree directories for OpenMEEG
SET (OPENMEEG_INCLUDE_DIRECTORIES
    ${OpenMEEG_BINARY_DIR}
    ${OpenMEEG_SOURCE_DIR}
    ${OpenMEEG_SOURCE_DIR}/src
    ${OpenMEEG_SOURCE_DIR}/libs/OpenMEEGMaths/include
    ${OpenMEEG_SOURCE_DIR}/libs/OpenMEEG/include
    ${OpenMEEG_BINARY_DIR}/libs/OpenMEEG/include
    ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES}
)

INCLUDE_DIRECTORIES(${OPENMEEG_INCLUDE_DIRECTORIES})

#---------------------------------------------------------------
# Test the C++ compiler features
#---------------------------------------------------------------

INCLUDE(CheckCXXFeatures)
CHECK_CXX_STANDARD_LIBRARY()

if (NOT HAVE_ISNORMAL_IN_NAMESPACE_STD)
    include(CheckSymbolExists)
    CHECK_SYMBOL_EXISTS(isnormal math.h HAVE_ISNORMAL_IN_MATH_H)
endif()

#---------------------------------------------------------------
# Setting sources and libraries
#---------------------------------------------------------------

SUBDIRS(libs src)

#######################################################################
# Configure files with settings for use by the build.
#######################################################################

CONFIGURE_FILE(${OpenMEEG_SOURCE_DIR}/OpenMEEGConfigure.h.in
               ${OpenMEEG_BINARY_DIR}/OpenMEEGConfigure.h)

IF(UNIX)
    INCLUDE(AddDemoTarget)
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

INCLUDE(TestOption)

#-----------------------------------------------
# Dart
#-----------------------------------------------

INCLUDE(Dart)

#-----------------------------------------------
# documentation
#-----------------------------------------------

OPTION(BUILD_DOCUMENTATION "Build Doxygen documentation" OFF)
OPTION(BUILD_TUTORIALS "Build Tutorials" OFF)

SUBDIRS(${OpenMEEG_SOURCE_DIR}/documentation)

#-----------------------------------------------
# Wrapping
#-----------------------------------------------

OPTION(ENABLE_PYTHON "Enable Python Wrapping" OFF)

IF (ENABLE_PYTHON AND NOT BUILD_SHARED_LIBS)
    MESSAGE(SEND_ERROR "Python wrapping does not work with static build, consider changing the BUILD_SHARED_LIBS option.")
ENDIF ()

IF(APPLE AND ENABLE_PYTHON)
    IF (NOT FORCE_BUILD_32BITS) # Builtin Python2.5 on Mac OS X in 32 bits...
        MESSAGE("You may have a 32bits Python installed.
                 In case of compilation problem try to force 32 bits compilation.")
    ENDIF()
ENDIF()

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

IF (CHECK_SYNTAX)
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

INCLUDE(Packaging)

#######################################################################
# Help other projects use OpenMEEG
# DEFINED FOR OpenMEEGConfig.cmake.in file
#######################################################################

#-----------------------------------------------
# Creating files for external projects
#-----------------------------------------------

# Include the libraries and their directories for OpenMEEG

SET (OPENMEEG_LIBRARY_DIRECTORIES
    ${OpenMEEG_BINARY_DIR}/src
    ${OpenMEEG_BINARY_DIR}/libs/OpenMEEGMaths/src
    ${OpenMEEG_BINARY_DIR}/libs/OpenMEEG/src
    ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES}
)

SET(OPENMEEG_SOURCE_DIRECTORIES ${OpenMEEG_SOURCE_DIR})

SET (OPENMEEG_LIBRARIES OpenMEEG OpenMEEGMaths ${OPENMEEG_LIBRARIES})

# This section uses the OPENMEEG_INCLUDE_DIRECTORIES variable.
# Copy the UseOpenMEEG.cmake file to the binary tree for backward compatability.
CONFIGURE_FILE(
    ${OpenMEEG_SOURCE_DIR}/macros/UseOpenMEEG.cmake.in
    ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake COPYONLY IMMEDIATE)

# The "use" file.
SET(OPENMEEG_USE_FILE ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake)

SET(OPENMEEG_HEADER_INSTALLDIR include/openmeeg)

# Configure OpenMEEGConfig.cmake for the build tree.
CONFIGURE_FILE(
    ${OpenMEEG_SOURCE_DIR}/macros/OpenMEEGConfig.cmake.in
    ${OpenMEEG_BINARY_DIR}/OpenMEEGConfig.cmake @ONLY IMMEDIATE)

# Generate another configure for installation
SET (OPENMEEG_INCLUDE_DIRECTORIES
    ${CMAKE_INSTALL_PREFIX}/${OPENMEEG_HEADER_INSTALLDIR}
    ${OPENMEEG_OTHER_INCLUDE_DIRECTORIES}
)
SET (OPENMEEG_LIBRARY_DIRECTORIES
    ${CMAKE_INSTALL_PREFIX}/lib
    ${OPENMEEG_OTHER_LIBRARY_DIRECTORIES}
)

# The "use" file.
SET(OPENMEEG_USE_FILE ${CMAKE_INSTALL_PREFIX}/lib/UseOpenMEEG.cmake)

# Configure OpenMEEGConfig.cmake for the install tree.
CONFIGURE_FILE(
    ${OpenMEEG_SOURCE_DIR}/macros/OpenMEEGConfig.cmake.in
    ${OpenMEEG_BINARY_DIR}/macros/OpenMEEGConfig.cmake @ONLY IMMEDIATE)

INSTALL(FILES ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake
    DESTINATION lib)

INSTALL(FILES ${OpenMEEG_BINARY_DIR}/macros/OpenMEEGConfig.cmake
    DESTINATION lib)

# install generated headers OpenMEEGConfigure and FC
INSTALL(FILES ${OpenMEEG_BINARY_DIR}/OpenMEEGConfigure.h
    DESTINATION ${OPENMEEG_HEADER_INSTALLDIR} COMPONENT Development)

IF(USE_ATLAS)
    INSTALL(FILES ${OpenMEEG_BINARY_DIR}/FC.h
        DESTINATION ${OPENMEEG_HEADER_INSTALLDIR} COMPONENT Development)
ENDIF()

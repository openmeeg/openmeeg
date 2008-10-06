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
SET(PACKAGE_VERSION_MINOR "0")
SET(PACKAGE_VERSION_PATCH "svn")

#   The next two lines are a kludge to make rpm generation work. If you change this,
#   please check that rpm generation still works.

SET(CPACK_PACKAGE_VERSION_MAJOR ${PACKAGE_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${PACKAGE_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${PACKAGE_VERSION_PATCH})

PROJECT(OpenMEEG)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

#SET(CMAKE_SKIP_RPATH ON)

#######################################################################
# Setting compilation options
#######################################################################

OPTION(BUILD_SHARED "Build shared libs" ON)
MARK_AS_ADVANCED(BUILD_SHARED)

IF ( BUILD_SHARED )
    SET(LIB_TYPE SHARED)
ELSE ( BUILD_SHARED )
    SET(LIB_TYPE STATIC)
ENDIF ( BUILD_SHARED )

IF(${CMAKE_C_COMPILER} MATCHES "gcc")
    SET(USE_GCC YES)
ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

IF(${CMAKE_C_COMPILER} MATCHES "icc")
    SET(USE_ICC YES)
ENDIF(${CMAKE_C_COMPILER} MATCHES "icc")

IF ( USE_GCC )
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
    IF ( APPLE )
        IF ( NOT XCODE ) # Test if not xcode
            IF ( NOT PYTHON_WRAP AND NOT USE_VTK)
                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
                SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
            ENDIF ( NOT PYTHON_WRAP AND NOT USE_VTK)
        ENDIF ( NOT XCODE  )
    ENDIF ( APPLE )
ENDIF ( USE_GCC )

# ==============================
# = Fix DLL search pb on WIN32 =
# ==============================

IF ( WIN32 )
    SET(LIBRARY_OUTPUT_PATH "${OpenMEEG_BINARY_DIR}/src")
    SET(EXECUTABLE_OUTPUT_PATH "${OpenMEEG_BINARY_DIR}/src")
ENDIF ( WIN32 )

# ==============
# = Endianess =
# ==============

INCLUDE(${CMAKE_ROOT}/Modules/TestBigEndian.cmake)
TEST_BIG_ENDIAN(WORDS_BIGENDIAN)

#######################################################################
# External Libraries
#######################################################################

#------------------------------------------------------------
# VTK library
#------------------------------------------------------------

OPTION(USE_VTK "Build the project using VTK" OFF)

IF (USE_VTK)
    FIND_PACKAGE(VTK)
    IF (VTK_FOUND)
        SET(USE_VTK 1)
        INCLUDE(${VTK_USE_FILE})
        SET (VTK_LIBRARIES
                # vtkRendering
                vtkGraphics
                # vtkHybrid
                # vtkImaging
                vtkIO
                # vtkFiltering
                # vtkGenericFiltering
                vtkCommon
                # vtkDICOMParser
                # vtkzlib
        )
        SET(OPENMEEG_OTHER_LIBRARIES
            ${OPENMEEG_OTHER_LIBRARIES}
            ${VTK_LIBRARIES}
        )
    ELSE(VTK_FOUND)
        MESSAGE(FATAL_ERROR "Please set VTK_DIR")
    ENDIF(VTK_FOUND)

ENDIF (USE_VTK)

SET(OPENMEEG_OTHER_LIBRARIES
        ${OPENMEEG_OTHER_LIBRARIES}
        OpenMEEGMaths
)

#------------------------------------------------------------
# MATIO library
#------------------------------------------------------------

OPTION(USE_MATIO "Build the project using MATIO" OFF)

IF (USE_MATIO)

    SET(SEARCHPATH
        /usr/lib/
        /usr/lib64/
        /usr/local/lib/atlas
        $ENV{HOME}/local/lib
    )

    FIND_LIBRARY(MATIO_LIB
                 NAMES matio
                 PATHS ${SEARCHPATH}
    )

    FIND_LIBRARY(ZLIB_LIB
                 NAMES z
                 PATHS ${SEARCHPATH}
    )

    SET(OPENMEEG_OTHER_LIBRARIES
        ${OPENMEEG_OTHER_LIBRARIES}
        ${MATIO_LIB}
        ${ZLIB_LIB}
    )

    FIND_PATH(MATIO_INCLUDE_PATH matio.h
              /usr/include/
              /usr/local/include/
              $ENV{HOME}/local/include/
    )
    INCLUDE_DIRECTORIES(${MATIO_INCLUDE_PATH})

ENDIF (USE_MATIO)


#######################################################################
# OpenMEEG
#######################################################################

#------------------------------------------------------------
# Atlas - MKL for blas and lapack
#------------------------------------------------------------

INCLUDE(UseAtlas.cmake)

#------------------------------------------------------------
# OpenMP
#------------------------------------------------------------

INCLUDE(UseOpenMP.cmake)

#------------------------------------------------------------
# MatLib
#------------------------------------------------------------

SET(MATLIB_PATH ${OpenMEEG_SOURCE_DIR}/libs/MatLib/MatLib)
INCLUDE_DIRECTORIES(${MATLIB_PATH})
SET(MATLIB_HEADERS
    LinOp.h
    MatLibConfig.h
    sparse_matrix.h
    fast_sparse_matrix.h
    generic_matrix.h
    matrix.h
    matrix.h
    sparse_matrix.h
    symmatrix.h
    symmatrix.h
    vector.h
    vector.h)
FOREACH ( HEADER ${MATLIB_HEADERS} )
    SET(MATLIB_FILES "${MATLIB_FILES} ${MATLIB_PATH}/${HEADER}")
ENDFOREACH ( HEADER )

#---------------------------------------------------------------
# Configure files with settings for use by the build.
#---------------------------------------------------------------
CONFIGURE_FILE (${OpenMEEG_SOURCE_DIR}/OpenMEEGConfigure.h.in
               ${OpenMEEG_BINARY_DIR}/OpenMEEGConfigure.h)

#---------------------------------------------------------------
# Include OpenMEEG Directories
#---------------------------------------------------------------

# Include the file with all source tree directories for OpenMEEG
#

SET (OPENMEEG_INCLUDE_DIRECTORIES
    ${OpenMEEG_BINARY_DIR}
    ${OpenMEEG_SOURCE_DIR}
    ${OpenMEEG_SOURCE_DIR}/src
)

INCLUDE_DIRECTORIES(
    ${OPENMEEG_INCLUDE_DIRECTORIES}
)

#---------------------------------------------------------------
# Setting sources and libraries
#---------------------------------------------------------------

SUBDIRS(libs src)

IF ( UNIX )
    ADD_CUSTOM_TARGET(copy_head1
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/Models/Head1
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/IntermediateFiles/Head1
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/Computations/Head1
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/Models/Head1 ${CMAKE_BINARY_DIR}/data/Models
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/IntermediateFiles/Head1 ${CMAKE_BINARY_DIR}/data/IntermediateFiles
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/Computations/Head1 ${CMAKE_BINARY_DIR}/data/Computations
    )
    ADD_CUSTOM_TARGET(copy_head2
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/Models/Head2
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/IntermediateFiles/Head2
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/Computations/Head2
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/Models/Head2 ${CMAKE_BINARY_DIR}/data/Models
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/IntermediateFiles/Head2 ${CMAKE_BINARY_DIR}/data/IntermediateFiles
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/Computations/Head2 ${CMAKE_BINARY_DIR}/data/Computations
    )
    ADD_CUSTOM_TARGET(copy_head3
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/Models/Head3
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/IntermediateFiles/Head3
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/Computations/Head3
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/Models/Head3 ${CMAKE_BINARY_DIR}/data/Models
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/IntermediateFiles/Head3 ${CMAKE_BINARY_DIR}/data/IntermediateFiles
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/Computations/Head3 ${CMAKE_BINARY_DIR}/data/Computations
    )
    ADD_CUSTOM_TARGET(demo_head1
        COMMAND make forinv
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/data/Computations/Head1
    )
    ADD_CUSTOM_TARGET(demo_head2
        COMMAND make forward
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/data/Computations/Head2
    )
    ADD_CUSTOM_TARGET(demo_head3
        COMMAND make forward
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/data/Computations/Head3
    )
    ADD_DEPENDENCIES(demo_head1 copy_head1)
    ADD_DEPENDENCIES(demo_head2 copy_head2)
    ADD_DEPENDENCIES(demo_head3 copy_head3)
ENDIF ( UNIX )

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

OPTION(BUILD_TEST "Build tests" OFF)

IF (BUILD_TEST)
    INCLUDE(CTest)
    ENABLE_TESTING()
    SUBDIRS(tests)
ENDIF (BUILD_TEST)

#-----------------------------------------------
# documentation
#-----------------------------------------------

OPTION(BUILD_DOCUMENTATION "Build Doxygen documentation" OFF)

IF ( BUILD_DOCUMENTATION )
    SET(OPENMEEG_DOXYGEN_HOME ${OpenMEEG_SOURCE_DIR}/documentation)
    SUBDIRS(${OPENMEEG_DOXYGEN_HOME})
ENDIF ( BUILD_DOCUMENTATION )

#-----------------------------------------------
# Wrapping
#-----------------------------------------------

OPTION(PYTHON_WRAP "Enable TCL Wrapping" OFF)

IF ( PYTHON_WRAP AND APPLE AND USE_ICC ) # Builtin Python2.5 on Mac OS X in 32 bits...
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
ENDIF ( PYTHON_WRAP AND APPLE AND USE_ICC )

SUBDIRS(Wrapping)

#-----------------------------------------------
# Code Syntax Checking with KWStyle
#-----------------------------------------------

OPTION(CHECK_SYNTAX "Checking code syntax using KWStyle" OFF)

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

OPTION(ENABLE_PACKAGING "Enable Packaging" ON)

IF(ENABLE_PACKAGING)

    INCLUDE(InstallRequiredSystemLibraries)

    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMEEG Project")
    SET(CPACK_PACKAGE_VENDOR "INRIA - Odyssee ENPC/INRIA/Ens Ulm")
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${OpenMEEG_SOURCE_DIR}/README.txt")
    SET(CPACK_RESOURCE_FILE_LICENSE "${OpenMEEG_SOURCE_DIR}/LICENSE.txt")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMEEG")

    IF(APPLE OR WIN32)

        IF(WIN32 AND NOT UNIX)
          # There is a bug in NSI that does not handle full unix paths properly. Make
          # sure there is at least one set of four (4) backlasshes.
          SET(CPACK_NSIS_DISPLAY_NAME "OpenMEEG Project")
          SET(CPACK_NSIS_HELP_LINK "https:\\\\\\\\gforge.inria.fr/projects/openmeeg/")
          SET(CPACK_NSIS_URL_INFO_ABOUT "https:\\\\\\\\gforge.inria.fr/projects/openmeeg/")
          SET(CPACK_NSIS_CONTACT "openmeeg-info@lists.gforge.inria.fr")
          SET(CPACK_NSIS_MODIFY_PATH ON)
        ELSE(WIN32 AND NOT UNIX)
          SET(CPACK_SOURCE_STRIP_FILES "")
        ENDIF(WIN32 AND NOT UNIX)
    ENDIF(APPLE OR WIN32)

    INCLUDE(CPack)

    OPTION(BUILD_RPM "Build RPM package" OFF)
    IF(UNIX AND BUILD_RPM) # linux
        INCLUDE(UseRPMTools)
        IF (RPMTools_FOUND)
            RPMTools_ADD_RPM_TARGETS(${PROJECT_NAME} "packaging/${PROJECT_NAME}.spec.in")
        ENDIF (RPMTools_FOUND)
    ENDIF(UNIX AND BUILD_RPM)

ENDIF(ENABLE_PACKAGING)

#######################################################################
# Help other projects use OpenMEEG
# DEFINED FOR OpenMEEGConfig.cmake.in file
#######################################################################

#-----------------------------------------------
# Include directories
#-----------------------------------------------
SET (OPENMEEG_INCLUDE_DIRECTORIES
        ${OPENMEEG_INCLUDE_DIRECTORIES}
)

#-----------------------------------------------
# Creating files for external projects
#-----------------------------------------------

# Copy the UseOpenMEEG.cmake file to the binary tree for backward compatability.
CONFIGURE_FILE(
    ${OpenMEEG_SOURCE_DIR}/UseOpenMEEG.cmake.in
    ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake COPYONLY IMMEDIATE)

# The "use" file.
SET(OPENMEEG_USE_FILE ${OpenMEEG_BINARY_DIR}/UseOpenMEEG.cmake)

# Configure OpenMEEGConfig.cmake for the build tree.
CONFIGURE_FILE(
    ${OpenMEEG_SOURCE_DIR}/OpenMEEGConfig.cmake.in
    ${OpenMEEG_BINARY_DIR}/OpenMEEGConfig.cmake @ONLY IMMEDIATE)

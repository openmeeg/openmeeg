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

SET(PACKAGE_VERSION_MAJOR "0")
SET(PACKAGE_VERSION_MINOR "8")
SET(PACKAGE_VERSION_PATCH "svn")

#   The next two lines are a kludge to make rpm generation work. If you change this,
#   please check that rpm generation still works.

SET(CPACK_PACKAGE_VERSION_MAJOR ${PACKAGE_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${PACKAGE_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${PACKAGE_VERSION_PATCH})

PROJECT(OpenMEEG)

ADD_DEFINITIONS(-D USE_CMAKE)

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
                # vtkGraphics
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

OPTION(USE_OMP "Use OpenMP" OFF)

OPTION(USE_ICC "Use Intel Compiler" OFF)

IF(USE_OMP)
    IF(UNIX)
        IF(USE_ICC)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -openmp")
        ELSE(USE_ICC)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
        ENDIF(USE_ICC)
    ELSE(UNIX)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp")
    ENDIF(UNIX)
    ADD_DEFINITIONS(-D USE_OMP)
ENDIF(USE_OMP)


#------------------------------------------------------------
# MatLib
#------------------------------------------------------------

INCLUDE_DIRECTORIES(${OpenMEEG_SOURCE_DIR}/libs/MatLib/MatLib)

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

SUBDIRS(src)

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
    ENABLE_TESTING()
    SUBDIRS(tests)
    SUBDIRS(libs)
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
# packaging
#-----------------------------------------------

OPTION(ENABLE_PACKAGING "Enable Packaging" ON)

IF(ENABLE_PACKAGING)

    INCLUDE(CPack)

    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMEEG Project")
    SET(CPACK_PACKAGE_VENDOR "Odyssee ENPC/INRIA/Ens Ulm")
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${OpenMEEG_SOURCE_DIR}/README")
    SET(CPACK_RESOURCE_FILE_LICENSE "${OpenMEEG_SOURCE_DIR}/Copyright.txt")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMEEG")

    IF(APPLE OR WIN32)

        INCLUDE(InstallRequiredSystemLibraries)

        IF(WIN32 AND NOT UNIX)
          # There is a bug in NSI that does not handle full unix paths properly. Make
          # sure there is at least one set of four (4) backlasshes.
          SET(CPACK_NSIS_DISPLAY_NAME "OpenMEEG Project")
          SET(CPACK_NSIS_HELP_LINK "https:\\\\\\\\gforge.inria.fr/projects/openmeeg/")
          SET(CPACK_NSIS_URL_INFO_ABOUT "https:\\\\\\\\gforge.inria.fr/projects/openmeeg/")
          SET(CPACK_NSIS_CONTACT "odyssee@sophia.inria.fr")
          SET(CPACK_NSIS_MODIFY_PATH ON)
        ELSE(WIN32 AND NOT UNIX)
          SET(CPACK_SOURCE_STRIP_FILES "")
        ENDIF(WIN32 AND NOT UNIX)

    ELSE(APPLE OR WIN32) # linux

        OPTION(BUILD_RPM "Build RPM package" OFF)
        IF (BUILD_RPM)
            INCLUDE(UseRPMTools)
            IF (RPMTools_FOUND)
                RPMTools_ADD_RPM_TARGETS(${PROJECT_NAME} "packaging/${PROJECT_NAME}.spec.in")
            ENDIF (RPMTools_FOUND)
        ENDIF (BUILD_RPM)

    ENDIF(APPLE OR WIN32)

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

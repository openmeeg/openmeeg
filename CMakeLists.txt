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


IF (WIN32)

    OPTION(USE_MKL "Build the project with MKL" ON)

    IF(USE_MKL)

        FIND_PATH(MKL_INCLUDE_PATH mkl.h
                    "C:/Program Files/Intel/MKL/8.1.1/include"
                    ../../mkl/include
        )
        INCLUDE_DIRECTORIES(${MKL_INCLUDE_PATH})

        FIND_PATH(MKL_LIB_PATH libguide.lib
                    "C:/Program Files/Intel/MKL/8.1.1/ia32/lib"
                    ../../mkl/ia32/lib
        )
        LINK_DIRECTORIES(${MKL_LIB_PATH})

        SET(OPENMEEG_OTHER_LIBRARIES
            ${OPENMEEG_OTHER_LIBRARIES}
            mkl_solver
            mkl_c
            libguide
        )

    ENDIF(USE_MKL)

ELSE (WIN32) # Apple or Linux

    OPTION(USE_ATLAS "Build the project using ATLAS" ON)

    IF(USE_ATLAS)
        FIND_PATH(ATLAS_LDFLAGS libatlas.so
                    /usr/lib/sse2
                    /usr/lib64/atlas
                 )
        LINK_DIRECTORIES(${ATLAS_LDFLAGS})
        LINK_LIBRARIES(lapack cblas blas)
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
        INCLUDE_DIRECTORIES( /usr/include/atlas )
    ELSE(USE_ATLAS)
        OPTION(USE_MKL "Build the project using MKL" ON)
        IF(USE_MKL)
            FIND_PATH(MKL_INCLUDE_PATH mkl.h
                        ~/intel/mkl/8.1/include
                        ~/Intel/MKL/8.1/include
                        /Library/Frameworks/Intel_MKL.framework/Headers
            )
            INCLUDE_DIRECTORIES(${MKL_INCLUDE_PATH})

            FIND_PATH(MKL_LIB_PATH libguide.a
                        ~/intel/mkl/8.1/lib/32
                        ~/Intel/MKL/8.1/lib/32
                        /Library/Frameworks/Intel_MKL.framework/Libraries/32
            )
            LINK_DIRECTORIES(${MKL_LIB_PATH})

            SET(OPENMEEG_OTHER_LIBRARIES
                ${OPENMEEG_OTHER_LIBRARIES}
                mkl_lapack
                mkl_ia32
                guide
            )
        ENDIF(USE_MKL)
    ENDIF(USE_ATLAS)

ENDIF (WIN32)

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
    ${OpenMEEG_SOURCE_DIR}/libs/LsLib
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

        OPTION(BUILD_RPM "Build RPM package" ON)
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

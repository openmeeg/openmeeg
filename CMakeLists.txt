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

#PROJECT(OpenMEEG C CXX Fortran)
PROJECT(OpenMEEG)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

#######################################################################
# Setting compilation options
#######################################################################

SET(LIB_TYPE SHARED)

IF (NOT WIN32)
    OPTION(BUILD_SHARED "Build shared libs" ON)
    MARK_AS_ADVANCED(BUILD_SHARED)

    IF ( BUILD_SHARED )
        SET(LIB_TYPE SHARED)
    ELSE ( BUILD_SHARED )
        SET(LIB_TYPE STATIC)
    ENDIF ( BUILD_SHARED )

    STRING(COMPARE NOTEQUAL "${BUILD_SHARED_STATUS}" "" BUILD_SHARED_STATUS_NOT_EMPTY)
    IF(BUILD_SHARED_STATUS_NOT_EMPTY)
        STRING(COMPARE NOTEQUAL "${BUILD_SHARED_STATUS}" "${BUILD_SHARED}" RESET)
        IF(${RESET})
            # MESSAGE("Reset")
            IF(UNIX AND NOT APPLE AND USE_ATLAS)
                SET(LAPACK_ATLAS_LIB NOTFOUND CACHE STRING "Atlas Lib" FORCE)
                SET(LAPACK_LIB NOTFOUND CACHE STRING "Atlas Lib" FORCE)
                SET(CBLAS_LIB NOTFOUND CACHE STRING "Atlas Lib" FORCE)
                SET(ATLAS_LIB NOTFOUND CACHE STRING "Atlas Lib" FORCE)
                SET(BLAS_LIB NOTFOUND CACHE STRING "Atlas Lib" FORCE)
                SET(GFORTRAN_LIB NOTFOUND CACHE STRING "Atlas Lib" FORCE)
            ENDIF()
        ENDIF()
    ENDIF()

    # Store in cache previous value of BUILD_SHARED
    SET(BUILD_SHARED_STATUS "${BUILD_SHARED}" CACHE STRING "Build shared status" FORCE)
    MARK_AS_ADVANCED(BUILD_SHARED_STATUS)
ENDIF ()

IF(${CMAKE_C_COMPILER} MATCHES "gcc")
    SET(USE_GCC YES)
ENDIF(${CMAKE_C_COMPILER} MATCHES "gcc")

IF(${CMAKE_C_COMPILER} MATCHES "icc")
    SET(USE_ICC YES)
ENDIF(${CMAKE_C_COMPILER} MATCHES "icc")

IF ( USE_GCC )
    SET(GCC_WARNING_OPTIONS "-Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${GCC_WARNING_OPTIONS}")
    SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${GCC_WARNING_OPTIONS}")
    #IF ( APPLE )
    #    IF ( NOT XCODE ) # Test if not xcode
    #        IF ( NOT PYTHON_WRAP AND NOT USE_VTK )
    #            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
    #            SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
    #        ENDIF ( NOT PYTHON_WRAP AND NOT USE_VTK )
    #    ENDIF ( NOT XCODE  )
    #ENDIF ( APPLE )
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
MARK_AS_ADVANCED(USE_VTK)

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

    SET(MATIO_LIBRARIES ${MATIO_LIB} ${ZLIB_LIB})

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

#IF (WIN32)       # WINDOWS
#        OPTION(USE_ATLAS "Build the project using ATLAS" OFF)
#        OPTION(USE_MKL "Build the project with MKL" ON)
#        MARK_AS_ADVANCED(USE_MKL)
#        MARK_AS_ADVANCED(USE_ATLAS)
#ELSE (WIN32)
#    IF (APPLE)   # MACOSX
#        OPTION(USE_ATLAS "Build the project using ATLAS" ON)
#        OPTION(USE_MKL "Build the project with MKL" OFF)
#    ELSE (APPLE) # LINUX
#        OPTION(USE_ATLAS "Build the project using ATLAS" ON)
#        OPTION(USE_MKL "Build the project with MKL" OFF)
#    ENDIF (APPLE)
#ENDIF (WIN32)
#
#IF ( USE_ATLAS AND USE_MKL )
#    MESSAGE("Cannot set both USE_ATLAS and USE_MKL to ON")
#ENDIF ( USE_ATLAS AND USE_MKL )
#
#OPTION(USE_STATIC_BLAS_LAPACK "Use static BLAS LAPACK libraries for easy deployment" OFF)
#MARK_AS_ADVANCED(USE_STATIC_BLAS_LAPACK)
#IF ( BUILD_SHARED )
#    SET(BLA_STATIC FALSE)
#ELSE ( BUILD_SHARED )
#    SET(BLA_STATIC TRUE)
#ENDIF ( BUILD_SHARED )
#
#IF ( USE_MKL )
#    SET(BLAS_HEADER mkl.h)
#ENDIF ( USE_MKL )
#
#IF ( USE_ATLAS )
#    SET(BLAS_HEADER cblas.h)
#ENDIF ( USE_ATLAS )
#
#SET(BLAS_SEARCHPATH)
#
## Apple
#IF ( APPLE )
#    IF ( USE_MKL)
#        #SET(ENV{BLA_VENDOR} "Intel(mkl)")
#        SET(ENV{BLA_VENDOR} "Intel")
#        SET(BLAS_SEARCHPATH /Library/Frameworks/Intel_MKL.framework/Versions/Current/include/)
#    ELSE ( USE_MKL )
#        #SET(ENV{BLA_VENDOR} "Apple")
#        #SET(ENV{BLA_VENDOR} "Generic")
#        SET(ENV{BLA_VENDOR})
#        SET(BLAS_SEARCHPATH /System/Library/Frameworks/vecLib.framework/Versions/A/Headers/)
#    ENDIF ( USE_MKL )
#ENDIF ( APPLE )
#
## Linux
#IF ( UNIX AND NOT APPLE )
#    SET(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:/usr/lib64/atlas/")
#    IF ( USE_MKL )
#        SET(ENV{BLA_VENDOR} "Intel")
#        SET(BLAS_SEARCHPATH /usr/local/include /usr/include)
#    ELSE ( USE_MKL )
#        SET(ENV{BLA_VENDOR} "ATLAS")
#        SET(BLAS_HEADER atlas/cblas.h)
#        SET(BLAS_SEARCHPATH /usr/local/include /usr/include)
#    ENDIF ( USE_MKL )
#ENDIF ( UNIX AND NOT APPLE )
#
## Win32
#IF ( WIN32 )
#    IF ( USE_MKL )
#        SET(ENV{BLA_VENDOR} "Intel")
#        SET(BLAS_SEARCHPATH
#            "C:/Program Files/Intel/MKL/9.1.027/include"
#            "C:/Program Files/Intel/MKL/8.1.1/include"
#        )
#    ELSE ( USE_MKL )
#        MESSAGE("MKL is required on WIN32")
#    ENDIF ( USE_MKL )
#ENDIF ( WIN32 )
#
#FIND_PACKAGE(BLAS)
#IF (BLAS_FOUND)
#    #SET(OPENMEEG_OTHER_LIBRARIES
#    #    ${OPENMEEG_OTHER_LIBRARIES}
#    #    ${BLAS_LIBRARIES}
#    #)
#ELSE(BLAS_FOUND)
#    MESSAGE(FATAL_ERROR "Could not find BLAS")
#ENDIF(BLAS_FOUND)
#
#FIND_PACKAGE(LAPACK)
#IF (LAPACK_FOUND)
#    FIND_PATH(BLAS_INCLUDE_DIR ${BLAS_HEADER}
#        ${BLAS_SEARCHPATH}
#    )
#    INCLUDE_DIRECTORIES(${BLAS_INCLUDE_DIR})
#ELSE(LAPACK_FOUND)
#    MESSAGE(FATAL_ERROR "Could not find LAPACK")
#ENDIF(LAPACK_FOUND)

INCLUDE(macros/UseAtlas.cmake)

#------------------------------------------------------------
# OpenMP
#------------------------------------------------------------

INCLUDE(macros/UseOpenMP.cmake)

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
    ADD_CUSTOM_TARGET(copy_nerve1
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/Models/nerve1
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/IntermediateFiles/nerve1
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/Computations/nerve1
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/Models/nerve1 ${CMAKE_BINARY_DIR}/data/Models
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/Computations/nerve1 ${CMAKE_BINARY_DIR}/data/Computations
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
    ADD_CUSTOM_TARGET(demo_nerve1
        COMMAND make gain
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/data/Computations/nerve1
    )
    ADD_DEPENDENCIES(demo_head1 copy_head1)
    ADD_DEPENDENCIES(demo_head2 copy_head2)
    ADD_DEPENDENCIES(demo_head3 copy_head3)
    ADD_DEPENDENCIES(demo_nerve1 copy_nerve1)
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
    MARK_AS_ADVANCED(BUILD_TESTING)
    SUBDIRS(tests)
ENDIF (BUILD_TEST)

IF (USE_GCC AND BUILD_TEST)
    OPTION(ENABLE_COVERAGE "Enable coverage" OFF)
    MARK_AS_ADVANCED(ENABLE_COVERAGE)
ENDIF ( USE_GCC AND BUILD_TEST)

IF (ENABLE_COVERAGE)
    IF (USE_GCC)
        SET(CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs -ftest-coverage")
        SET(CMAKE_CXX_FLAGS "-g -O0 ${GCC_WARNING_OPTIONS} ${CMAKE_EXE_LINKER_FLAGS}")
        SET(CMAKE_C_FLAGS "-g -O0 -Wall -W ${CMAKE_EXE_LINKER_FLAGS}")
    ELSE (USE_GCC)
        MESSAGE(SEND_ERROR "Coverage is only available with gcc.")
    ENDIF (USE_GCC)
ENDIF (ENABLE_COVERAGE)

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

OPTION(ENABLE_PACKAGING "Enable Packaging" ON)

IF (UNIX AND NOT APPLE) # LINUX
    OPTION(BUILD_RPM "Enable RPM Packaging" OFF)
ENDIF (UNIX AND NOT APPLE)

IF(ENABLE_PACKAGING OR BUILD_RPM)

    INCLUDE(InstallRequiredSystemLibraries)

    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMEEG Project")
    SET(CPACK_PACKAGE_VENDOR "INRIA - Odyssee ENPC/INRIA/Ens Ulm")
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${OpenMEEG_SOURCE_DIR}/README.txt")
    SET(CPACK_RESOURCE_FILE_LICENSE "${OpenMEEG_SOURCE_DIR}/LICENSE.txt")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMEEG")
    SET(CPACK_PACKAGE_CONTACT "openmeeg-info_at_lists.gforge.inria.fr")

    IF(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64)
        SET(CPACK_RPM_PACKAGE_ARCHITECTURE x86_64)
    ELSEIF(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "i686")
        SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386)
        SET(CPACK_RPM_PACKAGE_ARCHITECTURE i386)
    ENDIF(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")

    IF ( APPLE )
        IF ( USE_OMP )
            SET(CPACK_PACKAGE_FILE_NAME
                "OpenMEEG-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-MacOSX-Intel-gcc42-OpenMP")
        ELSE ( USE_OMP )
            SET(CPACK_PACKAGE_FILE_NAME
                "OpenMEEG-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-MacOSX-Intel")
        ENDIF ( USE_OMP )
    ENDIF ( APPLE )

    IF ( UNIX AND NOT APPLE )
        SET(CPACK_PACKAGE_FILE_NAME
            "OpenMEEG-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-Linux.${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
    ENDIF ( UNIX AND NOT APPLE )
    IF ( WIN32 )
        SET(CPACK_PACKAGE_FILE_NAME
            "OpenMEEG-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-win32-x86")
    ENDIF ( WIN32 )

    IF (WIN32)
        # There is a bug in NSIS that does not handle full unix paths properly. Make
        # sure there is at least one set of four (4) backlasshes.
        SET(CPACK_NSIS_DISPLAY_NAME "OpenMEEG Project")
        SET(CPACK_NSIS_HELP_LINK "https:\\\\\\\\gforge.inria.fr/projects/openmeeg/")
        SET(CPACK_NSIS_URL_INFO_ABOUT "https:\\\\\\\\gforge.inria.fr/projects/openmeeg/")
        SET(CPACK_NSIS_CONTACT "openmeeg-info@lists.gforge.inria.fr")
        SET(CPACK_NSIS_MODIFY_PATH ON)
    ENDIF (WIN32)

    SET(CPACK_SOURCE_STRIP_FILES "")

    INCLUDE(CPack)

    IF(UNIX AND BUILD_RPM) # linux
        IF (CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION LESS 8)
            INCLUDE(UseRPMTools)
            IF (RPMTools_FOUND)
                RPMTools_ADD_RPM_TARGETS(${PROJECT_NAME} "packaging/${PROJECT_NAME}.spec.in")
            ENDIF (RPMTools_FOUND)
        ELSE ()
            SET(CPACK_RPM_USER_BINARY_SPECFILE "packaging/${PROJECT_NAME}.spec.in")
        ENDIF ()
    ENDIF(UNIX AND BUILD_RPM)

ENDIF(ENABLE_PACKAGING OR BUILD_RPM)

IF ( ENABLE_PACKAGING AND WIN32 )
    INCLUDE(macros/UseWin32dlls.cmake)
ENDIF ( ENABLE_PACKAGING AND WIN32 )

#######################################################################
# Help other projects use OpenMEEG
# DEFINED FOR OpenMEEGConfig.cmake.in file
#######################################################################

#-----------------------------------------------
# Include directories
#-----------------------------------------------
SET (OPENMEEG_INCLUDE_DIRECTORIES ${OPENMEEG_INCLUDE_DIRECTORIES})

#-----------------------------------------------
# Creating files for external projects
#-----------------------------------------------

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

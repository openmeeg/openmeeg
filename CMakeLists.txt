# OpenMEEG
#
# Copyright (c) INRIA 2005-2017. All rights reserved.
# See LICENSE.txt for details.
#
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.

cmake_minimum_required(VERSION 2.8.12)
if (CMAKE_MAJOR_VERSION EQUAL 3)
    cmake_policy(SET CMP0048 NEW)
    cmake_policy(SET CMP0037 OLD)
    set(PROJECT_VERSION VERSION 2.4.0.0)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/OpenMEEG/cmake
                      ${CMAKE_CURRENT_SOURCE_DIR}/OpenMEEG/cmake/FortranCInterface ${CMAKE_MODULE_PATH})

project(OpenMEEG-superProject ${PROJECT_VERSION})

# Setup everything MODULE_PATH, check environment, basic targets, options
include(Setup)

# Define revision (TODO: Get revision from the OpenMEEG repo).
if ("${CMAKE_VERSION}" VERSION_LESS "3.0")
    SetRevision(2 4 0 dev)
endif()

# Add packages

# There are three possibilities:
# - if the option USE_SYSTEM_${external_project} is set, then the system version
#   will be used
# - if the option USE_SYSTEM_${external_project} is not set:
#   * there is a subdir with the package name, in which case this subdir is used
#     directly
#   * otherwise use download and compile locally the package as an external
#     module.

# See matio/CMakeLists.txt for details.
# Provide stdint.h for windows, if it is not already there.
# MSVC starts to include stddint.h only from MSVC10 !!!

set(MSINTTYPES)
if (WIN32)
    include(CheckIncludeFiles)
    check_include_files(stdint.h HAVE_C99_STDINT)
    if (NOT HAVE_C99_STDINT)
        set(MSINTTYPES msinttypes)
    endif()
endif()

if (NOT LAPACK_LIBRARIES AND USE_LAPACK)
    set(LAPACK clapack)
    option(CLAPACK_BUILD_TESTING "Build CLAPACK tests" OFF)
endif()

if (USE_VTK)
    set(BUILD_VTK VTK)
    if (VTK_FOUND)
        set(USE_SYSTEM_VTK True CACHE BOOL "Use the VTK from the system" FORCE)
    endif()
endif()

if (USE_OPENBLAS AND WIN32)
    include(installOpenBLAS)
endif()

subprojects(${MSINTTYPES} zlib ${BUILD_VTK} hdf5 matio ${LAPACK} OpenMEEG)

include(Packaging)

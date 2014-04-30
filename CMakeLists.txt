# OpenMEEG
#
# Copyright (c) INRIA 2005-2014. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.

cmake_minimum_required(VERSION 2.8.9)

# Set CMAKE_BUILD_TYPE to Release by default.

if (DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING 
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
else()
    set(CMAKE_BUILD_TYPE Release CACHE STRING 
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

project(OpenMEEG-superProject)

option(GIT_HTTPS "Use https protocol to get git repositories." ON)
set(GIT_PREFIX "git@")
if (GIT_HTTPS)
    set(GIT_PREFIX "https://")
endif()

# Setup everything (MODULE_PATH, check environment, basic targets

include(Setup)

# Define revision (TODO: Get revision from teh OpenMEEG repo).

SetRevision(2 3 0 dev)

# Add packages

# There are three possibilities:
# - if the option USE_SYSTEM_${external_project} is set, then the system version  
#   will be used (not yet fully tested)
# - if the option USE_SYSTEM_${external_project} is not set:
#   * there is a subdir with the package name, in which case this subdir is used 
#     directly
#   * otherwise use download and compile locally the package as an external 
#     module.

# See matio/CMakeLists.txt for details.
# Provide stdint.h for windows, if it is not already there.
# MSVC starts to include stddint.h only from MSVC10 !!!

set(MSINTTYPES)
if (WIN32 AND NOT CYGWIN_INSTALL_PATH)
    set(MSINTTYPES msinttypes)
#    set(INTTYPES_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/msinttypes)
#    set(CMAKE_REQUIRED_INCLUDES ${INTTYPES_INCLUDES})
endif()

subprojects(${MSINTTYPES} zlib hdf5 matio clapack OpenMEEG)

#include(Packaging)

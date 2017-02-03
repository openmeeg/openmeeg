##############################################################################
# OpenMEEG
#
# Copyright (c) INRIA 2015-2017. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.
################################################################################

#   Get distribution name and architecture

execute_process(COMMAND lsb_release -a
                COMMAND grep "^Distributor ID:" 
                COMMAND sed -e "s/Distributor ID:[ \t]*//ig"
                OUTPUT_VARIABLE DISTRIBUTOR_ID
                OUTPUT_STRIP_TRAILING_WHITESPACE)
  
execute_process(COMMAND lsb_release -a
                COMMAND grep "^Release:"
                COMMAND sed -e "s/Release:[ \t]*//ig"
                OUTPUT_VARIABLE RELEASE
                OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND arch 
                OUTPUT_VARIABLE ARCH 
                OUTPUT_STRIP_TRAILING_WHITESPACE)
  
set(PACKAGE_ARCH "${DISTRIBUTOR_ID}_${RELEASE}-${ARCH}")
set(PACKAGE_ARCH_SHORT "Linux")
 
#   Set the right package generator

set(CPACK_GENERATOR DEB)
if(${DISTRIBUTOR_ID} MATCHES fc|fedora|Fedora|Centos|centos|SUSE|Suse|suse)
    set(CPACK_GENERATOR RPM)
endif()

set(CPACK_GENERATOR "TGZ;${CPACK_GENERATOR}")

#   Remember the linux packaging source dir

set(CURRENT_SRC_DIR ${CMAKE_SOURCE_DIR}/packaging/linux)
set(CURRENT_BIN_DIR ${CMAKE_BINARY_DIR}/packaging/linux)

#   include settings specific to DEB and RPM

include(${CURRENT_SRC_DIR}/RPM.cmake)
include(${CURRENT_SRC_DIR}/DEB.cmake)

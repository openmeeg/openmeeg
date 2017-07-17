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

execute_process(COMMAND cat /etc/os-release
                COMMAND grep "^ID="
                COMMAND sed -e "s/ID=[ \t]*//ig"
                OUTPUT_VARIABLE DISTRIBUTOR_ID
                OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND cat /etc/os-release
                COMMAND grep "^VERSION_ID="
                COMMAND sed -e "s/VERSION_ID=[ \t]*//ig"
                OUTPUT_VARIABLE RELEASE
                OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND arch
                OUTPUT_VARIABLE ARCH
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# on debian /etc/os-release: string contain "..."
string(REPLACE "\"" "" DISTRIBUTOR_ID ${DISTRIBUTOR_ID})
string(REPLACE "\"" "" RELEASE ${RELEASE})

set(PACKAGE_ARCH "${DISTRIBUTOR_ID}_${RELEASE}-${ARCH}")
set(PACKAGE_ARCH_SHORT "Linux")
 
#   Set the right package generator

set(CPACK_GENERATOR DEB)
if(${DISTRIBUTOR_ID} MATCHES fc|fedora|Fedora|Centos|centos|SUSE|Suse|suse)
    set(CPACK_GENERATOR RPM)
endif()

set(CPACK_GENERATOR "TGZ;${CPACK_GENERATOR}")

#   Remember the linux packaging source dir

set(CURRENT_SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CURRENT_BIN_DIR ${PROJECT_BINARY_DIR}/packaging/linux)

#   include settings specific to DEB and RPM

include(${CURRENT_SRC_DIR}/RPM.cmake)
include(${CURRENT_SRC_DIR}/DEB.cmake)

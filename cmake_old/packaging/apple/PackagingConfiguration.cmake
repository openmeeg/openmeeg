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

set(CPACK_BINARY_TGZ ON)

set(CPACK_BINARY_DRAGNDROP OFF)
set(CPACK_BINARY_PACKAGEMAKER OFF)

#   Get distribution name and architecture

set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_GENERATOR "PackageMaker;TGZ")

set(PACKAGE_ARCH "MacOSX")
set(PACKAGE_ARCH_SHORT "${PACKAGE_ARCH}")

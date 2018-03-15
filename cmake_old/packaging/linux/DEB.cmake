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

set(CPACK_DEBIAN_PACKAGE_HOMEPAGE http://openmeeg.github.io)
set(CPACK_DEBIAN_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
set(CPACK_DEBIAN_PACKAGE_PROVIDES ${CPACK_PACKAGE_NAME})
set(CPACK_DEBIAN_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION})
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "${CPACK_PACKAGE_DESCRIPTION}")

#set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA ${CMAKE_BINARY_DIR}/packaging/linux/prerm;${CMAKE_BINARY_DIR}/packaging/linux/postinst)

SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libmatio-dev (>= 1.5.0)")

if (ENABLE_PYTHON)
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, python-numpy (>= 1.6), swig (>= 2.0)")
endif()
if (USE_VTK)
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, libvtk5-dev (>= 5.8.0) | libvtk6-dev (>= 6.0.0) | libvtk7-dev (>= 7.0.0)")
endif()

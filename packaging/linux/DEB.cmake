##############################################################################
# OpenMEEG
#
# Copyright (c) INRIA 2015. All rights reserved.
# See LICENSE.txt for details.
# 
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.
################################################################################

set(CPACK_GENERATOR "${CPACK_GENERATOR};DEB")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE http://openmeeg.gforge.inria.fr)
set(CPACK_DEBIAN_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
set(CPACK_DEBIAN_PACKAGE_PROVIDES ${CPACK_PACKAGE_NAME})
set(CPACK_DEBIAN_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION})

#set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA ${CMAKE_BINARY_DIR}/packaging/linux/prerm;${CMAKE_BINARY_DIR}/packaging/linux/postinst)

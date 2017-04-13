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

set(CPACK_RPM_PACKAGE_PROVIDES "${CPACK_PACKAGE_NAME} = ${CPACK_PACKAGE_VERSION}")
set(CPACK_RPM_PACKAGE_LICENSE "CeCILL-B")
set(CPACK_RPM_PACKAGE_ARCHITECTURE ${ARCH})

set(CPACK_RPM_PACKAGE_DESCRIPTION "${CPACK_PACKAGE_DESCRIPTION}")
set(CPACK_RPM_PACKAGE_GROUP "Applications/Medical")
set(CPACK_RPM_CHANGELOG_FILE "${PROJECT_SOURCE_DIR}/cmake/packaging/changelog.txt")

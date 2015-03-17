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

set(CPACK_RPM_PACKAGE_PROVIDES "${CPACK_PACKAGE_NAME} = ${CPACK_PACKAGE_VERSION}")
set(CPACK_RPM_PACKAGE_LICENSE "CeCILL-B")
set(CPACK_RPM_PACKAGE_ARCHITECTURE ${ARCH})

set(CPACK_RPM_PACKAGE_DESCRIPTION  "OpenMEEG is a package for forward/inverse problems of EEG/MEG. The forward problem uses the symmetric Boundary Element Method. The inverse problem uses a distributed approach (L2, L1 regularization). Developped within Odyssee (INRIA-ENPC-ENS).")
set(CPACK_RPM_PACKAGE_GROUP "Applications/Medical")


#set(CPACK_RPM_POST_INSTALL_SCRIPT_FILE  ${CMAKE_BINARY_DIR}/packaging/linux/postinst)
#set(CPACK_RPM_PRE_UNINSTALL_SCRIPT_FILE ${CMAKE_BINARY_DIR}/packaging/linux/prerm)

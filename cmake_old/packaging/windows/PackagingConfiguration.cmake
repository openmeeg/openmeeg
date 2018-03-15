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

set(DISTRIB windows)
set(CPACK_GENERATOR "${CPACK_GENERATOR}" "NSIS")
set(CPACK_PACKAGE_TYPE NSIS)
set(CPACK_PACKAGING_INSTALL_PREFIX "")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
set(CPACK_MONOLITHIC_INSTALL 1)

if(CMAKE_CL_64)
	SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")
	#  - Text used in the installer GUI
	SET(CPACK_NSIS_PACKAGE_NAME "${CPACK_PACKAGE_NAME} (Win64)")
	#  - Registry key used to store info about the installation
	SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION} (Win64)")
	SET(PACKAGE_ARCH Win64)
else()
    SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
	SET(CPACK_NSIS_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
	SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}")
	SET(PACKAGE_ARCH Win32)
endif()
set(PACKAGE_ARCH_SHORT "${PACKAGE_ARCH}")

#   There is a bug in NSIS that does not handle full unix paths properly. Make
#   sure there is at least one set of four (4) backlasshes.

set(CPACK_NSIS_DISPLAY_NAME "OpenMEEG")
set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\openmeeg.github.io")
set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\openmeeg.github.io")
set(CPACK_NSIS_CONTACT ${CPACK_PACKAGE_CONTACT})

#   Dealing with the icon.

set(ICON_PATH "${PROJECT_SOURCE_DIR}/cmake/packaging/openmeeg.ico")
set(CPACK_NSIS_MUI_ICON ${ICON_PATH})
set(CPACK_NSIS_MUI_UNIICON ${ICON_PATH})
# set(CPACK_NSIS_INSTALLED_ICON_NAME bin\\\\om_assemble.exe)
set(CPACK_NSIS_DELETE_ICONS_EXTRA " Delete '\$SMPROGRAMS\\\\$MUI_TEMP\\\\*.*' ")

#   Add OpenMEEG to the PATH and shortcut in the Startup menu and/or on the desktop.

set(CPACK_NSIS_MODIFY_PATH ON)
set(CPACK_PACKAGE_EXECUTABLES "om_assemble" "OpenMEEG (Ignore)")

#   Add a link to the application website in the Startup menu.

set(CPACK_NSIS_MENU_LINKS
    "OpenMEEG/doc/LICENSE.txt"        "License"
    "OpenMEEG/README.rst"             "README"
    "http://openmeeg.github.io"       "OpenMEEG homepage"
)

#   Run OpenMEEG after installation
#set(CPACK_NSIS_MUI_FINISHPAGE_RUN "om_assemble.exe")

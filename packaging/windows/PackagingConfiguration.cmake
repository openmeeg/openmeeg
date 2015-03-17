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

set(DISTRIB windows)
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
	SET(PACKAGE_ARCH x64)
else()
    SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
	SET(CPACK_NSIS_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
	SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${CPACK_PACKAGE_VERSION}")
	SET(PACKAGE_ARCH x86)
endif()

#   There is a bug in NSIS that does not handle full unix paths properly. Make
#   sure there is at least one set of four (4) backlasshes.

set(CPACK_NSIS_DISPLAY_NAME "OpenMEEG")
set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\openmeeg.gforge.inria.fr")
set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\openmeeg.gforge.inria.fr")
set(CPACK_NSIS_CONTACT ${CPACK_PACKAGE_CONTACT})

#   Dealing with the icon.

set(ICON_PATH "${PROJECT_SOURCE_DIR}/medInria/app/medInria/resources/medInria.ico")
set(CPACK_NSIS_MUI_ICON ${ICON_PATH})
set(CPACK_NSIS_MUI_UNIICON ${ICON_PATH})
set(CPACK_NSIS_INSTALLED_ICON_NAME bin\\\\medInria.exe)
set(CPACK_NSIS_DELETE_ICONS_EXTRA " Delete '\$SMPROGRAMS\\\\$MUI_TEMP\\\\*.*' ")

#   Add OpenMEEG to the PATH and shortcut in the Startup menu and/or on the desktop.

set(CPACK_NSIS_MODIFY_PATH ON)
set(CPACK_PACKAGE_EXECUTABLES "om_assemble" "OpenMEEG (Ignore)")
#set(CPACK_CREATE_DESKTOP_LINKS "medInria")

#   Add a link to the application website in the Startup menu.

set(CPACK_NSIS_MENU_LINKS
    "OpenMEEG/doc/LICENSE.txt"        "License"
    "OpenMEEG/README.rst"             "README"
    "http://openmeeg.gforge.inria.fr" "OpenMEEG homepage"
)

#   Run OpenMEEG after installation
#set(CPACK_NSIS_MUI_FINISHPAGE_RUN "medInria.exe")

include(UseWin32dlls)

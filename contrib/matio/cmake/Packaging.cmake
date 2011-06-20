#-----------------------------------------------
# packaging
#-----------------------------------------------

OPTION(ENABLE_PACKAGING "Enable Packaging" ON)

IF(CMAKE_C_COMPILER MATCHES gcc)
    EXEC_PROGRAM(${CMAKE_C_COMPILER}
        ARGS -dumpversion
        OUTPUT_VARIABLE PACKAGE_COMPILER)
    SET(PACKAGE_COMPILER gcc-${PACKAGE_COMPILER})
ELSE()
    SET(PACKAGE_COMPILER ${CMAKE_CXX_COMPILER})
ENDIF()

IF (UNIX AND NOT APPLE) # LINUX
    OPTION(BUILD_RPM "Enable RPM Packaging" OFF)
ENDIF()

IF(ENABLE_PACKAGING OR BUILD_RPM)

    INCLUDE(InstallRequiredSystemLibraries)

    CONFIGURE_FILE(${MATIO_SOURCE_DIR}/README  README.txt  COPYONLY)
    CONFIGURE_FILE(${MATIO_SOURCE_DIR}/COPYING COPYING.txt COPYONLY)

    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "matio library for matlab IOs")
    SET(CPACK_PACKAGE_VENDOR "Christopher Hulbert")
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${MATIO_BINARY_DIR}/README.txt")
    SET(CPACK_RESOURCE_FILE_LICENSE "${MATIO_BINARY_DIR}/COPYING.txt")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "matio")
    SET(CPACK_PACKAGE_CONTACT "cch@isl-inc.com")

    IF(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64)
        SET(CPACK_RPM_PACKAGE_ARCHITECTURE x86_64)
    ELSE()
        SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386)
        SET(CPACK_RPM_PACKAGE_ARCHITECTURE i386)
    ENDIF()

    SET(PACKAGE_NAME "matio-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")
    IF(UNIX)
        IF (APPLE)
            SET(PACKAGE_NAME ${PACKAGE_NAME}-MacOSX-Intel)
            IF(BUILD_UNIVERSAL)
                SET(PACKAGE_NAME ${PACKAGE_NAME}-Universal)
            ENDIF()
        ELSE()
            SET(PACKAGE_NAME ${PACKAGE_NAME}-Linux.${CPACK_RPM_PACKAGE_ARCHITECTURE})
        ENDIF()
    ELSE()
        SET(PACKAGE_NAME ${PACKAGE_NAME}-win32-x86)
    ENDIF()

    SET(PACKAGE_NAME ${PACKAGE_NAME}-${PACKAGE_COMPILER})

    IF (USE_OMP)
        SET(PACKAGE_NAME ${PACKAGE_NAME}-OpenMP)
    ENDIF()

    IF(BUILD_SHARED)
        IF(PYTHON_WRAP)
            SET(PACKAGE_NAME ${PACKAGE_NAME}-python)
        ENDIF()
        SET(PACKAGE_NAME ${PACKAGE_NAME}-shared)
    ELSE()
        SET(PACKAGE_NAME ${PACKAGE_NAME}-static)
    ENDIF()

    SET(CPACK_PACKAGE_FILE_NAME ${PACKAGE_NAME})

    IF (WIN32)
        # There is a bug in NSIS that does not handle full unix paths properly. Make
        # sure there is at least one set of four (4) backlasshes.
        SET(CPACK_NSIS_DISPLAY_NAME "matio library for matlab IOs")
        SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\sourceforge.net/projects/matio/support")
        SET(CPACK_NSIS_URL_INFO_ABOUT "https:\\\\\\\\gforge.inria.fr/projects/openmeeg/")
        SET(CPACK_NSIS_CONTACT "cch@isl-inc.com")
        SET(CPACK_NSIS_MODIFY_PATH ON)
    ENDIF()

    SET(CPACK_SOURCE_STRIP_FILES "")

    IF(UNIX AND NOT APPLE)
        SET(CPACK_GENERATOR "TGZ")
    ENDIF()

    IF(APPLE)
        SET(CPACK_GENERATOR "PackageMaker;TGZ")
    ENDIF()

    INCLUDE(CPack)

    IF(UNIX AND BUILD_RPM) # linux
        SET(CPACK_GENERATOR "${CPACK_GENERATOR};RPM")
        IF (CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION LESS 8)
            INCLUDE(UseRPMTools)
            IF (RPMTools_FOUND)
                RPMTools_ADD_RPM_TARGETS(${PROJECT_NAME} "${PROJECT_SOURCE_DIR}/packaging/${PROJECT_NAME}.spec.in")
            ENDIF()
        ELSE()
            SET(CPACK_RPM_PACKAGE_LICENSE "LGPL")
            SET(CPACK_RPM_PACKAGE_DESCRIPTION  "libmatio is an open-source library for reading/writing Matlab MAT files.  This
library is designed for use by programs/libraries that do not have access or
do not want to rely on Matlab's libmat shared library.")
            SET(CPACK_RPM_PACKAGE_GROUP "Libraries")
        ENDIF()
    ENDIF()

ENDIF()

#IF (ENABLE_PACKAGING AND WIN32)
#    INCLUDE(UseWin32dlls)
#ENDIF()

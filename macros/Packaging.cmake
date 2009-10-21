#-----------------------------------------------
# packaging
#-----------------------------------------------

OPTION(ENABLE_PACKAGING "Enable Packaging" ON)

IF (UNIX AND NOT APPLE) # LINUX
    OPTION(BUILD_RPM "Enable RPM Packaging" OFF)
ENDIF()

IF(ENABLE_PACKAGING OR BUILD_RPM)

    INCLUDE(InstallRequiredSystemLibraries)

    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMEEG Project")
    SET(CPACK_PACKAGE_VENDOR "INRIA - Odyssee ENPC/INRIA/Ens Ulm")
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${OpenMEEG_SOURCE_DIR}/README.txt")
    SET(CPACK_RESOURCE_FILE_LICENSE "${OpenMEEG_SOURCE_DIR}/LICENSE.txt")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMEEG")
    SET(CPACK_PACKAGE_CONTACT "openmeeg-info_at_lists.gforge.inria.fr")

    IF(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64)
        SET(CPACK_RPM_PACKAGE_ARCHITECTURE x86_64)
    ELSEIF()
        SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386)
        SET(CPACK_RPM_PACKAGE_ARCHITECTURE i386)
    ENDIF()

    IF (APPLE)
        IF (USE_OMP)
            SET(CPACK_PACKAGE_FILE_NAME
                "OpenMEEG-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-MacOSX-Intel-gcc42-OpenMP")
        ELSE()
            SET(CPACK_PACKAGE_FILE_NAME
                "OpenMEEG-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-MacOSX-Intel")
        ENDIF()
    ENDIF()

    IF (UNIX AND NOT APPLE)
        SET(CPACK_PACKAGE_FILE_NAME
            "OpenMEEG-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-Linux.${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
    ENDIF()
    IF (WIN32)
        SET(CPACK_PACKAGE_FILE_NAME
            "OpenMEEG-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-win32-x86")
    ENDIF()

    IF (WIN32)
        # There is a bug in NSIS that does not handle full unix paths properly. Make
        # sure there is at least one set of four (4) backlasshes.
        SET(CPACK_NSIS_DISPLAY_NAME "OpenMEEG Project")
        SET(CPACK_NSIS_HELP_LINK "https:\\\\\\\\gforge.inria.fr/projects/openmeeg/")
        SET(CPACK_NSIS_URL_INFO_ABOUT "https:\\\\\\\\gforge.inria.fr/projects/openmeeg/")
        SET(CPACK_NSIS_CONTACT "openmeeg-info@lists.gforge.inria.fr")
        SET(CPACK_NSIS_MODIFY_PATH ON)
    ENDIF()

    SET(CPACK_SOURCE_STRIP_FILES "")

    INCLUDE(CPack)

    IF(UNIX AND BUILD_RPM) # linux
        IF (CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION LESS 8)
            INCLUDE(UseRPMTools)
            IF (RPMTools_FOUND)
                RPMTools_ADD_RPM_TARGETS(${PROJECT_NAME} "packaging/${PROJECT_NAME}.spec.in")
            ENDIF()
        ELSE()
            SET(CPACK_RPM_USER_BINARY_SPECFILE "packaging/${PROJECT_NAME}.spec.in")
        ENDIF()
    ENDIF()
ENDIF()

IF (ENABLE_PACKAGING AND WIN32)
    INCLUDE(macros/UseWin32dlls.cmake)
ENDIF()

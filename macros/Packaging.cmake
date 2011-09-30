#-----------------------------------------------
# packaging
#-----------------------------------------------

OPTION(ENABLE_PACKAGING "Enable Packaging" ON)

IF(CMAKE_C_COMPILER MATCHES gcc)
    EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
        ARGS -dumpversion
        OUTPUT_VARIABLE PACKAGE_COMPILER)
    SET(PACKAGE_COMPILER gcc-${PACKAGE_COMPILER})
ELSE()
    SET(PACKAGE_COMPILER ${CMAKE_CXX_COMPILER})
ENDIF()

IF (UNIX AND NOT APPLE) # LINUX
    OPTION(BUILD_RPM "Enable RPM Packaging" OFF)
ENDIF()

# Install README
INSTALL(FILES LICENSE.txt README.txt
    DESTINATION share/doc/OpenMEEG)

IF(ENABLE_PACKAGING OR BUILD_RPM)

    INCLUDE(InstallRequiredSystemLibraries)

    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMEEG Project")
    SET(CPACK_PACKAGE_VENDOR "INRIA-ENPC")
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${OpenMEEG_SOURCE_DIR}/README.txt")
    SET(CPACK_RESOURCE_FILE_LICENSE "${OpenMEEG_SOURCE_DIR}/LICENSE.txt")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMEEG")
    SET(CPACK_PACKAGE_CONTACT "openmeeg-info_at_lists.gforge.inria.fr")

    IF(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
        SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64)
        SET(CPACK_RPM_PACKAGE_ARCHITECTURE x86_64)
        SET(NBITS 64)
    ELSE()
        SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386)
        SET(CPACK_RPM_PACKAGE_ARCHITECTURE i386)
        SET(NBITS 32)
    ENDIF()

    SET(PACKAGE_NAME "OpenMEEG-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")
    IF(UNIX)
        IF (APPLE)
            SET(PACKAGE_NAME ${PACKAGE_NAME}-MacOSX-Intel)
            IF(BUILD_UNIVERSAL)
                SET(PACKAGE_NAME ${PACKAGE_NAME}-Universal)
            ENDIF()
        ELSE()
            SET(PACKAGE_NAME ${PACKAGE_NAME}-Linux${NBITS}.${CPACK_DEBIAN_PACKAGE_ARCHITECTURE})
        ENDIF()
    ELSE()
	if (CMAKE_CL_64)
	    set (CPACK_SYSTEM_NAME win64)
            SET(PACKAGE_NAME ${PACKAGE_NAME}-win64-x86_64)
	else ()
	    set (CPACK_SYSTEM_NAME win32)
            SET(PACKAGE_NAME ${PACKAGE_NAME}-win32-x86)
	endif ()
    ENDIF()

    SET(PACKAGE_NAME ${PACKAGE_NAME}-${PACKAGE_COMPILER})

    IF (USE_OMP)
        SET(PACKAGE_NAME ${PACKAGE_NAME}-OpenMP)
    ENDIF()

    IF(BUILD_SHARED_LIBS)
        IF(ENABLE_PYTHON)
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
        SET(CPACK_NSIS_DISPLAY_NAME "OpenMEEG")
        SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\openmeeg.gforge.inria.fr")
        SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\openmeeg.gforge.inria.fr")
        SET(CPACK_NSIS_CONTACT "openmeeg-info@lists.gforge.inria.fr")
        SET(CPACK_NSIS_MODIFY_PATH ON)
        SET(CPACK_PACKAGE_EXECUTABLES "om_assemble" "OpenMEEG (Ignore)")
        SET(CPACK_NSIS_MENU_LINKS
            "doc/LICENSE.txt" "README"
            "http://openmeeg.gforge.inria.fr" "OpenMEEG homepage"
        )

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
        IF (CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION LESS 10)
            INCLUDE(UseRPMTools)
            IF (RPMTools_FOUND)
                RPMTools_ADD_RPM_TARGETS(${PROJECT_NAME} "${PROJECT_SOURCE_DIR}/packaging/${PROJECT_NAME}.spec.in")
            ENDIF()
        ELSE()
            SET(CPACK_RPM_PACKAGE_LICENSE "CeCILL-B")
            SET(CPACK_RPM_PACKAGE_DESCRIPTION  "OpenMEEG is a package for forward/inverse problems of EEG/MEG. The forward problem uses the symmetric Boundary Element Method. The inverse problem uses a distributed approach (L2, L1 regularization). Developped within Odyssee (INRIA-ENPC-ENS).")
            SET(CPACK_RPM_PACKAGE_GROUP "Applications/Medical")
        ENDIF()
    ENDIF()

ENDIF()

IF (ENABLE_PACKAGING AND WIN32)
    INCLUDE(UseWin32dlls)
ENDIF()

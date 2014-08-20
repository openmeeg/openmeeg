#-----------------------------------------------
# packaging
#-----------------------------------------------

option(ENABLE_PACKAGING "Enable Packaging" ON)

if (CMAKE_C_COMPILER MATCHES gcc)
    exec_program(${CMAKE_CXX_COMPILER}
        ARGS -dumpversion
        OUTPUT_VARIABLE PACKAGE_COMPILER)
    set(PACKAGE_COMPILER gcc-${PACKAGE_COMPILER})
else()
    set(PACKAGE_COMPILER ${CMAKE_CXX_COMPILER})
endif()

if (UNIX AND NOT APPLE) # LINUX
    option(BUILD_RPM "Enable RPM Packaging" OFF)
endif()

# Install README
install(FILES ${PROJECT_SOURCE_DIR}/OpenMEEG/LICENSE.txt  ${PROJECT_SOURCE_DIR}/OpenMEEG/README.rst
        DESTINATION ${INSTALL_DATA_DIR}/doc/OpenMEEG)

if (NOT PACKAGE_NAME)
    if (ENABLE_PACKAGING OR BUILD_RPM)

        include(InstallRequiredSystemLibraries)

        set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMEEG Project")
        set(CPACK_PACKAGE_VENDOR "INRIA-ENPC")
        set(CPACK_PACKAGE_DESCRIPTION_FILE "${OpenMEEG_SOURCE_DIR}/README.rst")
        set(CPACK_RESOURCE_FILE_LICENSE "${OpenMEEG_SOURCE_DIR}/LICENSE.txt")
        set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMEEG")
        set(CPACK_PACKAGE_CONTACT "openmeeg-info_at_lists.gforge.inria.fr")

        if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
            set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64)
            set(CPACK_RPM_PACKAGE_ARCHITECTURE x86_64)
            set(NBITS 64)
        else()
            set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386)
            set(CPACK_RPM_PACKAGE_ARCHITECTURE i386)
            set(NBITS 32)
        endif()

        set(PACKAGE_NAME "OpenMEEG-${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")
        if (UNIX)
            if (APPLE)
                set(PACKAGE_NAME ${PACKAGE_NAME}-MacOSX-Intel)
                if (BUILD_UNIVERSAL)
                    set(PACKAGE_NAME ${PACKAGE_NAME}-Universal)
                endif()
            else()
                set(PACKAGE_NAME ${PACKAGE_NAME}-Linux${NBITS}.${CPACK_DEBIAN_PACKAGE_ARCHITECTURE})
            endif()
        else()
            if (CMAKE_CL_64)
                set(CPACK_SYSTEM_NAME win64)
                set(PACKAGE_NAME ${PACKAGE_NAME}-win64-x86_64)
            else()
                set(CPACK_SYSTEM_NAME win32)
                set(PACKAGE_NAME ${PACKAGE_NAME}-win32-x86)
            endif()
        endif()

        set(PACKAGE_NAME ${PACKAGE_NAME}-${PACKAGE_COMPILER})

        if (USE_OMP)
            set(PACKAGE_NAME ${PACKAGE_NAME}-OpenMP)
        endif()

        if (BUILD_SHARED_LIBS)
            if (ENABLE_PYTHON)
                set(PACKAGE_NAME ${PACKAGE_NAME}-python)
            endif()
            set(PACKAGE_NAME ${PACKAGE_NAME}-shared)
        else()
            set(PACKAGE_NAME ${PACKAGE_NAME}-static)
        endif()

        set(CPACK_PACKAGE_FILE_NAME ${PACKAGE_NAME})

        if (WIN32)
            # There is a bug in NSIS that does not handle full unix paths properly. Make
            # sure there is at least one set of four (4) backlasshes.
            set(CPACK_NSIS_DISPLAY_NAME "OpenMEEG")
            set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\openmeeg.gforge.inria.fr")
            set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\openmeeg.gforge.inria.fr")
            set(CPACK_NSIS_CONTACT "openmeeg-info@lists.gforge.inria.fr")
            set(CPACK_NSIS_MODIFY_PATH ON)
            set(CPACK_PACKAGE_EXECUTABLES "om_assemble" "OpenMEEG (Ignore)")
            set(CPACK_NSIS_MENU_LINKS
                "doc/LICENSE.txt" "README.rst"
                "http://openmeeg.gforge.inria.fr" "OpenMEEG homepage"
            )

        endif()

        set(CPACK_SOURCE_STRIP_FILES "")

        if (UNIX AND NOT APPLE)
            set(CPACK_GENERATOR "TGZ")
        endif()

        if (APPLE)
            set(CPACK_GENERATOR "PackageMaker;TGZ")
        endif()

        include(CPack)

        if (UNIX AND BUILD_RPM) # linux
            set(CPACK_GENERATOR "${CPACK_GENERATOR};RPM")
            if (CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION LESS 10)
                include(UseRPMTools)
                if (RPMTools_FOUND)
                    RPMTools_ADD_RPM_TARGETS(${PROJECT_NAME} "${PROJECT_SOURCE_DIR}/packaging/${PROJECT_NAME}.spec.in")
                endif()
            else()
                set(CPACK_RPM_PACKAGE_LICENSE "CeCILL-B")
                set(CPACK_RPM_PACKAGE_DESCRIPTION  "OpenMEEG is a package for forward/inverse problems of EEG/MEG. The forward problem uses the symmetric Boundary Element Method. The inverse problem uses a distributed approach (L2, L1 regularization). Developped within Odyssee (INRIA-ENPC-ENS).")
                set(CPACK_RPM_PACKAGE_GROUP "Applications/Medical")
            endif()
        endif()
    endif()
endif()

if (ENABLE_PACKAGING AND WIN32)
    include(UseWin32dlls)
endif()

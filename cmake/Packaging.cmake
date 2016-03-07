#-----------------------------------------------
# packaging
#-----------------------------------------------

option(ENABLE_PACKAGING "Enable Packaging" ON)

if (ENABLE_PACKAGING)

    set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/packaging ${CMAKE_MODULE_PATH})

    set(CPACK_PACKAGE_NAME "OpenMEEG")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMEEG Project")
    set(CPACK_PACKAGE_VENDOR "INRIA-ENPC")
    set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
    #set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/OpenMEEG/README.rst")
    #set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/OpenMEEG/LICENSE.txt")
    set(CPACK_PACKAGE_CONTACT "openmeeg-info@lists.gforge.inria.fr")

    set(CPACK_SET_DESTDIR true)
    set(CPACK_INSTALL_PREFIX "Packaging")
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMEEG")
    set(CPACK_SOURCE_STRIP_FILES "")

    if (UNIX)
        set(SYSTEMDIR linux)
        if (APPLE)
            set(SYSTEMDIR apple)
        endif()
    else()
        set(SYSTEMDIR windows)
    endif()
    include(${SYSTEMDIR}/PackagingConfiguration)

    if (USE_OMP)
        set(PACKAGE_OPTIONS OpenMP)
    endif()

    if (BUILD_SHARED_LIBS)
        if (ENABLE_PYTHON)
            set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-python)
        endif()
        set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-shared)
    else()
        set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-static)
    endif()

    if (PACKAGE_OPTIONS)
        set(PACKAGE_OPTIONS "-${PACKAGE_OPTIONS}")
    endif()

    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${PACKAGE_ARCH}${PACKAGE_OPTIONS}")

    include(InstallRequiredSystemLibraries)
    include(CPack)

endif()

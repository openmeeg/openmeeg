#-----------------------------------------------
# packaging
#-----------------------------------------------

option(ENABLE_PACKAGING "Enable Packaging" ON)


if (ENABLE_PACKAGING)
    set(CPACK_GENERATOR "TGZ")

    # set variables
    set(CPACK_PACKAGE_NAME "OpenMEEG")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMEEG Project")
    set(CPACK_PACKAGE_DESCRIPTION "A C++ package for low-frequency bio-electromagnetism solving forward problems in the field of EEG and MEG.")
    set(CPACK_PACKAGE_VENDOR "INRIA-ENPC")
    set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})

    set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.rst")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")

    set(CPACK_PACKAGE_CONTACT "openmeeg-info@lists.gforge.inria.fr")
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMEEG")
    set(CPACK_SOURCE_STRIP_FILES "")

    set(PACKAGE_ARCH_SHORT "Linux")
    if(WIN32)
        SET(PACKAGE_ARCH_SHORT "Win64")
        if(NOT CMAKE_CL_64)
            set(PACKAGE_ARCH_SHORT "Win32")
        endif()
    elseif(APPLE)
        set(PACKAGE_ARCH_SHORT "MacOSX")
    endif()

    # set(PACKAGE_OPTIONS ${BLASLAPACK_IMPLEMENTATION})
    # set(PACKAGE_OPTIONS)

    # if (USE_OMP)
    #     set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-OpenMP)
    # endif()

    # if (USE_VTK)
    #     set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-vtk)
    # endif()

    # if (USE_CGAL)
    #     set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-cgal)
    # endif()

    # if (ENABLE_PYTHON)
    #     set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-python)
    # endif()

    # set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${PACKAGE_ARCH_SHORT}-${PACKAGE_OPTIONS}")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${PACKAGE_ARCH_SHORT}")

    # Following https://blog.quickmediasolutions.com/2017/11/24/using-windeployqt-with-cpack.html
    set(CMAKE_INSTALL_UCRT_LIBRARIES TRUE)
    set(CMAKE_INSTALL_OPENMP_LIBRARIES TRUE)
    include(InstallRequiredSystemLibraries)
    include(CPack)

endif()

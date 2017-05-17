#-----------------------------------------------
# packaging
#-----------------------------------------------

option(ENABLE_PACKAGING "Enable Packaging" ON)

if (ENABLE_PACKAGING)
    # set(CPACK_RPM_PACKAGE_DEBUG TRUE)

    set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/OpenMEEG/cmake/packaging ${CMAKE_MODULE_PATH})

    # load variables
    include(PackagingOpenMEEG)

    set(CPACK_GENERATOR "TGZ")

    # if we want to generate all the sub-project packages:
    #set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR}/matio/build;matio;ALL;/" "${CMAKE_CURRENT_BINARY_DIR}/OpenMEEG/build;OpenMEEG;ALL;/")

    # call cpack for OpenMEEG
    set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR}/OpenMEEG/build;OpenMEEG;ALL;/")

    set(PACKAGE_OPTIONS ${BLASLAPACK_IMPLEMENTATION})

    if (BUILD_SHARED_LIBS)
        set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-shared)
    else()
        set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-static)
    endif()

    if (USE_OMP)
        set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-OpenMP)
    endif()

    if (USE_VTK)
        set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-vtk)
    endif()

    if (USE_CGAL)
        set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-cgal)
    endif()

    if (ENABLE_PYTHON)
        set(PACKAGE_OPTIONS ${PACKAGE_OPTIONS}-python)
    endif()

    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${PACKAGE_ARCH_SHORT}-${PACKAGE_OPTIONS}")

    include(InstallRequiredSystemLibraries)
    include(CPack)

endif()

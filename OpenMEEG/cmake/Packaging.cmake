#-----------------------------------------------
# packaging
#-----------------------------------------------

option(ENABLE_PACKAGING "Enable Packaging" OFF)

if (ENABLE_PACKAGING)

    set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/packaging ${CMAKE_MODULE_PATH})

    set(PACKAGE_COMPILER ${CMAKE_CXX_COMPILER})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        exec_program(${CMAKE_CXX_COMPILER}
            ARGS -dumpversion
            OUTPUT_VARIABLE PACKAGE_COMPILER)
        set(PACKAGE_COMPILER Clang-${PACKAGE_COMPILER})
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        exec_program(${CMAKE_CXX_COMPILER}
            ARGS -dumpversion
            OUTPUT_VARIABLE PACKAGE_COMPILER)
        set(PACKAGE_COMPILER gcc-${PACKAGE_COMPILER})
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(PACKAGE_COMPILER "icc")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        set(PACKAGE_COMPILER "MSVC")
    endif()

    if (UNIX AND NOT APPLE) # LINUX
        option(BUILD_LINUX_PACKAGE "Enable RPM or Debian Packaging" ON)
    endif()

    if (ENABLE_PACKAGING OR BUILD_LINUX_PACKAGE)

        include(InstallRequiredSystemLibraries)

        # load variables
        include(PackagingOpenMEEG)

        if (UNIX)
            set(SYSTEMDIR linux)
            if (APPLE)
                set(SYSTEMDIR apple)
            endif()
        else()
            set(CPACK_SET_DESTDIR false)
            set(CPACK_INSTALL_PREFIX "")
            set(SYSTEMDIR windows)
        endif()

        include(${SYSTEMDIR}/PackagingConfiguration)

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

        include(CPack)

    endif()
endif()

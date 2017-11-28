#   Download VTK from Kitware web site and install it.

if (WIN32)

    message(STATUS "Downloading Kitware VTK installer, this may take a while...")

    set(VTK_INSTALLER_ARCHIVE vtkpython-7.1.1-Windows-64bit.exe)
    file(DOWNLOAD "http://www.vtk.org/files/release/7.1/${VTK_INSTALLER_ARCHIVE}" ${CMAKE_BINARY_DIR}/${VTK_INSTALLER_ARCHIVE}
         STATUS result)
    list(GET result 0 error_code)
    if (NOT ${error_code} STREQUAL "0")
        message(FATAL_ERROR "Could not download MKL install script. If no network connexion please provide MKL_DIR or environment {MKLDIR}")
    endif()

    #   Install VTK

    message(STATUS "Installing VTK, this may take a while...")

    execute_process(COMMAND cmd /c start /wait ${CMAKE_BINARY_DIR}/${VTK_INSTALLER_ARCHIVE} /S
                    OUTPUT_FILE ${CMAKE_BINARY_DIR}/install-vtk.out
                    ERROR_FILE ${CMAKE_BINARY_DIR}/install-vtk.err
                    RESULT_VARIABLE vtk_install_result)

    message("[[install-vtk.out]]")
    message("${OUTPUT_FILE}")

    message("[[install-vtk.err]]")
    message("${ERROR_FILE}")

    if (NOT ${vtk_install_result} STREQUAL "0")
        message(FATAL_ERROR "Could not install VTK: please look at files install-vtk.{out,err} or provide VTK_DIR.")
    endif()

    message("Looking for VTK...")

    set(VTK_COMPONENTS vtkIOXML vtkIOLegacy)
    find_package(VTK ${REQUIRED} COMPONENTS ${VTK_COMPONENTS)

    if (NOT VTK_FOUND)
        message(FATAL_ERROR "VTK seems to be incorrectly installed.")
    endif()
    message("VTK_DIR: ${VTK_DIR}")
endif()

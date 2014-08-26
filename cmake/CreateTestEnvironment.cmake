#   Create the proper environment variables for platforms that need those.

if (WIN32 OR APPLE)
    set(lib_separator ":")
    set(dll_var "DYLD_LIBRARY_PATH")
    if (WIN32)
        set(lib_separator ";")
        set(dll_var "PATH")
    endif()

    #   Python dir

    if (ENABLE_PYTHON)
        get_property(PYTHON_OPENMEEG_MODULE TARGET "_openmeeg" PROPERTY LOCATION)
        get_filename_component(PYTHON_OPENMEEG_MODULE_DIR ${PYTHON_OPENMEEG_MODULE} DIRECTORY)
    endif()

    #   Hdf5 dir

    set(TestConfigFile "${CMAKE_BINARY_DIR}/TestConfig.cmake")
    get_property(HDF5_LIB TARGET hdf5 PROPERTY LOCATION)
    get_filename_component(HDF5_ROOT_DIR ${HDF5_DIR} DIRECTORY)

    set(DLL_DIRS "${ZLIB_ROOT}/bin" "${HDF5_ROOT_DIR}/bin" "${matio_ROOT_DIR}/bin" "${LAPACK_DLL_DIR}"
                 "${PYTHON_OPENMEEG_MODULE_DIR}" "${VTK_LIBRARY_DIRS}" "${CGAL_LIBRARY_DIRS}" "${NIFTI_DIR}")
    foreach (dir ${DLL_DIRS})
        set(LIBRARY_PATHS "${dir}${lib_separator}${LIBRARY_PATHS}")
    endforeach()
    file(WRITE ${TestConfigFile} "set(ENV{${dll_var}} \"${LIBRARY_PATHS}\$ENV{${dll_var}}\")\n")
    set_directory_properties(PROPERTIES TEST_INCLUDE_FILE "${TestConfigFile}")
endif()


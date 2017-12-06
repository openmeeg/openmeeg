#   Create the proper environment variables for platforms that need those.

if (WIN32 OR APPLE)
    set(lib_separator ":")
    set(dll_var "DYLD_LIBRARY_PATH")
    set(subdir "lib")
    if (WIN32)
        set(lib_separator ";")
        set(dll_var "PATH")
        set(subdir "bin")
    endif()

    #   Python dir

    if (ENABLE_PYTHON)
        get_property(PYTHON_OPENMEEG_MODULE TARGET "_openmeeg" PROPERTY LOCATION)
        get_filename_component(PYTHON_OPENMEEG_MODULE_DIR ${PYTHON_OPENMEEG_MODULE} DIRECTORY)
        if (WIN32)
            get_filename_component(PYTHON_OPENMEEG_MODULE_DIR ${PYTHON_OPENMEEG_MODULE_DIR} DIRECTORY)
            set(PYTHON_OPENMEEG_MODULE_DIR "${PYTHON_OPENMEEG_MODULE_DIR}/${CMAKE_BUILD_TYPE}")
        endif()
    endif()

    #   Hdf5 dir

    set(TestConfigFile "${CMAKE_BINARY_DIR}/TestConfig.cmake")

    set(HDF5_ROOT_DIR ${HDF5_LIBRARY_DIRS})
    if (TARGET hdf5)
        get_property(HDF5_LIB TARGET hdf5 PROPERTY LOCATION)
        get_filename_component(HDF5_LIB ${HDF5_LIB} DIRECTORY)
        get_filename_component(HDF5_ROOT_DIR ${HDF5_LIB} DIRECTORY)
    endif()

    set(DLL_DIRS "${ZLIB_ROOT}/${subdir}" "${HDF5_ROOT_DIR}/${subdir}" "${matio_ROOT_DIR}/${subdir}"
                 "${CMAKE_MSVCIDE_RUN_PATH}" "${LAPACK_DLL_DIR}" ${MKL_LIBRARY_DIR}
                 "${PYTHON_OPENMEEG_MODULE_DIR}" "${VTK_LIBRARY_DIRS}" "${CGAL_LIBRARY_DIRS}" "${NIFTI_LIBRARY_DIR}")
    message("[[${DLL_DIRS}]]")
    message("[[${MKL_LIBRARY_DIR}]]")
    foreach(i ${MKL_LIBRARY_DIR})
        set(bb ${bb} ${i}/*)
    endforeach()
    file(GLOB_RECURSE aa ${bb})
    message("[[[${aa}]]]")
    foreach (dir ${DLL_DIRS})
        FILE(TO_NATIVE_PATH "${dir}" dir)
        string(REPLACE "\\" "\\\\" dir ${dir})
        set(LIBRARY_PATHS "${dir}${lib_separator}${LIBRARY_PATHS}")
    endforeach()
    file(WRITE ${TestConfigFile} "${CONFIG}set(ENV{${dll_var}} \"${LIBRARY_PATHS}\$ENV{${dll_var}}\")\n")
    file(APPEND ${TestConfigFile} "${CONFIG}set(ENV{PYTHONPATH} \"${LIBRARY_PATHS}\$ENV{${dll_var}}\")\n")
    set_directory_properties(PROPERTIES TEST_INCLUDE_FILE "${TestConfigFile}")
endif()


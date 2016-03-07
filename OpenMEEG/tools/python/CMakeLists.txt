find_package(PythonInterp)

if (PYTHON_EXECUTABLE)
    exec_program("${PYTHON_EXECUTABLE}"
        ARGS "-c \"from distutils.sysconfig import get_python_lib; print get_python_lib(1)\""
        OUTPUT_VARIABLE PYTHON_SITE_ARCH
        RETURN_VALUE PYTHON_SITEARCH_NOT_FOUND)
    string(REGEX REPLACE "^.*/(lib|lib64)(/python.*)$" "\\1\\2" PYTHON_REL_SITE_ARCH ${PYTHON_SITE_ARCH})

    install(FILES om_display.py
        DESTINATION ${PYTHON_REL_SITE_ARCH}
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
endif()

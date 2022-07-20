# Adapted from FindNumpy.cmake (MIT)

if (Pytest_FIND_REQUIRED)
    find_package(PythonInterp REQUIRED)
else()
    find_package(PythonInterp)
endif()

if (NOT PYTHONINTERP_FOUND)
    set(PYTEST_FOUND FALSE)
    return()
endif()

execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" "import pytest; print(pytest.__version__)"
    RESULT_VARIABLE _PYTEST_SEARCH_SUCCESS
    OUTPUT_VARIABLE PYTEST_VERSION
    ERROR_VARIABLE _PYTEST_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

if (NOT _PYTEST_SEARCH_SUCCESS MATCHES 0)
    if (Pytest_FIND_REQUIRED)
        message(FATAL_ERROR "pytest import failure:\n${_PYTEST_ERROR_VALUE}")
    endif()
    set(PYTEST_FOUND FALSE)
    return()
endif()

find_package_message(PYTEST "Found Pytest: version \"${PYTEST_VERSION}\"" "${PYTEST_VERSION}")

set(PYTEST_FOUND TRUE)

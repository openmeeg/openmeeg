#######################################################################
# SWIG Embedding
#######################################################################

#   For now force the use of Python 2.{5,6,7}, numpy is not yet well adapted to python3.
#   Reconsider this choice later.

set(Python_ADDITIONAL_VERSIONS 2.5 2.6 2.7)
find_package(PythonInterp)

if (PYTHON_EXECUTABLE)
    exec_program("${PYTHON_EXECUTABLE}"
      ARGS "-c \"from distutils.sysconfig import get_python_lib; print get_python_lib(1)\""
      OUTPUT_VARIABLE PYTHON_SITE_ARCH
      RETURN_VALUE PYTHON_SITEARCH_NOT_FOUND)
    string(REGEX REPLACE "^.*/(lib|lib64)(/python.*)$" "\\1\\2" PYTHON_REL_SITE_ARCH ${PYTHON_SITE_ARCH})
endif()

include(FindNumpy.cmake)

set(MODULE_NAME openmeeg)
set(SWIG_SRC_FILE_NAME ${MODULE_NAME}.i)
set(CMAKE_SWIG_FLAGS "-c++ -v" )

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
find_package(PythonLibs)

include_directories(${PYTHON_INCLUDE_PATH} ${PYTHON_NUMPY_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_SWIG_FLAGS "-v")

file(GLOB SWIG_SRC_FILE ${SWIG_SRC_FILE_NAME})

set_source_files_properties(${SWIG_SRC_FILE_NAME} ${SWIG_SRC_FILE} PROPERTIES CPLUSPLUS ON)

set (SWIG_MODULE_${MODULE_NAME}_EXTRA_DEPS ${SWIG_SRC_FILE})
if (BUILD_DOCUMENTATION)
    set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} -DDOCSTRINGS)
    add_custom_target(docstrings
                      COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/doxy2swig.py ${PROJECT_BINARY_DIR}/documentation/reference/xml/index.xml docstrings.i)
    add_dependencies(docstrings OpenMEEG-doxygen_outputs)
    add_dependencies(OpenMEEG-documentation docstrings)
    set (SWIG_MODULE_${MODULE_NAME}_EXTRA_DEPS ${SWIG_MODULE_${MODULE_NAME}_EXTRA_DEPS} ${PROJECT_BINARY_DIR}/Wrapping/src/docstrings.i)
endif()

SWIG_ADD_MODULE(${MODULE_NAME} python ${SWIG_SRC_FILE_NAME} "")

if (PYTHON_LIBRARIES)
    target_link_libraries("_${MODULE_NAME}" ${PYTHON_LIBRARIES} OpenMEEG)
endif()

install(TARGETS  "_${MODULE_NAME}"
        DESTINATION ${PYTHON_REL_SITE_ARCH}
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

install(FILES  ${OpenMEEG_BINARY_DIR}/Wrapping/src/${MODULE_NAME}.py
        DESTINATION ${PYTHON_REL_SITE_ARCH}
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

# Test Python wrapping

if (BUILD_TESTING)
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${OpenMEEG_SOURCE_DIR}/Wrapping/src/test_python.py ${OpenMEEG_BINARY_DIR}/Wrapping/src/)
    add_test(openmeeg_python
             ${PYTHON_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/Wrapping/src/test_python.py --path ${OpenMEEG_SOURCE_DIR}/data)
endif()

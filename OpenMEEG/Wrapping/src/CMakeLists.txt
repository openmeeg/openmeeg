#######################################################################
# SWIG Embedding
#######################################################################

find_package(PythonInterp REQUIRED)
message(STATUS "PYTHON VERSION: ${PYTHON_VERSION_STRING} ${PYTHON_EXECUTABLE}")

find_package(PythonLibs REQUIRED)
find_package(Numpy REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

set(MODULE_NAME openmeeg)
set(SWIG_SRC_FILE_NAME ${MODULE_NAME}.i)
set(CMAKE_SWIG_FLAGS -v -O)

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

message("-:- ${SWIG_SRC_FILE} :: ${SWIG_SRC_FILE_NAME} :_:")
SWIG_ADD_MODULE(${MODULE_NAME} python ${SWIG_SRC_FILE} "")
SWIG_LINK_LIBRARIES(${MODULE_NAME} ${PYTHON_LIBRARIES} OpenMEEG ${LAPACK_LIBRARIES} ${VTK_LIBRARIES})

execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c" "from distutils.sysconfig import get_python_lib; print get_python_lib(1)"
             OUTPUT_VARIABLE PYTHON_SITE_ARCH
             RESULT_VARIABLE PYTHON_SITEARCH_NOT_FOUND)
#extract correctly python location: e.g lib64/python2.7/site-packages
string(REGEX MATCH "lib.*packages" PYTHON_REL_SITE_ARCH ${PYTHON_SITE_ARCH})

install(TARGETS  "_${MODULE_NAME}"
        LIBRARY DESTINATION ${PYTHON_REL_SITE_ARCH}
        RUNTIME DESTINATION ${PYTHON_REL_SITE_ARCH}
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

install(FILES  ${OpenMEEG_BINARY_DIR}/Wrapping/src/${MODULE_NAME}.py
        DESTINATION ${PYTHON_REL_SITE_ARCH}
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

# Test Python wrapping
if (BUILD_TESTING)
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${OpenMEEG_SOURCE_DIR}/Wrapping/src/test_python.py ${OpenMEEG_BINARY_DIR}/Wrapping/src/)
    add_test(openmeeg_python ${PYTHON_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/Wrapping/src/test_python.py --path ${OpenMEEG_SOURCE_DIR}/data)
    if (USE_MKL AND BUILD_SHARED_LIBS AND NOT MKL_USE_sdl)
        # due to this bug in the MKL
        message(WARNING "due to a bug in the MKL when used from python: see https://software.intel.com/en-us/forums/intel-distribution-for-python/topic/628976 "
            "workaround: set environment variable:\nexport LD_PRELOAD=\"${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:"
            "${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so\"")
        set_tests_properties(openmeeg_python PROPERTIES ENVIRONMENT LD_PRELOAD=${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so)
    endif()
endif()

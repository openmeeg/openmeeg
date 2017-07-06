#######################################################################
# SWIG Embedding
#######################################################################

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Numpy REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

set(CMAKE_SWIG_FLAGS -v -O)
set(SWIG_SRC_FILE ${CMAKE_CURRENT_SOURCE_DIR}/openmeeg.i)
set_source_files_properties(${SWIG_SRC_FILE} PROPERTIES CPLUSPLUS ON)

set(SWIG_MODULE_OPENMEEG_EXTRA_DEPS ${SWIG_SRC_FILE})
if (BUILD_DOCUMENTATION)
    set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} -DDOCSTRINGS)
    add_custom_target(docstrings
                      COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/doxy2swig.py ${PROJECT_BINARY_DIR}/documentation/reference/xml/index.xml docstrings.i)
    add_dependencies(docstrings OpenMEEG-doxygen_outputs)
    add_dependencies(OpenMEEG-documentation docstrings)
    set (SWIG_MODULE_OPENMEEG_EXTRA_DEPS ${SWIG_MODULE_OPENMEEG_EXTRA_DEPS} ${PROJECT_BINARY_DIR}/wrapping/src/docstrings.i)
endif()

if ("${CMAKE_VERSION}" VERSION_LESS "3.8")
    swig_add_module(openmeeg python ${SWIG_SRC_FILE} "")
else()
    swig_add_library(openmeeg LANGUAGE python SOURCES ${SWIG_SRC_FILE})
endif()
swig_link_libraries(openmeeg ${PYTHON_LIBRARIES} OpenMEEG ${LAPACK_LIBRARIES} ${VTK_LIBRARIES})

execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c" "from distutils.sysconfig import get_python_lib; print(get_python_lib(1))"
                OUTPUT_VARIABLE PYTHON_SITE_ARCH
                RESULT_VARIABLE PYTHON_SITEARCH_NOT_FOUND)
#extract correctly python location: e.g lib64/python2.7/site-packages
string(REGEX MATCH "[Ll]ib.*packages" PYTHON_REL_SITE_ARCH ${PYTHON_SITE_ARCH})

install(TARGETS  "_openmeeg"
        LIBRARY DESTINATION ${PYTHON_REL_SITE_ARCH}
        RUNTIME DESTINATION ${PYTHON_REL_SITE_ARCH}
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

install(FILES  ${OpenMEEG_BINARY_DIR}/wrapping/src/openmeeg.py
        DESTINATION ${PYTHON_REL_SITE_ARCH}
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

# Test Python wrapping
if (BUILD_TESTING)
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${OpenMEEG_SOURCE_DIR}/wrapping/src/test_python.py ${OpenMEEG_BINARY_DIR}/wrapping/src/)
    add_test(openmeeg_python ${PYTHON_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/wrapping/src/test_python.py --path ${OpenMEEG_SOURCE_DIR}/data)
    if (USE_MKL AND BUILD_SHARED_LIBS AND NOT MKL_USE_sdl)
        # due to this bug in the MKL
        message(WARNING "Due to a bug in the MKL when used from python: see https://software.intel.com/en-us/forums/intel-distribution-for-python/topic/628976 "
            "workaround: set environment variable:\nexport LD_PRELOAD=\"${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:"
            "${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so\"")
        set_tests_properties(openmeeg_python PROPERTIES ENVIRONMENT LD_PRELOAD=${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so)
    endif()
endif()

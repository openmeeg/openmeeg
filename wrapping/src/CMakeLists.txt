#######################################################################
# SWIG Embedding
#######################################################################
#
# SWIG generation is driven by the following SWIG macros:
#
#   - swig_add_module
#   - swig_add_library
#   - swig_link_libraries
#
# those macros are driven by CMAKE_SWIG_FLAGS, and SWIG_MODULE_XXXXXXX_EXTRA_DEPS
# where XXXXXXX is the name of the target object in swig_add_library(XXXXXXX ...)

set(Python_ADDITIONAL_VERSIONS ${PYTHON_VERSION})
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Numpy REQUIRED)

# XXXX: this include_directories should be removed and transitive properties
#       of OpenMEEG::OpenMEEG and OpenMEEG::OpenMEEGMaths should be used instead
include_directories(${PYTHON_INCLUDE_PATH}
                    ${NUMPY_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_SOURCE_DIR}
                    ${PROJECT_SOURCE_DIR}/OpenMEEG/include/
                    ${PROJECT_SOURCE_DIR}/OpenMEEGMaths/include/
)

find_package(SWIG REQUIRED)
# load SWIG macros. The path to UseSWIG.cmake is provided by SWIG_USE_FILE which is set in FindSWIG
if (SWIG_USE_FILE STREQUAL "")
  message(FATAL_ERROR "unexpected: the variable SWIG_USE_FILE is the empty string. Did you run FindSWIG.cmake?")
else()
  include(${SWIG_USE_FILE})
endif()


if (BUILD_DOCUMENTATION)
    # In order to generate python docstrings we use doxy2swig (which is bundled in our codebase)
    # doxy2swig takes xml output from doxygen and generates a .i file that can be include in openmeeg.i
    set(DOXY2SWIG_EXECUTABLE ${PROJECT_SOURCE_DIR}/wrapping/src/doxy2swig.py)
    set(DOXY_INDEX_XML ${PROJECT_BINARY_DIR}/documentation/reference/xml/index.xml)
    add_custom_command(
        OUTPUT docstrings.i
        COMMAND ${PYTHON_EXECUTABLE} ${DOXY2SWIG_EXECUTABLE}
        ARGS    --function-signature --type-info --overloaded-functions --attribute-list --constructor-list
                ${DOXY_INDEX_XML} ${CMAKE_CURRENT_BINARY_DIR}/docstrings.i
        DEPENDS OpenMEEG-doxygen_outputs
        COMMENT "Generating Swig documentation" VERBATIM
    )

    list(APPEND CMAKE_SWIG_FLAGS "-DDOCSTRINGS")
    add_custom_target(docstrings DEPENDS docstrings.i)
    add_dependencies(OpenMEEG-documentation docstrings)
endif()

list(APPEND CMAKE_SWIG_FLAGS -v -O)
set_source_files_properties(openmeeg.i PROPERTIES CPLUSPLUS ON)
swig_add_library(openmeeg LANGUAGE python SOURCES openmeeg.i)
swig_link_libraries(openmeeg ${PYTHON_LIBRARIES} OpenMEEG::OpenMEEG)

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

set(PYTHON_TEST_FILES
  "test_matrix.py"
  "test_mesh.py"
  "test_sensors.py"
  "test_vector.py"
  )

set(PYTHON_TEST_FILES_REQUIRING_DATA
  "test_python.py"
  "test_python2.py"
  )

# Test Python wrapping
#
# XXXX: we copy all python scripts from source to binary so that openmeeg.py is in the path. There should be a manner of just adding
#       it to the python path for the testing environment and allow to execute the python scripts from source.
if (BUILD_TESTING)
    # Add tests requiring data
    foreach( T ${PYTHON_TEST_FILES_REQUIRING_DATA})
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${OpenMEEG_SOURCE_DIR}/wrapping/src/${T} ${OpenMEEG_BINARY_DIR}/wrapping/src/)
      add_test(openmeeg_python_${T} ${PYTHON_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/wrapping/src/${T} --path ${OpenMEEG_SOURCE_DIR}/data)
      if (USE_MKL AND BUILD_SHARED_LIBS AND NOT MKL_USE_sdl)
          # due to this bug in the MKL
          message(WARNING "Due to a bug in the MKL when used from python: see https://software.intel.com/en-us/forums/intel-distribution-for-python/topic/628976 "
              "workaround: set environment variable:\nexport LD_PRELOAD=\"${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:"
              "${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so\"")
          set_tests_properties(openmeeg_python_${T} PROPERTIES ENVIRONMENT LD_PRELOAD=${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so)
      endif()
    endforeach()

    # Add python tests not requiring data
    foreach( T ${PYTHON_TEST_FILES})
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${OpenMEEG_SOURCE_DIR}/wrapping/src/${T} ${OpenMEEG_BINARY_DIR}/wrapping/src/)
      add_test(openmeeg_python_${T} ${PYTHON_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/wrapping/src/${T})
      if (USE_MKL AND BUILD_SHARED_LIBS AND NOT MKL_USE_sdl)
          # due to this bug in the MKL
          message(WARNING "Due to a bug in the MKL when used from python: see https://software.intel.com/en-us/forums/intel-distribution-for-python/topic/628976 "
              "workaround: set environment variable:\nexport LD_PRELOAD=\"${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:"
              "${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so\"")
          set_tests_properties(openmeeg_python_${T} PROPERTIES ENVIRONMENT LD_PRELOAD=${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so)
      endif()
    endforeach()

    if(BUILD_REFERENCE_DOC)
        execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${OpenMEEG_SOURCE_DIR}/wrapping/src/test_doc.py ${OpenMEEG_BINARY_DIR}/wrapping/src/)
        add_test(openmeeg_python_doc ${PYTHON_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/wrapping/src/test_doc.py)
    endif(BUILD_REFERENCE_DOC)

endif()

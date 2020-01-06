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

find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

# XXXX: this include_directories should be removed and transitive properties
#       of OpenMEEG::OpenMEEG and OpenMEEG::OpenMEEGMaths should be used instead

include_directories(${Python3_INCLUDE_DIRS}
                    ${Python3_NumPy_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_SOURCE_DIR}
                    ${PROJECT_SOURCE_DIR}/OpenMEEG/include/
                    ${PROJECT_SOURCE_DIR}/OpenMEEGMaths/include/
)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13")
    cmake_policy(SET CMP0078 OLD) # Use old policy as our required cmake version (3.9) is older than 3.13 (For the NEW we need to chage the target name).
endif()

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.14")
    cmake_policy(SET CMP0086 NEW) # Use NEW policy as we do not use SWIG_MODULE_NAME
endif()

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
        COMMAND ${Python3_EXECUTABLE} ${DOXY2SWIG_EXECUTABLE}
        ARGS    --function-signature --type-info --overloaded-functions --attribute-list --constructor-list
                ${DOXY_INDEX_XML} ${CMAKE_CURRENT_BINARY_DIR}/docstrings.i
        DEPENDS OpenMEEG-doxygen_outputs
        COMMENT "Generating Swig documentation" VERBATIM
    )

    list(APPEND CMAKE_SWIG_FLAGS "-DDOCSTRINGS")
    add_custom_target(docstrings DEPENDS docstrings.i)
    add_dependencies(OpenMEEG-documentation docstrings)
    set(SWIG_MODULE_openmeeg_EXTRA_DEPS docstrings)
endif()

list(APPEND CMAKE_SWIG_FLAGS -v -O)
set(SWIG_MODULE_openmeeg_EXTRA_DEPS numpy.i make_geometry.py ${SWIG_MODULE_openmeeg_EXTRA_DEPS})
set_source_files_properties(openmeeg.i PROPERTIES CPLUSPLUS ON GENERATED_COMPILE_DEFINITIONS SWIG_PYTHON_SILENT_MEMLEAK)
swig_add_library(openmeeg LANGUAGE python SOURCES openmeeg.i)

if (APPLE)
  swig_link_libraries(openmeeg OpenMEEG::OpenMEEG)
else()
  swig_link_libraries(openmeeg ${Python3_LIBRARIES} OpenMEEG::OpenMEEG)
endif()

#extract correctly python location: e.g lib64/python2.7/site-packages

string(REGEX MATCH "[Ll]ib.*packages" PYTHON_REL_SITE_ARCH ${Python3_SITEARCH})

install(TARGETS  "_openmeeg"
        LIBRARY DESTINATION ${PYTHON_REL_SITE_ARCH}
        RUNTIME DESTINATION ${PYTHON_REL_SITE_ARCH}
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

if (APPLE)
  set_target_properties(_openmeeg PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

install(FILES  ${OpenMEEG_BINARY_DIR}/wrapping/src/openmeeg.py
        DESTINATION ${PYTHON_REL_SITE_ARCH}
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

# Test Python wrapping
#
# XXXX: we copy all python scripts from source to binary so that openmeeg.py is in the path. There should be a manner of just adding
#       it to the python path for the testing environment and allow to execute the python scripts from source.

if (BUILD_TESTING)

    set(PYTHON_TEST_FILES test_vector.py test_matrix.py test_sensors.py test_geom_import.py)
    set(PYTHON_TEST_FILES_REQUIRING_DATA test_python.py test_python2.py test_mesh.py test_make_geometry.py)

    # Add tests requiring data

    foreach (test ${PYTHON_TEST_FILES_REQUIRING_DATA})
        execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${OpenMEEG_SOURCE_DIR}/wrapping/src/${test} ${OpenMEEG_BINARY_DIR}/wrapping/src/)
        add_test(openmeeg_python_${test} ${Python3_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/wrapping/src/${test} --path ${OpenMEEG_SOURCE_DIR}/data)
        if (USE_MKL AND BUILD_SHARED_LIBS AND NOT MKL_USE_sdl)
            # due to this bug in the MKL
            message(WARNING "Due to a bug in the MKL when used from python: see https://software.intel.com/en-us/forums/intel-distribution-for-python/topic/628976 "
                "workaround: set environment variable:\nexport LD_PRELOAD=\"${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:"
                "${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so\"")
            set_tests_properties(openmeeg_python_${test} PROPERTIES ENVIRONMENT LD_PRELOAD=${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so)
        endif()
    endforeach()

    # Add python tests not requiring data

    foreach (test ${PYTHON_TEST_FILES})
        execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${OpenMEEG_SOURCE_DIR}/wrapping/src/${test} ${OpenMEEG_BINARY_DIR}/wrapping/src/)
        add_test(openmeeg_python_${test} ${Python3_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/wrapping/src/${test})
        if (USE_MKL AND BUILD_SHARED_LIBS AND NOT MKL_USE_sdl)
            # due to this bug in the MKL
            message(WARNING "Due to a bug in the MKL when used from python: see https://software.intel.com/en-us/forums/intel-distribution-for-python/topic/628976 "
                "workaround: set environment variable:\nexport LD_PRELOAD=\"${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:"
                "${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so\"")
            set_tests_properties(openmeeg_python_${test} PROPERTIES ENVIRONMENT LD_PRELOAD=${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so)
        endif()
    endforeach()

    if (BUILD_REFERENCE_DOC)
        execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${OpenMEEG_SOURCE_DIR}/wrapping/src/test_doc.py ${OpenMEEG_BINARY_DIR}/wrapping/src/)
        add_test(openmeeg_python_doc ${Python3_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/wrapping/src/test_doc.py)
    endif()

endif()

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

# CMake changed the swig target name in version 3.13

if (${CMAKE_VERSION} VERSION_GREATER "3.13")
    cmake_policy(SET CMP0078 NEW)
    set(SWIG_TARGET_NAME "openmeeg")
else()
    set(SWIG_TARGET_NAME "_openmeeg")
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

list(APPEND CMAKE_SWIG_FLAGS -v -O)
set(SWIG_MODULE_openmeeg_EXTRA_DEPS openmeeg/numpy.i ${SWIG_MODULE_openmeeg_EXTRA_DEPS})
set_source_files_properties(openmeeg/openmeeg.i PROPERTIES CPLUSPLUS ON GENERATED_COMPILE_DEFINITIONS SWIG_PYTHON_SILENT_MEMLEAK)
# This puts openmeeg.py in openmeeg/
swig_add_library(${SWIG_TARGET_NAME} LANGUAGE python OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/openmeeg SOURCES openmeeg/openmeeg.i)
# This puts _openmeeg.so in openmeeg/
set_target_properties(${SWIG_TARGET_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY openmeeg/)
if ( MSVC )
    set_target_properties(${SWIG_TARGET_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_DEBUG openmeeg/)
    set_target_properties(${SWIG_TARGET_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELEASE openmeeg/)
endif()

if (APPLE)
    swig_link_libraries(openmeeg OpenMEEG::OpenMEEG)
else()
    swig_link_libraries(openmeeg ${Python3_LIBRARIES} OpenMEEG::OpenMEEG)
endif()

set(PYTHON_SITE_ARCH ${Python3_SITEARCH})
if (${PYTHON_INSTALL_RELATIVE})
    string(REGEX MATCH "[Ll]ib.*packages" PYTHON_SITE_ARCH ${PYTHON_SITE_ARCH})
endif()
file(TO_CMAKE_PATH "${PYTHON_SITE_ARCH}" PYTHON_SITE_ARCH)
string(REGEX REPLACE "\n$" "" PYTHON_SITE_ARCH "${PYTHON_SITE_ARCH}")
message(STATUS "Python installation direcotry: ${PYTHON_SITE_ARCH}")

# Copy our .py files over, dealing with Windows wildcard issues
file(GLOB PYTHON_COPY_FILENAMES ${CMAKE_SOURCE_DIR}/wrapping/python/openmeeg/*.py)
add_custom_command(
    TARGET ${SWIG_TARGET_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${PYTHON_COPY_FILENAMES} ${CMAKE_CURRENT_BINARY_DIR}/openmeeg/)
add_custom_command(
    TARGET ${SWIG_TARGET_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/wrapping/python/MANIFEST.in
    ${CMAKE_SOURCE_DIR}/wrapping/python/setup.py
    ${CMAKE_SOURCE_DIR}/wrapping/python/README.rst
    ${CMAKE_CURRENT_BINARY_DIR})
add_custom_command(
    TARGET ${SWIG_TARGET_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/wrapping/python/openmeeg/tests
    ${CMAKE_CURRENT_BINARY_DIR}/openmeeg/tests)

# Installation just copies the directory
install(DIRECTORY  ${OpenMEEG_BINARY_DIR}/wrapping/python/openmeeg
        DESTINATION ${PYTHON_SITE_ARCH}
        DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

if (APPLE)
    set_target_properties(${SWIG_TARGET_NAME} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()
if (WIN32)
    add_custom_command(TARGET ${SWIG_TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:${SWIG_TARGET_NAME}> ${CMAKE_BINARY_DIR}
        COMMAND_EXPAND_LISTS
    )
    add_custom_command(TARGET ${SWIG_TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:OpenMEEG> ${CMAKE_BINARY_DIR}
        COMMAND_EXPAND_LISTS
    )
    add_custom_command(TARGET ${SWIG_TARGET_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:OpenMEEGMaths> ${CMAKE_BINARY_DIR}
        COMMAND_EXPAND_LISTS
    )
endif()

# Test Python wrapping

if (BUILD_TESTING)

    set(PYTHON_TEST_FILES test_vector.py test_matrix.py test_sensors.py test_geom_import.py)

    # TODO: This conditional shouldn't be necessary, but we get:
    # 1/2 Test #689: openmeeg_python_test_python.py ..........***Failed    0.29 sec
    # Could not read the geometry file: D:\a\openmeeg\openmeeg\data\Head1\Head1.geom
    # 2/2 Test #692: openmeeg_python_test_make_geometry.py ...***Failed    0.30 sec
    # Could not read the geometry file: D:\a\openmeeg\openmeeg\data\Head1\Head1.geom
    if(MSVC)
        set(PYTHON_TEST_FILES_REQUIRING_DATA test_python2.py test_mesh.py)
    else()
        set(PYTHON_TEST_FILES_REQUIRING_DATA test_python.py test_python2.py test_mesh.py test_make_geometry.py)
    endif()

    # Add tests requiring data

    foreach (test ${PYTHON_TEST_FILES_REQUIRING_DATA})
        add_test(openmeeg_python_${test} ${Python3_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/wrapping/python/openmeeg/tests/${test} --path ${OpenMEEG_SOURCE_DIR}/data)
        if (USE_MKL AND BUILD_SHARED_LIBS AND NOT MKL_USE_sdl)
            # due to this bug in the MKL
            message(WARNING "Due to a bug in the MKL when used from python: see https://software.intel.com/en-us/forums/intel-distribution-for-python/topic/628976 "
                "workaround: set environment variable:\nexport LD_PRELOAD=\"${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:"
                "${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so\"")
            set_tests_properties(openmeeg_python_${test} PROPERTIES ENVIRONMENT LD_PRELOAD=${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so)
        endif()
        set_tests_properties(openmeeg_python_${test} PROPERTIES ENVIRONMENT PYTHONPATH=${OpenMEEG_BINARY_DIR}/wrapping/python)
    endforeach()

    # Add python tests not requiring data

    foreach (test ${PYTHON_TEST_FILES})
        add_test(openmeeg_python_${test} ${Python3_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/wrapping/python/openmeeg/tests/${test})
        if (USE_MKL AND BUILD_SHARED_LIBS AND NOT MKL_USE_sdl)
            # due to this bug in the MKL
            message(WARNING "Due to a bug in the MKL when used from python: see https://software.intel.com/en-us/forums/intel-distribution-for-python/topic/628976 "
                "workaround: set environment variable:\nexport LD_PRELOAD=\"${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:"
                "${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so\"")
            set_tests_properties(openmeeg_python_${test} PROPERTIES ENVIRONMENT LD_PRELOAD=${MKL_ROOT_DIR}/lib/intel64/libmkl_core.so:${MKL_ROOT_DIR}/lib/intel64/libmkl_sequential.so)
        endif()
        set_tests_properties(openmeeg_python_${test} PROPERTIES ENVIRONMENT PYTHONPATH=${OpenMEEG_BINARY_DIR}/wrapping/python)
    endforeach()

    if (BUILD_REFERENCE_DOC)
        add_test(openmeeg_python_doc ${Python3_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/wrapping/python/openmeeg/tests/test_doc.py)
        set_tests_properties(openmeeg_python_doc PROPERTIES ENVIRONMENT PYTHONPATH=${OpenMEEG_BINARY_DIR}/wrapping/python)
    endif()

endif()

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



list(APPEND CMAKE_SWIG_FLAGS -v -O)
set_source_files_properties(openmeeg.i PROPERTIES CPLUSPLUS ON)
swig_add_library(openmeeg LANGUAGE python SOURCES openmeeg.i)
swig_link_libraries(openmeeg ${PYTHON_LIBRARIES} OpenMEEG::OpenMEEGMaths)

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

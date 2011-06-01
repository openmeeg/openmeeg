#######################################################################
# SWIG Embedding
#######################################################################

IF (ENABLE_PYTHON)

    INCLUDE(FindPythonInterp)

    IF (PYTHON_EXECUTABLE)
        EXEC_PROGRAM("${PYTHON_EXECUTABLE}"
          ARGS "-c \"from distutils.sysconfig import get_python_lib; print get_python_lib(1)\""
          OUTPUT_VARIABLE PYTHON_SITE_ARCH
          RETURN_VALUE PYTHON_SITEARCH_NOT_FOUND)
        STRING(REGEX REPLACE "^.*/(lib|lib64)(/python.*)$" "\\1\\2" PYTHON_REL_SITE_ARCH ${PYTHON_SITE_ARCH})
    ENDIF()

    INCLUDE(FindNumpy.cmake)

    SET (MODULE_NAME openmeeg)
    SET (SWIG_SRC_FILE_NAME ${MODULE_NAME}.i)
    SET (CMAKE_SWIG_FLAGS  "-c++" )

    FIND_PACKAGE(SWIG REQUIRED)
    INCLUDE(${SWIG_USE_FILE})
    FIND_PACKAGE(PythonLibs)

    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${PYTHON_NUMPY_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

    SET(CMAKE_SWIG_FLAGS "")

    FILE(GLOB SWIG_SRC_FILE ${SWIG_SRC_FILE_NAME})

    SET_SOURCE_FILES_PROPERTIES(${SWIG_SRC_FILE_NAME} ${SWIG_SRC_FILE} PROPERTIES CPLUSPLUS ON)

    SET (SWIG_MODULE_${MODULE_NAME}_EXTRA_DEPS ${SWIG_SRC_FILE})
    IF (BUILD_DOCUMENTATION)
        SET(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} -DDOCSTRINGS)
        ADD_CUSTOM_TARGET(docstrings
                          COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/doxy2swig.py ${PROJECT_BINARY_DIR}/documentation/reference/xml/index.xml docstrings.i)
        ADD_DEPENDENCIES(docstrings OpenMEEG-doxygen_outputs)
        ADD_DEPENDENCIES(OpenMEEG-documentation docstrings)
        SET (SWIG_MODULE_${MODULE_NAME}_EXTRA_DEPS ${SWIG_MODULE_${MODULE_NAME}_EXTRA_DEPS} ${PROJECT_BINARY_DIR}/Wrapping/src/docstrings.i)
    ENDIF()

    SWIG_ADD_MODULE(${MODULE_NAME} python ${SWIG_SRC_FILE_NAME} "")

    IF (PYTHON_LIBRARIES)
        TARGET_LINK_LIBRARIES ("_${MODULE_NAME}" ${PYTHON_LIBRARIES} OpenMEEG)
    ENDIF()

    INSTALL(TARGETS  "_${MODULE_NAME}"
            DESTINATION ${PYTHON_REL_SITE_ARCH}
            PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)

    INSTALL(FILES  ${OpenMEEG_BINARY_DIR}/Wrapping/src/${MODULE_NAME}.py
            DESTINATION ${PYTHON_REL_SITE_ARCH}
            PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

    # Test Python wrapping
    IF(BUILD_TESTING)
        EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${OpenMEEG_SOURCE_DIR}/data/demo.py ${OpenMEEG_BINARY_DIR}/Wrapping/src/)
        ADD_TEST(openmeeg_python
                 ${PYTHON_EXECUTABLE} ${OpenMEEG_BINARY_DIR}/Wrapping/src/demo.py --path ${OpenMEEG_SOURCE_DIR}/data)
    ENDIF()

ENDIF()

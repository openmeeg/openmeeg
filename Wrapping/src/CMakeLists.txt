#######################################################################
# SWIG Embedding
#######################################################################

IF (PYTHON_WRAP)

    SET (MODULE_NAME openmeeg)
    SET (SWIG_SRC_FILE_NAME ${MODULE_NAME}.i)
    SET (CMAKE_SWIG_FLAGS  "-c++" )

    IF (APPLE)
        # Darwin ports python
        SET(PYTHON_INCLUDE_PATH /opt/local/include/python2.5/)
        SET(PYTHON_LIBRARIES /opt/local/lib/libpython2.5.dylib)

        # Buildin python
        SET(PYTHON_INCLUDE_PATH /usr/include/python2.5)
        SET(PYTHON_LIBRARIES /usr/lib/libpython2.5.dylib)
    ELSE (APPLE)
        FIND_PACKAGE(PythonLibs)
    ENDIF (APPLE)

    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ${OPENMEEG_INCLUDE_DIRECTORIES})

    FILE(GLOB SWIG_SRC_FILE ${SWIG_SRC_FILE_NAME})

    SET_SOURCE_FILES_PROPERTIES(${SWIG_SRC_FILE_NAME} ${SWIG_SRC_FILE} PROPERTIES CPLUSPLUS ON)

    SET (SWIG_MODULE_${MODULE_NAME}_EXTRA_DEPS ${SWIG_SRC_FILE})
    IF (BUILD_DOCUMENTATION)
        SET(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} -DDOCSTRINGS)
        ADD_CUSTOM_COMMAND(OUTPUT docstrings.i
                           COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/doxy2swig.py ${PROJECT_BINARY_DIR}/documentation/reference/xml/index.xml docstrings.i
                           DEPENDS documentation)
        SET (SWIG_MODULE_${MODULE_NAME}_EXTRA_DEPS ${SWIG_MODULE_${MODULE_NAME}_EXTRA_DEPS} docstrings.i)
    ENDIF (BUILD_DOCUMENTATION)

    SWIG_ADD_MODULE(${MODULE_NAME} python ${SWIG_SRC_FILE_NAME} "")

    SET_SOURCE_FILES_PROPERTIES(${SWIG_SRC_FILE_NAME} PROPERTIES CPLUSPLUS ON)
    SET_SOURCE_FILES_PROPERTIES(${SWIG_SRC_FILE_NAME} PROPERTIES SWIG_FLAGS "-includeall")

    IF (PYTHON_LIBRARIES)
        TARGET_LINK_LIBRARIES ("_${MODULE_NAME}" ${PYTHON_LIBRARIES} ${OPENMEEG_OTHER_LIBRARIES} OpenMEEG)
    ENDIF (PYTHON_LIBRARIES)

    INSTALL(TARGETS  "_${MODULE_NAME}"
            DESTINATION lib/python
            PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ)

    INSTALL(FILES  ${OpenMEEG_BINARY_DIR}/Wrapping/src/${MODULE_NAME}.py
            DESTINATION lib/python
            PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ)

ENDIF (PYTHON_WRAP)


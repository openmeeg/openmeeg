#######################################################################
# SWIG Embedding
#######################################################################

#MESSAGE("VARS: ${CMAKE_REQUIRED_LIBRARIES}")

IF(PYTHON_WRAP)
    INCLUDE(AddFileDependencies)

    SET(CMAKE_SWIG_FLAGS  "-c++" )

    SET ( INTERFACE_FILES openmeeg.i )
    #FILE( GLOB INTERFACE_DEPENDENCIES *.i )
    FILE( GLOB INTERFACE_DEPENDENCIES openmeeg.i )

    SET_SOURCE_FILES_PROPERTIES(${INTERFACE_FILES}
        ${INTERFACE_DEPENDENCIES}
        PROPERTIES
            CPLUSPLUS ON
    )

    SET (MODULE_NAME openmeeg)
    SET (LIBRARY_NAME_PYTHON "${MODULE_NAME}")
    
    SET( SWIG_MODULE_${LIBRARY_NAME_PYTHON}_EXTRA_DEPS ${INTERFACE_DEPENDENCIES} )

    SWIG_ADD_MODULE(${LIBRARY_NAME_PYTHON} python ${INTERFACE_FILES} "" )

    IF ( APPLE )
        # Darwin ports python
        SET(PYTHON_INCLUDE_PATH /opt/local/include/python2.5/)
        SET(PYTHON_LIBRARIES /opt/local/lib/libpython2.5.dylib)

        # Buildin python
        SET(PYTHON_INCLUDE_PATH /usr/include/python2.5)
        SET(PYTHON_LIBRARIES /usr/lib/libpython2.5.dylib)
    ELSE ( APPLE )
        FIND_PACKAGE(PythonLibs)
    ENDIF ( APPLE )

    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ${OPENMEEG_INCLUDE_DIRECTORIES})

    SET_SOURCE_FILES_PROPERTIES(openmeeg.i PROPERTIES CPLUSPLUS ON)
    SET_SOURCE_FILES_PROPERTIES(openmeeg.i PROPERTIES SWIG_FLAGS "-includeall")

    IF (PYTHON_LIBRARIES)
        #MESSAGE("${PYTHON_LIBRARIES}")
        TARGET_LINK_LIBRARIES ("_${LIBRARY_NAME_PYTHON}" ${PYTHON_LIBRARIES} ${OPENMEEG_OTHER_LIBRARIES} OpenMEEG)
    ENDIF (PYTHON_LIBRARIES)

    INSTALL(TARGETS  "_${MODULE_NAME}"
            DESTINATION lib
            PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ
    )

    INSTALL(FILES  ${OpenMEEG_BINARY_DIR}/Wrapping/src/${MODULE_NAME}.py
            DESTINATION lib/python
            PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ
    )

ENDIF(PYTHON_WRAP)


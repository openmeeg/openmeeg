MACRO(ADD_DEMO_TARGET HEAD BINARY_TARGET)
    SET(COPY_TARGET "copy${HEAD}")
    ADD_CUSTOM_TARGET(${COPY_TARGET}
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/Models/${HEAD}
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/IntermediateFiles/${HEAD}
        COMMAND mkdir -p ${CMAKE_BINARY_DIR}/data/Computations/${HEAD}
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/Models/${HEAD} ${CMAKE_BINARY_DIR}/data/Models
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/IntermediateFiles/${HEAD} ${CMAKE_BINARY_DIR}/data/IntermediateFiles
        COMMAND cp -rf ${CMAKE_SOURCE_DIR}/data/Computations/${HEAD} ${CMAKE_BINARY_DIR}/data/Computations)
    SET(DEMO_TARGET "demo${HEAD}")
    ADD_CUSTOM_TARGET(${DEMO_TARGET}
        COMMAND make ${BINARY_TARGET}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/data/Computations/${HEAD})
    ADD_DEPENDENCIES(${DEMO_TARGET} ${COPY_TARGET})
ENDMACRO()

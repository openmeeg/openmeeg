INCLUDE(CompareVersionStrings)
FUNCTION(HDF5_VERIFY_VERSION HDF5MINVERS HDF5_COMP_RES)
    SET(DEFS ${HDF5_DEFINITIONS})
    FOREACH(i ${HDF5_INCLUDE_DIRS})
        SET(DEFS ${DEFS} -I ${i})
    ENDFOREACH()
    TRY_RUN(RRESULT CRESULT
            ${CMAKE_BINARY_DIR}/cmake
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/HDF5Version.c
            COMPILE_DEFINITIONS ${DEFS}
            RUN_OUTPUT_VARIABLE HDF5VERS)
    IF (NOT ${CRESULT})
        MESSAGE(ERROR "Unable to compile a simple hdf5 program. Check your installation.")
    ENDIF()
    IF (NOT ${RRESULT} EQUAL 0)
        MESSAGE(ERROR "Executing a simple hdf5 program.")
    ENDIF()
    COMPARE_VERSION_STRINGS(${HDF5MINVERS} ${HDF5VERS} RVAR)
    SET(${HDF5_COMP_RES} ${RVAR} PARENT_SCOPE)
ENDFUNCTION()

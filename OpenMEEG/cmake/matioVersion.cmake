#   If the matio config file already provided the variable HDF5_VERSION_STRING,
#   there is nothing to do.

if (NOT matio_VERSION_STRING)
    function(matio_GET_VERSION VARNAME)
        try_run(RRESULT CRESULT
                ${CMAKE_BINARY_DIR}/cmake
                ${CMAKE_CURRENT_SOURCE_DIR}/cmake/matioVersion.c
                CMAKE_FLAGS -DINCLUDE_DIRECTORIES::STRING=${matio_INCLUDE_DIRS}
                RUN_OUTPUT_VARIABLE matioVERS)
        if (NOT ${CRESULT})
            message(FATAL " Unable to compile a simple matio program. Check your installation.")
        endif()
        if (NOT ${RRESULT} EQUAL 0)
            message(FATAL " Executing a simple matio program.")
        endif()
        set(${VARNAME} ${matioVERS} PARENT_SCOPE)
    endfunction()
    matio_GET_VERSION(matio_VERSION)
endif()

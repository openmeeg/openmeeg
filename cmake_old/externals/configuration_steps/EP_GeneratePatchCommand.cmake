## #############################################################################
## Check which patch has to be applied
## #############################################################################

function(ep_GeneratePatchCommand ep OutVar)
    set(${ep}_PATCHES_TO_APPLY)

    #   If the source directory of the project does not exist all patches must be applied.

    foreach (patch ${ARGN})
        set(PATCHES_TO_APPLY ${PATCHES_TO_APPLY} ${CMAKE_SOURCE_DIR}/cmake/externals/patches/${patch})
    endforeach()

    #   If the source directory of the project already exists, prune the patch list
    #   and remove those that cannot be applied,

    if (EXISTS ${CMAKE_SOURCE_DIR}/${ep})
        foreach (patch ${PATCHES_TO_APPLY})
            execute_process(COMMAND git apply --ignore-whitespace  --check ${patch}
                            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/${ep}
                            RESULT_VARIABLE   PATCH_OK
                            OUTPUT_QUIET
                            ERROR_QUIET)
            if (PATCH_OK EQUAL 0)
                set(PATCH_LIST ${PATCH_LIST} ${patch})
            endif()
        endforeach()
        set(PATCHES_TO_APPLY ${PATCH_LIST})
    endif()

    set(PATCH_COMMAND)
    if (NOT "${PATCHES_TO_APPLY}" STREQUAL "")
        set(PATCH_COMMAND PATCH_COMMAND git apply --ignore-whitespace ${PATCHES_TO_APPLY})
    endif()

    set(${OutVar} ${PATCH_COMMAND} PARENT_SCOPE)
endfunction()

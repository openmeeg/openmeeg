macro(EP_AddCustomTargets ep)

    set(GIT_COMMAND git checkout master && git pull --ff-only)
    set(SVN_COMMAND svn update)

    if (DEFINED tag)
        set(GIT_COMMAND git fetch ${url} && git stash && git pull origin ${tag} && git stash pop || true)
        set(SVN_COMMAND svn switch ${url} && svn update -r ${tag})
    endif()

    if (EXISTS ${source_dir}/.git) # Git update 
        add_custom_target(update-${ep} 
            COMMAND ${GIT_COMMAND}
            WORKING_DIRECTORY ${source_dir}
            COMMENT "Updating '${ep}' with '${GIT_COMMAND}'")
        set(update-${ep} ON PARENT_SCOPE)

    elseif (EXISTS ${source_dir}/.svn ) ## Svn update 
        add_custom_target(update-${ep} 
            COMMAND ${SVN_COMMAND}
            WORKING_DIRECTORY ${source_dir}
            COMMENT "Updating '${ep}' with '${SVN_COMMAND}'")
        set(update-${ep} ON PARENT_SCOPE)
    endif()

    ## build

    foreach (dependency ${${ep}_dependencies})
        set(build-${ep}_dependences build-${dependency} ${build-${ep}_dependences})
    endforeach()

    add_custom_target(build-${ep} 
        COMMAND cmake --build . --config ${CMAKE_BUILD_TYPE}
        WORKING_DIRECTORY ${binary_dir}
        COMMENT "build '${ep}' with 'cmake --build . --config ${CMAKE_BUILD_TYPE}'"
        DEPENDS ${build-${ep}_dependences})
    set(build-${ep} ON PARENT_SCOPE)
    set(BUILD_TARGETS ${BUILD_TARGETS} build-${ep} PARENT_SCOPE)

    add_custom_target(clean-${ep}
                      COMMAND ${CMAKE_COMMAND} --build . --target clean
                      COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/${ep}/install
                      WORKING_DIRECTORY ${binary_dir})
    set(CLEAN_TARGETS ${CLEAN_TARGETS} clean-${ep} PARENT_SCOPE)

    set(INSTALL_DIRS ${INSTALL_DIRS} ${PROJECT_BINARY_DIR}/${ep}/install PARENT_SCOPE)
                      
    set(extra_parms ${ARGN})
    list(LENGTH extra_parms num_extra_parms)
    if (${num_extra_parms} GREATER 0)
        list(GET extra_parms 0 optional_arg)
        if (BUILD_TESTING AND "${optional_arg}" STREQUAL "TESTS")
            add_custom_target(test-${ep}
                              COMMAND ${CMAKE_COMMAND} --build . --target test
                              DEPENDS build-${ep}
                              WORKING_DIRECTORY ${binary_dir})
            set(TEST_TARGETS ${TEST_TARGETS} test-${ep} PARENT_SCOPE)
        endif()
    endif()
endmacro()

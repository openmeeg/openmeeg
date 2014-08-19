include(Options)

set(GET_REVISIONS_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/externals/get_revisions)

# Add CMAKE_MODULE_PATH to superProjectConfig.cmake

set(${PROJECT_NAME}_CONFIG_FILE "${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake")
file(WRITE ${${PROJECT_NAME}_CONFIG_FILE}
"set(CMAKE_MODULE_PATH 
    ${CMAKE_MODULE_PATH}
    \${CMAKE_MODULE_PATH})\n\n
set(USE_GITHUB_SSH ${USE_GITHUB_SSH})\n"
)
    
# Add path of the get_revisions module to superProjectConfig.cmake
file(APPEND ${${PROJECT_NAME}_CONFIG_FILE}
"set(GET_REVISIONS_MODULE_PATH 
    ${GET_REVISIONS_MODULE_PATH})\n\n"
) 

set(CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/externals
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/externals/configuration_steps
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/externals/projects_modules
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/environment
#    ${CMAKE_CURRENT_SOURCE_DIR}/launchers
#    ${CMAKE_CURRENT_SOURCE_DIR}/packaging
#    ${CMAKE_CURRENT_SOURCE_DIR}/packaging/windows
#    ${CMAKE_CURRENT_SOURCE_DIR}/packaging/linux
#    ${CMAKE_CURRENT_SOURCE_DIR}/packaging/apple
    ${CMAKE_MODULE_PATH}
)

include(SetRevision)
include(GithubOption)
include(SetTargets)
include(CheckEnvironment)
include(ConfigureExternalProjects)
include(InstallPaths)
include(UninstallTargets)
include(TestingTargets)

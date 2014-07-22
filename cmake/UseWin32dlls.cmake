IF (WIN32 AND ENABLE_PACKAGING)

    #FIND_LIBRARY(MSVCP80
    #             NAME msvcp80
    #             PATHS "${OpenMEEG_SOURCE_DIR}/win32addons"
    #             NO_DEFAULT_PATH
    #             NO_CMAKE_ENVIRONMENT_PATH
    #             NO_CMAKE_PATH
    #             NO_SYSTEM_ENVIRONMENT_PATH
    #             NO_CMAKE_SYSTEM_PATH)
    #
    #FIND_LIBRARY(MSVCR80
    #             NAME msvcr80
    #             PATHS "${OpenMEEG_SOURCE_DIR}/win32addons"
    #             NO_DEFAULT_PATH
    #             NO_CMAKE_ENVIRONMENT_PATH
    #             NO_CMAKE_PATH
    #             NO_SYSTEM_ENVIRONMENT_PATH
    #             NO_CMAKE_SYSTEM_PATH)
    #
    #FIND_PATH(WIN32_MANIFEST
    #          NAME Microsoft.VC80.CRT.manifest
    #          PATHS "${OpenMEEG_SOURCE_DIR}/win32addons"
    #          NO_DEFAULT_PATH
    #          NO_CMAKE_ENVIRONMENT_PATH
    #          NO_CMAKE_PATH
    #          NO_SYSTEM_ENVIRONMENT_PATH
    #          NO_CMAKE_SYSTEM_PATH)

    SET(WIN32_MANIFEST "${OpenMEEG_SOURCE_DIR}/win32addons/Microsoft.VC80.CRT.manifest")
    SET(MSVCP80 "${OpenMEEG_SOURCE_DIR}/win32addons/msvcp80.dll")
    SET(MSVCR80 "${OpenMEEG_SOURCE_DIR}/win32addons/msvcr80.dll")

    ADD_CUSTOM_TARGET(copy_dlls ALL
        COMMAND ${CMAKE_COMMAND} -E make_directory ${OpenMEEG_BINARY_DIR}/win32depends/
        COMMAND ${CMAKE_COMMAND} -E copy ${MSVCP80} ${OpenMEEG_BINARY_DIR}/win32depends/
        COMMAND ${CMAKE_COMMAND} -E copy ${MSVCR80} ${OpenMEEG_BINARY_DIR}/win32depends/
        COMMAND ${CMAKE_COMMAND} -E copy ${WIN32_MANIFEST} ${OpenMEEG_BINARY_DIR}/win32depends/
    )

    INSTALL(DIRECTORY ${OpenMEEG_BINARY_DIR}/win32depends/ DESTINATION bin
              PATTERN "${OpenMEEG_BINARY_DIR}/win32depends/*"
              PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                          GROUP_EXECUTE GROUP_READ)

ENDIF()

IF ( WIN32 AND ENABLE_PACKAGING )

    FIND_LIBRARY(MSVCP80
                 NAMES msvcp80)

    FIND_LIBRARY(MSVCR80
                 NAMES msvcr80)

    FIND_PATH(WIN32_MANIFEST Microsoft.VC80.CRT.manifest)

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

ENDIF ( WIN32 AND ENABLE_PACKAGING )
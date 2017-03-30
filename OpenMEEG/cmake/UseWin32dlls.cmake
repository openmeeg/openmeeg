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

set(WIN32_MANIFEST "${PROJECT_SOURCE_DIR}/win32addons/Microsoft.VC80.CRT.manifest")
set(MSVCP80 "${PROJECT_SOURCE_DIR}/win32addons/msvcp80.dll")
set(MSVCR80 "${PROJECT_SOURCE_DIR}/win32addons/msvcr80.dll")

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/win32depends/)

foreach(lib ${MSVCP80} ${MSVCR80} ${WIN32_MANIFEST})
    file(COPY ${lib} DESTINATION ${PROJECT_BINARY_DIR}/win32depends/)
endforeach()

install(DIRECTORY ${PROJECT_BINARY_DIR}/win32depends/ DESTINATION bin
           PATTERN "${PROJECT_BINARY_DIR}/win32depends/*"
           PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                       GROUP_EXECUTE GROUP_READ)

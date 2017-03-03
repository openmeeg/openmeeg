#-----------------------------------------------
# Code Syntax Checking with KWStyle
#-----------------------------------------------

option(USE_KWSTYLE "Checking code syntax using KWStyle" OFF)
mark_as_advanced(USE_KWSTYLE)

if (USE_KWSTYLE)
    add_custom_target(check_syntax
        COMMAND KWStyle -xml ${PROJECT_SOURCE_DIR}/OpenMEEGConfig.kws.xml -html ${PROJECT_BINARY_DIR}/KWStyleCheck -D KWStyleFilesToCheck.txt -v
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endif()


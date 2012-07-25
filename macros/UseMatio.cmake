#------------------------------------------------------------
# MATIO library
#------------------------------------------------------------

OPTION(USE_SYSTEM_MATIO "Build the project using an already installed MATIO" OFF)

IF (USE_SYSTEM_MATIO)
    FIND_LIBRARY(MATIO_LIBRARIES matio)
    IF (NOT MATIO_LIBRARIES)
        MESSAGE(WARNING,"System's matio library not found. Using local version instead.")
    ENDIF()
ENDIF()

IF (NOT MATIO_LIBRARIES)
    SET(MATIO_LIBRARIES matio)
    SET(USE_MATIO 1)
ENDIF()

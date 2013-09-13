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

    #   Not sure to understand why those should be here (this should be made by the matio
    #   subdir). But they are.

    INCLUDE_DIRECTORIES(${OpenMEEG_SOURCE_DIR}/contrib/matio/src)
    INCLUDE_DIRECTORIES(${OpenMEEG_SOURCE_DIR}/contrib/matio/zlib)
    INCLUDE_DIRECTORIES(${OpenMEEG_BINARY_DIR}/contrib/matio/src)
    INCLUDE_DIRECTORIES(${OpenMEEG_BINARY_DIR}/contrib/matio/zlib)

ENDIF()

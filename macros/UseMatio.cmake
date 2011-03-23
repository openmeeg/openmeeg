#------------------------------------------------------------
# MATIO library
#------------------------------------------------------------

OPTION(USE_SYSTEM_MATIO "Build the project using an already installed MATIO" OFF)

IF (USE_SYSTEM_MATIO)
    FIND_LIBRARY(MATIO_LIBRARIES matio)
ENDIF()

IF (NOT MATIO_LIBRARIES)
    INCLUDE_DIRECTORIES(${OpenMEEG_SOURCE_DIR}/contrib/matio/src)
    INCLUDE_DIRECTORIES(${OpenMEEG_SOURCE_DIR}/contrib/matio/zlib)
    INCLUDE_DIRECTORIES(${OpenMEEG_BINARY_DIR}/contrib/matio/src)
    INCLUDE_DIRECTORIES(${OpenMEEG_BINARY_DIR}/contrib/matio/zlib)

    IF(WIN32)
        INCLUDE_DIRECTORIES(${OpenMEEG_SOURCE_DIR}/contrib/matio/contrib/Windows)
    ENDIF()

    SET(MATIO_LIBRARIES matio)
ENDIF()

IF (NOT MATIO_LIBRARIES)
    MESSAGE(ERROR "There is a problem detecting Matio")
ELSE()
    SET(USE_MATIO 1)
ENDIF()

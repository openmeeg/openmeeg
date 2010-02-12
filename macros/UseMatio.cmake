#------------------------------------------------------------
# MATIO library
#------------------------------------------------------------

OPTION(USE_MATIO "Build the project using MATIO" OFF)

IF (USE_MATIO)

    SET(SEARCHPATH
        /usr/lib/
        /usr/lib64/
        /usr/local/lib/atlas
        $ENV{HOME}/local/lib
    )

    FIND_LIBRARY(MATIO_LIB NAMES matio PATHS ${SEARCHPATH})
    FIND_LIBRARY(ZLIB_LIB NAMES z PATHS ${SEARCHPATH})

    SET(MATIO_LIBRARIES ${MATIO_LIB} ${ZLIB_LIB})

    FIND_PATH(MATIO_INCLUDE_PATH matio.h
              /usr/include/
              /usr/local/include/
              $ENV{HOME}/local/include/
    )
    INCLUDE_DIRECTORIES(${MATIO_INCLUDE_PATH})

ENDIF()

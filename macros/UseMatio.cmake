#------------------------------------------------------------
# MATIO library
#------------------------------------------------------------

option(USE_SYSTEM_MATIO "Build the project using an already installed MATIO" OFF)

if (USE_SYSTEM_MATIO)
    find_library(MATIO_LIBRARIES matio)
    if (NOT MATIO_LIBRARIES)
        message(WARNING,"System's matio library not found. Using local version instead.")
    endif()
endif()

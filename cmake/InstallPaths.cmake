#   Offer the user the choice of overriding the installation directories

get_property(LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)

if ("X${LIB64}" STREQUAL "XTRUE" AND NOT APPLE)
    set(LIBSUFFIX 64)
endif()

set(INSTALL_LIB_DIR     lib${LIBSUFFIX} CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin             CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include         CACHE PATH "Installation directory for header files")
set(INSTALL_DATA_DIR    share           CACHE PATH "Installation directory for data files")

if (WIN32)
    set(INSTALL_DATA_DIR ${CMAKE_PROJECT_NAME} CACHE PATH "Installation directory for data files")
endif()

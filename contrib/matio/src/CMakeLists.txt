SET(HEADERS matio_private.h mat4.h mat5.h matio.h matioConfig.h)
#matio_io.c
SET(SRC_FILES mat.c mat4.c mat5.c io.c endian.c matvar_cell.c matvar_struct.c read_data.c snprintf.c)
IF (Z_LIBRARY)
    SET(SRC_FILES ${SRC_FILES} inflate.c)
ENDIF()

IF (WIN32)
    SET(SRC_FILES ${SRC_FILES} matio.rc)
ENDIF()

IF (ENABLE_FORTRAN)
    INCLUDE(FortranCInterface)
    FortranCInterface_HEADER(fortran.h MACRO_NAMESPACE "FC_")
    ADD_SUBDIRECTORY(fortran)
ENDIF()

IF (MAT73)
    SET(SRC_FILES ${SRC_FILES} mat73.c)
ENDIF()

ADD_LIBRARY(matio ${LIB_TYPE} ${SRC_FILES} ${HEADERS})

IF (Z_LIBRARY)
    TARGET_LINK_LIBRARIES(matio ${Z_LIBRARY})
ENDIF()

IF (MAT73)
    TARGET_LINK_LIBRARIES(matio ${HDF5_LIB})
ENDIF()

CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/matioConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/matioConfig.h"
    IMMEDIATE @ONLY)

CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/matio_pubconf.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/matio_pubconf.h"
    IMMEDIATE @ONLY)

# ================
# = INSTALLATION =
# ================

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/matio_pubconf.h ${CMAKE_CURRENT_SOURCE_DIR}/matio.h  ${CMAKE_CURRENT_BINARY_DIR}/fortran.h
        DESTINATION include)

INSTALL(TARGETS matio
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

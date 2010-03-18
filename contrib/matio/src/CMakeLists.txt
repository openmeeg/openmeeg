SET(HEADERS matio_private.h mat4.h mat5.h matio.h matioConfig.h)
SET(SRC_FILES mat.c mat4.c mat5.c io.c endian.c read_data.c snprintf.c)
IF (Z_LIB)
    SET(SRC_FILES ${SRC_FILES} inflate.c)
ENDIF()

IF (ENABLE_FORTRAN)
    INCLUDE(FortranCInterface)
    FortranCInterface_HEADER(fortran.h MACRO_NAMESPACE "FC_" FROM_LIBRARY blas[daxpy])
    SET(SIZEOF_VOID_P ${SIZEOF_VOID_PTR})
    Configure_File(fortran/matio_t.inc.in fortran/matio_t.inc IMMEDIATE @ONLY)
    set(SRC_FILES ${SRC_FILES} fortran/matio_internal.c fortran/matio.f90)
    #set(HEADERS ${HEADERS} fortran/create.f90 fortran/read_data.f90 fortran/write_data.f90 fortran/write.f90)
    INCLUDE_DIRECTORIES(fortran ${CMAKE_CURRENT_BINARY_DIR}/fortran)
ENDIF()

IF(WIN32)
	SET(LIB_TYPE STATIC)
ENDIF()

ADD_LIBRARY(matio ${LIB_TYPE} ${SRC_FILES} ${HEADERS})

IF (Z_LIB)
    TARGET_LINK_LIBRARIES(matio ${Z_LIB})
ENDIF()

CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/matioConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/matioConfig.h"
    IMMEDIATE @ONLY)

# ================
# = INSTALLATION =
# ================

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/matioConfig.h ${CMAKE_CURRENT_SOURCE_DIR}/matio.h
        DESTINATION include)

INSTALL(TARGETS matio
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

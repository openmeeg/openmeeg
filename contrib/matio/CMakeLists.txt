SET(PACKAGE_VERSION_MAJOR "1")
SET(PACKAGE_VERSION_MINOR "5")
SET(PACKAGE_VERSION_PATCH "0")

SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

# ===========
# = Options =
# ===========

OPTION(ENABLE_FORTRAN "Enable fortran interface to mat library")
SET(EXTRA_LANGAGE)
IF (ENABLE_FORTRAN)
    SET(EXTRA_LANGAGE Fortran)
    SET(MATIO_FORTRAN TRUE)
ENDIF()

PROJECT(MATIO C ${EXTRA_LANGAGE})

SET(MATIO_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(MATIO_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# It may seem that ${CMAKE_ARCH_TRIPLET} might appear in a future cmake release.
# For now, do for the best.
SET(MATIO_PLATFORM ${CMAKE_SYSTEM_PROCESSOR}-unknown-${CMAKE_SYSTEM_NAME})

INCLUDE(FindLibrary)

# =======================
# = HANDLE DEFINITTIONS =
# =======================

# ===========
# = Options =
# ===========

INCLUDE(DefaultFileVersion)

OPTION(ENABLE_PROFILE "Turn on profile")
IF (ENABLE_PROFILE)
ENDIF()

OPTION(EXTENDED_SP "Enable extended sparse matrix data types not supported in Matlab" ON)
IF (EXTENDED_SP)
    SET(EXTENDED_SPARSE 1)
ENDIF()

INCLUDE (Documentation OPTIONAL)
IF (BUILD_DOCUMENTATION)
    ADD_SUBDIRECTORY(documentation)
ENDIF()

# =================
# = Check headers =
# =================

#   For public header

# Provide stdint.h for windows, if it is not already there.
# MSVC starts to include stddint.h only from MSVC10 !!!
# Append after in order to use the system includes if they are provided.

IF (WIN32 AND NOT CYGWIN_INSTALL_PATH)
    SET(INTTYPES_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/contrib/msinttypes-r26)
    SET(CMAKE_REQUIRED_INCLUDES ${INTTYPES_INCLUDES})
    INCLUDE_DIRECTORIES(AFTER SYSTEM ${INTTYPES_INCLUDES})
ENDIF()

INCLUDE(CheckIncludeFile)
INCLUDE(CheckHeaderSTDC)

CHECK_INCLUDE_FILE(inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE(stdint.h   HAVE_STDINT_H)
CHECK_INCLUDE_FILE(strings.h  HAVE_STRINGS_H)
CHECK_INCLUDE_FILE(ctype.h    HAVE_CTYPE_H)
CHECK_INCLUDE_FILE(stdlib.h   HAVE_STDLIB_H)
CHECK_INCLUDE_FILE(string.h   HAVE_STRING_H)
#CHECK_INCLUDE_FILE(malloc.h   HAVE_MALLOC_H)
CHECK_HEADER_STDC()

# ===================
# = Check type size =
# ===================

INCLUDE(CheckTypeSize)
SET(TYPES double float long int short char size_t)
FOREACH(TYPE ${TYPES})
    string(TOUPPER ${TYPE} TYPE_UPPER)
    CHECK_TYPE_SIZE(${TYPE} SIZEOF_${TYPE_UPPER})
ENDFOREACH()
CHECK_TYPE_SIZE("long long" LONG_LONG)
SET(HAVE_LONG_LONG ${HAVE_${LONG_LONG}})
SET(SIZEOF_LONG_LONG ${LONG_LONG})
CHECK_TYPE_SIZE("long double" LONG_DOUBLE)
SET(HAVE_LONG_DOUBLE ${HAVE_${LONG_DOUBLE}})
CHECK_TYPE_SIZE(void* SIZEOF_VOID_P)

SET(TYPES uint8_t uint16_t uint32_t uint64_t int8_t int16_t int32_t int64_t)
FOREACH(TYPE ${TYPES})
    string(TOUPPER ${TYPE} TYPE_UPPER)
    CHECK_TYPE_SIZE(${TYPE} ${TYPE_UPPER})
    SET(HAVE_MAT_${TYPE_UPPER} ${HAVE_${TYPE_UPPER}})
    SET(SIZEOF_${TYPE_UPPER} ${${TYPE_UPPER}})
    IF(NOT ${${TYPE_UPPER}})
        MESSAGE(FATAL_ERROR "Unknown type ${TYPE}")
    ENDIF()
ENDFOREACH()

# ===================
# = Check functions =
# ===================

INCLUDE(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(vsnprintf  HAVE_VSNPRINTF)
CHECK_FUNCTION_EXISTS(snprintf   HAVE_SNPRINTF)
CHECK_FUNCTION_EXISTS(vasprintf  HAVE_VASPRINTF)
CHECK_FUNCTION_EXISTS(asprintf   HAVE_ASPRINTF)
CHECK_FUNCTION_EXISTS(strcasecmp HAVE_STRCASECMP)
CHECK_FUNCTION_EXISTS(getopt     HAVE_GETOPT)

INCLUDE(VA_COPY)
VA_COPY()

IF (WIN32)
    IF (NOT LIBRARY_OUTPUT_PATH)
        SET(LIBRARY_OUTPUT_PATH    "${CMAKE_CURRENT_BINARY_DIR}")
    ENDIF()
    IF (NOT EXECUTABLE_OUTPUT_PATH)
        SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}")
    ENDIF()
ENDIF()

IF (NOT HAVE_GETOPT)
    ADD_SUBDIRECTORY(getopt)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/getopt)
    SET(GETOPT_LIB gnu)
ENDIF()

# =========
# = Z Lib =
# =========

OPTION(USE_SYSTEM_ZLIB "Use installed zlib if possible" ON)
MARK_AS_ADVANCED(USE_SYSTEM_ZLIB)

STRING(COMPARE NOTEQUAL "${USE_SYSTEM_ZLIB_STATUS}" "" USE_SYSTEM_ZLIB_STATUS_NOT_EMPTY)
IF (USE_SYSTEM_ZLIB_STATUS_NOT_EMPTY)
    STRING(COMPARE NOTEQUAL "${USE_SYSTEM_ZLIB_STATUS}" "${USE_SYSTEM_ZLIB}" RESET)
ENDIF()

# Store in cache previous value of BUILD_SHARED_LIBS
SET(USE_SYSTEM_ZLIB_STATUS "${USE_SYSTEM_ZLIB}" CACHE INTERNAL "Previous system zlib status" FORCE)

IF (${RESET})
    UNSET(SYSTEM_ZLIB)
    UNSET(SYSTEM_ZLIB CACHE)
    UNSET(Z_LIBRARY)
    UNSET(Z_LIBRARY CACHE)
ENDIF()

IF (USE_SYSTEM_ZLIB)
    MATIO_FIND_LIBRARY(SYSTEM_ZLIB z)
    IF (SYSTEM_ZLIB)
        SET(Z_LIBRARY ${SYSTEM_ZLIB} CACHE STRING "TEST")
    ELSE()
        SET(USE_LOCAL_ZLIB 1)
    ENDIF()
ENDIF()

IF (NOT USE_SYSTEM_ZLIB OR USE_LOCAL_ZLIB)
    SET(ZLIBDIR contrib/zlib-1.2.6)
    IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ZLIBDIR})
        ADD_SUBDIRECTORY(${ZLIBDIR})
        INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/${ZLIBDIR} ${CMAKE_CURRENT_BINARY_DIR}/${ZLIBDIR})
        IF (NOT BUILD_SHARED_LIBS)
            SET(TARGET_SUFFIX -static)
        ENDIF()
        SET(Z_LIBRARY zlib)
    ELSE()
        MESSAGE(FATAL_ERROR "zlib is required by matio")
    ENDIF()
ENDIF()

SET(HAVE_ZLIB 1)
SET(Z_PREFIX)

OPTION(MAT73 "Matlab 7.3 support" ON)

STRING(COMPARE NOTEQUAL "${MAT73_STATUS}" "" MAT73_STATUS_NOT_EMPTY)
IF (MAT73_STATUS_NOT_EMPTY)
    STRING(COMPARE NOTEQUAL "${MAT73_STATUS}" "${MAT73}" RESET)
ENDIF()

# Store in cache previous value of BUILD_SHARED_LIBS
SET(MAT73_STATUS "${MAT73}" CACHE INTERNAL "Previous matlab 7.3 support status" FORCE)

IF (${RESET})
    UNSET(HDF5_FOUND)
    UNSET(HDF5_FOUND CACHE)
ENDIF()

IF (MAT73)
    FIND_PACKAGE(HDF5 1.8)
    IF (HDF5_FOUND)
        SET(HAVE_HDF5 1)
        INCLUDE(HDF5Version)
        HDF5_VERIFY_VERSION(1.8.0 HDF5_OK)
        IF (HDF5_OK EQUAL 1)
            MESSAGE(WARNING "hdf5 library must be 1.8.0 or higher")
            SET(HDF5_FOUND FALSE)
            SET(HDF5_INCLUDE_DIRS)
            SET(HDF5_LIBRARIES)
            UNSET(HAVE_HDF5)
        ENDIF()
    ENDIF()
    IF (NOT HDF5_FOUND)
        SET(HDF5LIBDIR contrib/hdf5-1.8.8)
        IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${HDF5LIBDIR})

            # Tell hdf5 that we are manually overriding certain settings

            SET(HDF5_EXTERNALLY_CONFIGURED 1)

            # Avoid duplicating names of installed libraries

            SET(HDF5_EXTERNAL_LIB_PREFIX "matio")

            ADD_SUBDIRECTORY(${HDF5LIBDIR})
            SET(HDF5_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${HDF5LIBDIR}/src ${CMAKE_CURRENT_BINARY_DIR}/${HDF5LIBDIR})
            SET(HDF5_LIBRARIES hdf5)
            SET(HAVE_HDF5 1)
        ELSE()
            MESSAGE(WARNING "hdf5 not found, mat7.3 format disabled")
            SET(MAT73 OFF CACHE BOOL "Matlab 7.3 support disabled because hdf5 was not found" FORCE)
        ENDIF()
    ENDIF()
    IF (MAT73 AND HAVE_HDF5)
        INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
        IF (NOT BUILD_SHARED_LIBS)
            SET(TARGET_SUFFIX -static)
        ENDIF()
        SET(HDF5_LIB ${HDF5_LIBRARIES})
    ENDIF()
ENDIF()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tools)

# ==============
# = pkg-config =
# ==============

SET(prefix ${CMAKE_INSTALL_PREFIX})
SET(VERSION ${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH})
SET(exec_prefix "\${prefix}")
SET(libdir "\${exec_prefix}/lib")
SET(includedir "\${prefix}/include")
SET(ZLIB_LIBS ${Z_LIBRARY})

CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/matio.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/matio.pc"
    IMMEDIATE @ONLY)

#INSTALL(EXPORT matio DESTINATION )

# ===========
# = Testing =
# ===========

OPTION(MATIO_BUILD_TESTING "Build matio tests" ON)
IF (MATIO_BUILD_TESTING)
    IF (IS_DIRECTORY ${CMAKE_CURREMT_SOURCE_DIR}/test/datasets)
        INCLUDE(CTest)
        ENABLE_TESTING()
        ADD_SUBDIRECTORY(test)
    ENDIF()
ENDIF()

# =============
# = Uninstall =
# =============

CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall-${PROJECT-NAME}
    "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake")

#-----------------------------------------------
# packaging
#-----------------------------------------------

INCLUDE(Packaging)

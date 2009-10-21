#-----------------------------------------------
# tests
#-----------------------------------------------

OPTION(BUILD_TEST "Build tests" OFF)

IF (BUILD_TEST)
    INCLUDE(CTest)
    ENABLE_TESTING()
    MARK_AS_ADVANCED(BUILD_TESTING)
    SUBDIRS(tests)
ENDIF()

IF (USE_GCC AND BUILD_TEST)
    OPTION(ENABLE_COVERAGE "Enable coverage" OFF)
    MARK_AS_ADVANCED(ENABLE_COVERAGE)
ENDIF()

IF (ENABLE_COVERAGE)
    IF (USE_GCC)
        SET(CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs -ftest-coverage")
        SET(CMAKE_CXX_FLAGS "-g -O0 ${GCC_WARNING_OPTIONS} ${CMAKE_EXE_LINKER_FLAGS}")
        SET(CMAKE_C_FLAGS "-g -O0 -Wall -W ${CMAKE_EXE_LINKER_FLAGS}")
    ELSE()
        MESSAGE(SEND_ERROR "Coverage is only available with gcc.")
    ENDIF()
ENDIF()

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
        SET(COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_FLAGS} -lgcov")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS}")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COVERAGE_FLAGS}")
    ELSE()
        MESSAGE(SEND_ERROR "Coverage is only available with gcc.")
    ENDIF()
ENDIF()

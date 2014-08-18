#-----------------------------------------------
# tests
#-----------------------------------------------

option(BUILD_TESTING "Build tests" OFF)

if (BUILD_TESTING)

    if (WIN32)
        file(WRITE ${CMAKE_BINARY_DIR}/TestConfig.cmake "set(ENV{PATH} \"${INSTALL_BIN_DIR}:\$ENV{PATH}\")\n")
        set_directory_properties(PROPERTIES TEST_INCLUDE_FILE "${CMAKE_BINARY_DIR}/TestConfig.cmake")
    endif()

    set(CTEST_BUILD_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
    if (USE_MKL)
        set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-MKL")
    elseif (USE_ATLAS)
        set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-Atlas")
    else()
        set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-BlasLapack")
    endif()
    if (BUILD_SHARED_LIBS)
        set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-shared")
    else()
        set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-static")
    endif()
    set(BUILDNAME ${CTEST_BUILD_NAME})
    include(CTest)
    enable_testing()
    mark_as_advanced(BUILD_TESTING)
    subdirs(tests)
endif()

if (USE_GCC AND BUILD_TESTING)
    option(ENABLE_COVERAGE "Enable coverage" OFF)
    mark_as_advanced(ENABLE_COVERAGE)
endif()

if (ENABLE_COVERAGE)
    if (USE_GCC)
        set(COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_FLAGS} -lgcov")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS}")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COVERAGE_FLAGS}")
    else()
        message(SEND_ERROR "Coverage is only available with gcc.")
    endif()
endif()

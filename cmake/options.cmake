set(BUILD_TESTING True)
enable_testing()

option(USE_VTK "Use VTK" OFF)

option(ENABLE_COVERAGE "Enable coverage" OFF)

option(ENABLE_PYTHON "Enable python bindings" OFF)
set(PYTHON_VERSION 3 CACHE STRING "Python version to use: 2, 2.x, 3, 3.x, or empty")

# Documentation configuration
option(BUILD_DOCUMENTATION "Build doxygen documentation when building all" OFF)
mark_as_advanced(BUILD_DOCUMENTATION)
include(CMakeDependentOption)
cmake_dependent_option(BUILD_REFERENCE_DOC "Build reference documentation" ON
                       "BUILD_DOCUMENTATION" OFF)

if (ENABLE_COVERAGE)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_FLAGS} -lgcov")
    elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(COVERAGE_FLAGS "-coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_FLAGS}")
    else()
        message(SEND_ERROR "Coverage is only available with gcc or clang.")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COVERAGE_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COVERAGE_FLAGS}")
endif()

# Installation options
include(GNUInstallDirs)

include(CTest)

option(USE_VTK   "Use VTK"   OFF)
option(USE_GIFTI "Use GIFTI" OFF)
option(USE_CGAL  "Use CGAL"  OFF)
option(BUILD_SHARED_LIBS "Build Shared Libraries" ON)

option(ENABLE_COVERAGE "Enable coverage" OFF)
option(ENABLE_APPS "Enable app creation" ON)

option(ENABLE_PYTHON "Enable python bindings" OFF)
set(PYTHON_VERSION 3 CACHE STRING "Python version to use: 3, 3.x, or empty")
option(PYTHON_INSTALL_RELATIVE "Make Python install path relative to install-prefix (instead of using Python3_SITEARCH directly)" ON)
option(PYTHON_COPY_RUNTIME_DLLS "Copy runtime DLLs to the cmake build path" OFF)
option(PYTHON_FORCE_EXT_SUFFIX "Force Python extension suffix" OFF)

if (ENABLE_PYTHON)
    option(ENABLE_PYTHON_APPS "Enable python bindings based apps" OFF)
endif()

option(ENABLE_WERROR "Turn on -Werror" OFF)
option(TEST_HEAD3 "Run tests on Head 3" OFF)
option(TEST_SLOW_PYTHON "Run slow tests on Python" OFF)
set(PACKAGE_ARCH_SUFFIX "" CACHE STRING "Package architecture suffix")
set(EXTRA_INSTALL_LIBRARIES "" CACHE STRING "Extra library files to install")

# Documentation configuration

option(BUILD_DOCUMENTATION "Build doxygen documentation when building all" OFF)

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

if (ENABLE_WERROR)
    if (MSVC)
        # warning level 4 and all warnings as errors
        add_compile_options(/W4)
        set(WERROR_COMPILE_OPTION /WX)
    else()
        # lots of warnings and all warnings as errors
        add_compile_options(-Wall -Wextra -pedantic)
        set(WERROR_COMPILE_OPTION -Werror)
    endif()
endif()

# Deal with:
# warning C4530: C++ exception handler used, but unwind semantics are not enabled. Specify /EHsc
if (MSVC)
    add_compile_options(/EHsc)
endif()
if(CMAKE_VS_WINDOWS_TARGET_PLATFORM_VERSION)
    message(STATUS "Selected Windows SDK version ${CMAKE_VS_WINDOWS_TARGET_PLATFORM_VERSION}")
endif()
# Installation options
include(GNUInstallDirs)

# Tell CMake to use rpath with the libs we build

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Set CMAKE_BUILD_TYPE to Release by default.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()
# Set the possible values of build type for cmake-gui
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")

option(BUILD_TOOLS "Build tools" ON)
option(BUILD_TESTING "Build the testing tree" ON)
option(BUILD_DOCUMENTATION "Build the Doxygen documentation" OFF)
option(BUILD_TUTORIALS "Build Tutorials" OFF)
option(ENABLE_PYTHON "Enable Python wrapping" OFF)
option(ENABLE_COVERAGE "Enable Coverage" OFF)

if(APPLE)
    option(APPLE_STANDALONE "Link with static libs to deliver a full standalone package" OFF)
endif()

include(UseOpenMP)
include(UseVTK)
include(UseGifti)
include(UseCGAL)
include(ProgressBar)
include(BlasLapackOption)

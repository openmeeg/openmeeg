option(GIT_HTTPS "Use https protocol to get git repositories." ON)
set(GIT_PREFIX "git@")
if (GIT_HTTPS)
    set(GIT_PREFIX "https://")
endif()

#   Various OpenMEEG options that will be forwarded.

include(OpenMEEGOptions)

#   Various matio options that will be forwarded.

option(MATIO_BUILD_TESTING "Build matio tests" OFF)
option(MATLAB_TESTING "Enable matlab read tests (requires a function matlab)" OFF)

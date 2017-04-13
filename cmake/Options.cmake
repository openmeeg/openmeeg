option(GIT_HTTPS "Use https protocol to get git repositories." ON)
set(GIT_PREFIX "git@")
if (GIT_HTTPS)
    set(GIT_PREFIX "https://")
endif()

# Option: do we want a static or a dynamic build ?
option(BUILD_SHARED_LIBS "Build shared libs for all subprojects" ON)
mark_as_advanced(BUILD_SHARED_LIBS)

#   Various OpenMEEG options that will be forwarded.

include(OpenMEEGOptions)

#   Various matio options that will be forwarded.

option(MATIO_BUILD_TESTING "Build matio tests" OFF)
option(MATLAB_TESTING "Enable matlab read tests (requires a function matlab)" OFF)


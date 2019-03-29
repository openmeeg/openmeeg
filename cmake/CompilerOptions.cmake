#---------------------------------------------------------------
# Test and setup the C++ compiler features
#---------------------------------------------------------------
# XXX Instead of calling CHECK_CXX_STANDARD_LIBRARY with a bunch of
#   unneeded options. The needed functionalities should be checked
#   here

include(CheckCXXFeatures)
CHECK_CXX_STANDARD_LIBRARY()

if (NOT HAVE_ISNORMAL_IN_NAMESPACE_STD)
    include(CheckSymbolExists)
    check_symbol_exists(isnormal math.h HAVE_ISNORMAL_IN_MATH_H)
endif()

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build.")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "None" "Debug" "Release" "RelWithDebInfo" "MinSizeRel")

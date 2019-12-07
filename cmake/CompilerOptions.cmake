#---------------------------------------------------------------
# Test and setup the C++ compiler features
#---------------------------------------------------------------

include(CheckCXXFeatures)

# Workaround the missing support of shared_ptr for arrays.

check_cxx_feature(HAVE_SHARED_PTR_ARRAY_EXTENSION
                  shared_ptr_array_extension.cpp
                  "has C++17 extension of shared_ptr for arrays")

set(SHARED_PTR_DEFINITIONS HAVE_SHARED_PTR_ARRAY_SUPPORT)
if (NOT HAVE_SHARED_PTR_ARRAY_EXTENSION)
    message("Missing support of shared_ptr for arrays. Using boost instead.")
    find_package(Boost 1.53.0 REQUIRED COMPONENTS shared_ptr)
    unset(SHARED_PTR_DEFINITIONS)
endif()

check_cxx_feature(HAVE_ISNORMAL_IN_NAMESPACE_STD
                  isnormal_in_namespace_std.cpp
                  "has isnormal function in namespace std")

if (NOT HAVE_ISNORMAL_IN_NAMESPACE_STD)
    include(CheckSymbolExists)
    check_symbol_exists(isnormal math.h HAVE_ISNORMAL_IN_MATH_H)
endif()

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build.")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "None" "Debug" "Release" "RelWithDebInfo" "MinSizeRel")

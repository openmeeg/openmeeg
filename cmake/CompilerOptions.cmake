if(MSVC)
    # to allow the use of and, or instead of && ||
    add_compile_options("/Za")
endif()

#---------------------------------------------------------------
# Test and setup the C++ compiler features
#---------------------------------------------------------------
# XXX Instead of calling CHECK_CXX_STANDARD_LIBRARY with a bunch of
#   unneeded options. The needed functionalities should be checked
#   here

set(CMAKE_CXX_STANDARD 11)  # use c++11
include(CheckCXXFeatures)
CHECK_CXX_STANDARD_LIBRARY()

if (NOT HAVE_ISNORMAL_IN_NAMESPACE_STD)
    include(CheckSymbolExists)
    check_symbol_exists(isnormal math.h HAVE_ISNORMAL_IN_MATH_H)
endif()

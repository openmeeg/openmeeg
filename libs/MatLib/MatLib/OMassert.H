#pragma once

#include <iostream>
#include <stdexcept>

//  Same as cassert but make sure that assert 
//  are never #defined away by cmake (a property
//  of cmake with makefile generation only ?)
//  And throw an exception that can be caught nicely.

//  Error is the same thing as assert but it is not removed by NDEBUG

#define om_xstr(s) str(s)
#define om_str(s) #s

#if defined(_MSC_VER)
    #define __PRETTY_FUNCTION__ __FUNCSIG__
#endif

#define om_error(expr) \
  ((expr)								\
   ? static_cast<void>(0)						\
   : Assert(om_str(expr),__FILE__,__LINE__,__PRETTY_FUNCTION__))

//  Remove om_assert wgen NDEBUG is defined.
//  Not clear this is useful...

#define IGNORE_NDEBUG
#if defined NDEBUG && !defined IGNORE_NDEBUG
#define om_assert(expr) static_cast<void>(0)
#else
#define om_assert(expr) om_error(expr)
#endif

inline void Assert(const char* assertion,const char* file,const unsigned line,const char* function) {
    std::cerr << file << ':' << line << ':' << function << ": Assertion `" << assertion << "' failed." << std::endl;
    if (!assertion)
        throw std::invalid_argument(assertion);
}

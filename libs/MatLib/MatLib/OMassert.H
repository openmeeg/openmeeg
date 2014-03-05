#pragma once

#include <iostream>
#include <stdexcept>

//  Same as cassert but make sure that assert 
//  are never #defined away by cmake (a property
//  of cmake with makefile generation only ?)
//  And throw an exception that can be caught nicely.

#define om_assert(expr)\
  ((expr)								\
   ? static_cast<void>(0)						\
   : Assert(__STRING(expr),__FILE__,__LINE__,__PRETTY_FUNCTION__))

inline void Assert(const char* assertion,const char* file,const unsigned line,const char* function) {
    std::cerr << file << ':' << line << ':' << function << ": Assertion `" << assertion << "' failed." << std::endl;
    if (!assertion)
        throw std::invalid_argument(assertion);
}

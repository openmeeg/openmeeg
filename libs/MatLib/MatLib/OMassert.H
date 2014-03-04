#pragma once

//  Same as cassert but make sure that assert 
//  are never #defined away by cmake (a property
//  of cmake with makefile generation only ?

#undef NDEBUG
#include <cassert>

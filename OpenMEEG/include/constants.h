// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

#if __cpp_lib_math_constants
#include <numbers>
#endif

namespace OpenMEEG {

    #if __cpp_lib_math_constants
    constexpr double Pi = std::numbers::pi;
    #else
    constexpr double Pi = 3.14159265358979323846;
    #endif
    constexpr double MagFactor = 1e-7;
    constexpr double K = 1.0/(4*Pi);
}

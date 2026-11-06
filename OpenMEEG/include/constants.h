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

#include <numbers>

namespace OpenMEEG {

    constexpr double Pi = std::numbers::pi;
    constexpr double MagFactor = 1e-7;
    constexpr double K = 1.0/(4*Pi);
}

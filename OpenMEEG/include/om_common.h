// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#if defined(__cpp_concepts) && __cpp_concepts >= 201907L && __has_include(<concepts>)
#include <concepts>
#define HAS_CONCEPTS 1
#else
#include <type_traits>
#endif

namespace OpenMEEG {

    // A vector of string is called Strings.

    using Strings = std::vector<std::string>;

    // Compare doubles and floats.

    #ifdef HAS_CONCEPTS
    template <std::floating_point T>
    #else
    template <typename T, typename = std::enable_if_t<std::is_floating_point<T>::value>>
    #endif
    bool almost_equal(const T x,const T y,const T tol=T(1000)*std::numeric_limits<T>::epsilon()) {

        // Exact match (also covers +inf == +inf and -inf == -inf)

        if (x==y)
            return true;

        // NaN, or infinity vs. anything different: never equal

        if (!std::isfinite(x) || !std::isfinite(y))
            return false;

        // Relative test (avoiding a problem around zero).

        return std::abs(x-y)<=tol*std::max({1.0,std::abs(x),std::abs(y)});
    }
}

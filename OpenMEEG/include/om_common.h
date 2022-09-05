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
#include <OMassert.H>
#include <limits>

namespace OpenMEEG {

    // a vector of string is called Strings
    using Strings = std::vector<std::string>;

    // how to compare doubles and floats
    template <typename T>
    bool almost_equal(const T x,const T y,const double eps=1e3) {
        return (std::abs(x-y)<std::numeric_limits<T>::epsilon()*std::abs(x+y)*eps) || std::abs(x-y)<std::numeric_limits<T>::min();
    }
}

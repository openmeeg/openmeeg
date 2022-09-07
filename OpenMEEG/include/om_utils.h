// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <string>
#include <cmath>
#include <random>
#include <iostream>
#include <chrono>
#include <sstream>
#include <algorithm>
#include <cctype>

#include "OpenMEEGConfigure.h"

namespace OpenMEEG {

    inline void dispEllapsed(const std::chrono::duration<double> elapsed_seconds) {
        std::cout <<  "-------------------------------------------" << std::endl
                  <<  "| Elapsed Time: " << elapsed_seconds.count() << " s." << std::endl
                  <<  "-------------------------------------------" << std::endl;
    }
}

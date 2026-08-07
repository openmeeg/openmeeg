// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the GPL open source license.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace OpenMEEG::OMCGAL {

    using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
    using FT      = Kernel::FT;
    using Point3  = Kernel::Point_3;
    using Sphere3 = Kernel::Sphere_3;
}

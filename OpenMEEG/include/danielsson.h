// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <limits>
#include <cmath>

#include <om_common.h>
#include <vertex.h>
#include <mesh.h>
#include <interface.h>
#include <geometry.h>

namespace OpenMEEG {

    double dist_point_cell(const Vect3&,const Triangle&,Vect3&, bool&);
    OPENMEEG_EXPORT std::tuple<double,const Triangle&,const Mesh&>  dist_point_interface(const Vect3&,const Interface&,Vect3&);
    OPENMEEG_EXPORT std::tuple<double,const Triangle&,const Mesh&,const Interface&> dist_point_geom(const Vect3&,const Geometry&,Vect3&);
}

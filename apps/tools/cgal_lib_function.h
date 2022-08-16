// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

// for verbosity
#define CGAL_MESH_3_VERBOSE

#include <mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/config.h>
#include <CGAL/version.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::FT FT;
typedef CGAL::Polyhedron_3<K>  Polyhedron;
typedef Polyhedron::HalfedgeDS HDS;

namespace OpenMEEG {

    namespace {
        template <typename C3t3>
        Mesh CGAL_to_OM(C3t3 c3t3);
    }

    Mesh cgal_mesh_function(double sphere_radius, double hemisphere, double radius_bound, double distance_bound);
}

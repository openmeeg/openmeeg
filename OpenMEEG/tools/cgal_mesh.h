/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

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

    #if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(4,6,0)
    namespace {
        class Build_Mesh2Polyhedron : public CGAL::Modifier_base<HDS> {
            const Mesh *m;
        public:
            Build_Mesh2Polyhedron(const Mesh& m): m(&m) {}
            void operator()(HDS& target);
        };

        Polyhedron Mesh2Polyhedron(const Mesh& m);
    }

    Mesh cgal_refine(const Mesh& m_in, double radius_bound, double distance_bound, const char* sizing_field);

    Mesh cgal_decimate(const Mesh& m_in, unsigned nb_points);
    #endif

    #ifdef CGAL_ImageIO
    Mesh cgal_mesh_3Dlabeled_image(const char* input_filename, double radius_bound, double distance_bound);
    Mesh cgal_mesh_3Dlevelset_image(const char* input_filename,  double levelset_value, bool positive_inside, double radius_bound, double distance_bound);
    #endif

}

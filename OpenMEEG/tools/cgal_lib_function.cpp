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

#include <cgal_lib.h>

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// sphere function
// for spheres
class FT_to_point_Sphere_wrapper : public std::unary_function<Point_3, FT>
{
    double sqrd;
    public:
    FT_to_point_Sphere_wrapper(FT sqrd_) : sqrd(sqrd_) {}
    FT operator()(Point_3 p) const { return (std::pow(p.x(),2)+std::pow(p.y(),2)+std::pow(p.z(),2)-sqrd); }
};

// Hemisphere function
class FT_to_point_HemiSphere_wrapper : public std::unary_function<Point_3, FT>
{
    double sqrd;
    public:
    FT_to_point_HemiSphere_wrapper(FT sqrd_) : sqrd(sqrd_) {}
    FT operator()(Point_3 p) const {
        double d_sphere = (std::pow(p.x(), 2) + std::pow(p.y(), 2) + std::pow(p.z(), 2) - sqrd);
        double d_total;
        if ( p.z() > 0 ) {
            d_total = ( d_sphere > 0 ) ? d_sphere : -std::min(-d_sphere, std::pow(p.z(), 2));
        } else {
            d_total = std::min(std::pow(p.x(), 2), sqrd ) + std::min(std::pow(p.y(), 2), sqrd ) + std::pow(p.z(), 2);
        }
        return d_total;
    }
};

// typedef to mesh a function
typedef FT_to_point_Sphere_wrapper SphereFunction;
typedef CGAL::Implicit_mesh_domain_3<SphereFunction, K> SphereDomain;
typedef FT_to_point_HemiSphere_wrapper HemiSphereFunction;
typedef CGAL::Implicit_mesh_domain_3<HemiSphereFunction, K> HemiSphereDomain;
typedef CGAL::Mesh_triangulation_3<SphereDomain>::type TrS;
typedef CGAL::Mesh_complex_3_in_triangulation_3<TrS> C3t3S;
// Criteria
typedef CGAL::Mesh_criteria_3<TrS> Mesh_criteriaS;

namespace OpenMEEG {
    /// decimate safely a mesh
    Mesh cgal_mesh_function(double sphere_radius, double hemisphere_radius, double radius_bound, double distance_bound)
    {
        // defining the sphere domain
        SphereFunction spherefunction(std::pow(sphere_radius, 2));
        SphereDomain sdomain(spherefunction, K::Sphere_3(CGAL::ORIGIN, std::pow(1.1*sphere_radius, 2)), 1e-6); // with its bounding sphere
        // defining the hemisphere domain
        HemiSphereFunction hemispherefunction(std::pow(hemisphere_radius, 2));
        HemiSphereDomain hdomain(hemispherefunction, K::Sphere_3(TrS::Point(0, 0, hemisphere_radius/2.), std::pow(1.1*hemisphere_radius, 2)), 1e-6); // with its bounding sphere

        // Mesh criteria
        Mesh_criteriaS criteria(facet_angle=30, facet_size=radius_bound, facet_distance=distance_bound);

        // meshing domain
        C3t3S c3t3;

        if ( sphere_radius > 0.0001 ) {
            c3t3 = CGAL::make_mesh_3<C3t3S>(sdomain, criteria, no_exude(), no_perturb());
        } else {
            // if you want want to add initial points on the hemisphere circle (for a better definition),
            // have a look here (it probably needs to construct the facets also ).
# if 0
            std::pair<Tr::Point,unsigned> p[init_points];
            for ( unsigned iip = 0; iip < init_points; ++iip) {
                p[iip] = std::make_pair(Tr::Point(hemisphere_radius*std::cos(2.*M_PI/init_points*iip), hemisphere_radius*std::sin(2.*M_PI/init_points*iip) , 0),0);
            }
            c3t3.insert_surface_points(&p[0],&p[init_points-1]);
            CGAL::refine_mesh_3<C3t3>(c3t3, hdomain, criteria, no_exude(), no_perturb());
#else
            c3t3 = CGAL::make_mesh_3<C3t3S>(hdomain, criteria, no_exude(), no_perturb());
#endif
        }

        return CGAL_to_OM(c3t3);
    }
}

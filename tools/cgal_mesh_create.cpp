/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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

#include "mesh.h"
#include <assert.h>
#include "options.h"
#include "cgal_mesh.h"

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

using namespace OpenMEEG;

int main(int argc, char **argv) {
    command_usage("Create a BEM mesh from either a very fine mesh, a 3D image (.inr format v4), or sphere (no input file):");
    const char * input_filename  = command_option("-i",(const char *) NULL,"Input image or mesh");
    const bool   inout           = command_option("-inout",false,"Inside out the image ?");
    const double radius_bound    = command_option("-fs",1e-1,"facet radius bound of elements");
    const double distance_bound  = command_option("-fd",1e-1,"facet distance bound to the input surface");
    const char * output_filename = command_option("-o",(const char *) NULL,"Output Mesh");
    const double sphere_radius   = command_option("-r", 0., "radius of the sphere");
    const double hemisphere_radius = command_option("-hr", 0., "radius of the hemisphere");
    const unsigned init_points   = command_option("-ip", 10, "initial number of points");

    if ( command_option("-h",(const char *)0,0) ) { 
        return 0; 
    }
    if ( output_filename == NULL ) {
        std::cerr << "Set an output filename" << std::endl;
        return 0;
    }

    Gray_level_image *image;

    Tr tr;           // 3D-Delaunay triangulation

    C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation

    // Mesh criteria
    double angle_bound = 30.; // angle bound for the triangular elements (not less than 30 degree)
    if ( input_filename != NULL ) {
        std::string extension = getNameExtension(input_filename);
        if ( extension == "inr" || extension == "inr4" ) {
            bool positive_inside  = inout;
            double value_outside  = 1.;
            double levelset_value = 0.;
            image = new Gray_level_image(input_filename,levelset_value,positive_inside,value_outside); 
            std::cout << "Input INR image:\n dimension: " << image->xdim() << "x"<< image->ydim() << "x"<< image->zdim() << std::endl;
            // Carefully choose bounding sphere: the center must be inside the
            // surface defined by 'image' and the radius must be high enough so that
            // the sphere actually bounds the whole image.
            Point_3 bounding_sphere_center(image->xdim()/2., image->ydim()/2., image->zdim()/2.);
            Geom_traits::FT bounding_sphere_squared_radius = image->xdim()*image->ydim()*2.;
            Geom_traits::Sphere_3 bounding_sphere(bounding_sphere_center, bounding_sphere_squared_radius);
            // definition of the surface, with 10^-5 as relative precision
            GreySurface_3 surface(*image, bounding_sphere, 1e-8);
            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle_bound, radius_bound, distance_bound);
            CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag()); // make the surfacic mesh with a manifold criteria
        } else {
            Mesh m_in(input_filename, false);
            std::cout << "Input surface:\n nb of points: " << m_in.nb_vertices() << "\t nb of triangles:\t" << m_in.nb_triangles() << std::endl;
            // Create input polyhedron
            Polyhedron polyhedron;
            for ( Mesh::const_iterator tit = m_in.begin(); tit != m_in.end(); ++tit) {
                const Vertex p1 = tit->s1();
                const Vertex p2 = tit->s2();
                const Vertex p3 = tit->s3();
                polyhedron.make_triangle(Point_3(p1.x(), p1.y(), p1.z()), Point_3(p2.x(), p2.y(), p2.z()), Point_3(p3.x(), p3.y(), p3.z()));
            }
            // Create domain
            Polyhedral_mesh_domain poly_mesh(polyhedron);

            unsigned nb_initial_points = 5;
            unsigned i = 0; 
            for ( Polyhedron::Vertex_iterator it = polyhedron.vertices_begin(); it != polyhedron.vertices_end(), i < nb_initial_points; it++, i++) {
                if ( i%m_in.nb_vertices() / 5 == 0 ) {
                    tr.insert(it->point());
                }
            }
            // Mesh generation
            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle_bound, radius_bound, distance_bound);
            CGAL::make_surface_mesh(c2t3, poly_mesh, poly_mesh, criteria, CGAL::Manifold_tag()); // make the surfacic mesh with a manifold criteria
        }
    } else { // implicit functions (sphere or hemisphere)
        if ( sphere_radius > 0.0001 ) {
            SphereFunction spherefunction(std::pow(sphere_radius, 2));
            // defining the surface
            SphereSurface_3 surface(spherefunction,             // pointer to function
                    Geom_traits::Sphere_3(CGAL::ORIGIN, std::pow(2.*sphere_radius, 2))); // bounding sphere
            // defining meshing criteria
            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle_bound,  // angular bound
                    radius_bound,  // radius bound
                    distance_bound); // distance bound
            // meshing surface
            CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
        } else if ( hemisphere_radius > 0.0001 ) { // bug and too many point close to the origin (0, 0, 0) if BB is centered there
            HemiSphereFunction hspherefunction(std::pow(hemisphere_radius, 2));
            // defining the surface
            HemiSphereSurface_3 surface(hspherefunction,             // pointer to function
                    Geom_traits::Sphere_3(Tr::Point(0, 0, hemisphere_radius/2.), std::pow(2.*hemisphere_radius, 2))); // bounding sphere
            // defining meshing criteria
            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle_bound,  // angular bound
                    radius_bound,  // radius bound
                    distance_bound); // distance bound

            // Setting initial points on the circle
            unsigned init_points2 = init_points;
            if ( 2.*M_PI*hemisphere_radius/init_points > radius_bound ) {
                init_points2 = unsigned(2.*M_PI*hemisphere_radius / radius_bound);
                std::cout << "Increasing nb of initial points from " << init_points << " to " << init_points2 << std::endl;
            }

            for ( unsigned iip = 1; iip <= init_points2; ++iip) {
                Tr::Point p = Tr::Point(hemisphere_radius*std::cos(2.*M_PI/init_points2*iip), hemisphere_radius*std::sin(2.*M_PI/init_points2*iip) , 0);
                tr.insert(p);
            }
            // meshing surface
            CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag(), 10); // 10 is nb of initial point (chosen by CGAL!...)
            // project the disk point on the x-y plane (i.e assign a 0 z-coordinate strictly)
            for ( Tr::Finite_vertices_iterator vit = c2t3.triangulation().finite_vertices_begin(), end = c2t3.triangulation().finite_vertices_end(); vit != end; ++vit)
            {
                Point_3& p = vit->point();
                double d_sphere = (std::pow(p.x(), 2) + std::pow(p.y(), 2) + std::pow(p.z(), 2) - hemisphere_radius); 

                if ( ( std::abs(p.z()) < std::min(distance_bound, radius_bound) ) && ( d_sphere < std::min(distance_bound, radius_bound)) ) {
                    p = Point_3(p.x(), p.y(), 0);
                }
            }
        }
    }

    Mesh m_out = CGAL_to_OM(c2t3);

    m_out.update();
    m_out.save(output_filename);
    m_out.info();
    return 0;
}

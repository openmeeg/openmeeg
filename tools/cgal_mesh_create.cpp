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

#include "mesh3.h"
#include <assert.h>
#include "options.h"
#include "cgal_mesh_include.h"

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

int main(int argc, char **argv) {
    command_usage("Create a BEM mesh from either a very fine mesh, a 3D image (.inr format v4), or sphere (no input file):");
    const char *input_filename = command_option("-i",(const char *) NULL,"Input image or mesh");
    const double radius_bound = command_option("-fs",1e-1,"facet radius bound of elements");
    const double distance_bound = command_option("-fd",1e-1,"facet distance bound to the input surface");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output Mesh");
    const double sphere_radius = command_option("-r", 1., "radius of the sphere");

    if (command_option("-h",(const char *)0,0)) return 0;
    if (output_filename == NULL) {
        std::cerr << "Set an output filename" << std::endl;
        return 0;
    }

    Gray_level_image *image;

    Tr tr;           // 3D-Delaunay triangulation
    C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation

    // Mesh criteria
    double angle_bound = 30.; // angle bound for the triangular elements (not less than 30 degree)
    if (input_filename != NULL) {
        std::string extension = getNameExtension(input_filename);
        if ( extension == "inr" || extension == "inr4") {
            bool positive_inside  = false;
            double value_outside  = 1.;
            double levelset_value = 0.;
            image = new Gray_level_image(input_filename,levelset_value,positive_inside,value_outside); 
            std::cout << "Input INR image:\n dimension: " << image->xdim() << "x"<< image->ydim() << "x"<< image->zdim() << std::endl;
            // Carefully choosen bounding sphere: the center must be inside the
            // surface defined by 'image' and the radius must be high enough so that
            // the sphere actually bounds the whole image.
            Point_3 bounding_sphere_center(image->xdim()/2., image->ydim()/2., image->zdim()/2.);
            Geom_traits::FT bounding_sphere_squared_radius = image->xdim()*image->ydim()*2.;
            Geom_traits::Sphere_3 bounding_sphere(bounding_sphere_center, bounding_sphere_squared_radius);
            // definition of the surface, with 10^-5 as relative precision
            Surface_3 surface(*image, bounding_sphere, 1e-8);
            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle_bound, radius_bound, distance_bound);
            CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag()); // make the surfacic mesh with a manifold criteria
        } 
        else 
        {
            Mesh m_in;
            m_in.load(input_filename, false);
            std::cout << "Input surface:\n nb of points: " << m_in.nbPts() << "\t nb of triangles:\t" << m_in.nbTrgs() << std::endl;
            // Create input polyhedron
            Polyhedron polyhedron;
            for (int it = 0; it < m_in.nbTrgs(); it++) {
                Vect3 p1 = m_in.point(m_in.triangle(it).s1()), p2 = m_in.point(m_in.triangle(it).s2()), p3 = m_in.point(m_in.triangle(it).s3());
                polyhedron.make_triangle(Point_3(p1.x(), p1.y(), p1.z()),Point_3(p2.x(), p2.y(), p2.z()), Point_3(p3.x(), p3.y(), p3.z()));
            }

            Input_surface input_surf(polyhedron);

            unsigned nb_initial_points = 5;
            unsigned i = 0; 
            for(Polyhedron::Vertex_iterator it = polyhedron.vertices_begin(); it != polyhedron.vertices_end(), i < nb_initial_points; it++, i++) {
                if (i%m_in.nbPts()/5 == 0)
                    tr.insert(it->point());
            }
            // Mesh generation
            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angle_bound, radius_bound, distance_bound);
            CGAL::make_surface_mesh(c2t3, input_surf, input_surf, criteria, CGAL::Manifold_tag()); // make the surfacic mesh with a manifold criteria
        }
    }
    else 
    {
        if (sphere_radius > 0.) {
            SphereFunction spherefunction(std::pow(sphere_radius,2));
            // defining the surface
            SphereSurface_3 surface(spherefunction,             // pointer to function
                    Geom_traits::Sphere_3(CGAL::ORIGIN, std::pow(2.*sphere_radius, 2))); // bounding sphere

            // defining meshing criteria
            CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                    radius_bound,  // radius bound
                    distance_bound); // distance bound
            // meshing surface
            CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
        }
    }

    Mesh m_out = CGAL_to_OM(c2t3);

    m_out.make_links();
    m_out.update_triangles();
    m_out.recompute_normals(); // Compute normals since off files don't have any !
    if (!m_out.has_correct_orientation()) {
        m_out.flip_faces();
        std::cout << "Global orientation problem corrected" << std::endl;
    }
    m_out.has_correct_orientation();
    m_out.save(output_filename);
    m_out.info();
    return 0;
}

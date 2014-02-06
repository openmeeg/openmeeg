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
#include <CGAL/Polyhedral_mesh_domain_3.h>

using namespace OpenMEEG;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

typedef CGAL::Polyhedron_3<K> Polyhedron;
// Domain 
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

typedef K::Point_3 Point_3;
typedef K::FT FT;

int main(int argc, char **argv) {
    command_usage("Re-mesh a mesh:");
    const char * input_filename  = command_option("-i",(const char *) NULL,"Input image or mesh");
    const double radius_bound    = command_option("-fs",1e-1,"facet radius bound of elements");
    const double distance_bound  = command_option("-fd",1e-1,"facet distance bound to the input surface");
    const char * output_filename = command_option("-o",(const char *) NULL,"Output Mesh");

    if ( command_option("-h",(const char *)0,0) ) { 
        return 0; 
    }
    if ( output_filename == NULL ) {
        std::cerr << "Set an output filename" << std::endl;
        return 0;
    }

    // Mesh input
    Mesh m_in(input_filename, false);
    std::cout << "Input surface:\n nb of points: " << m_in.nb_vertices() << "\t nb of triangles:\t" << m_in.nb_triangles() << std::endl;

    // Mesh criteria
    Mesh_criteria criteria(facet_angle=30, facet_size=radius_bound, facet_distance=distance_bound);

    // Create input polyhedron
    Polyhedron polyhedron;
    for ( Mesh::const_iterator tit = m_in.begin(); tit != m_in.end(); ++tit) {
        const Vertex p1 = tit->s1();
        const Vertex p2 = tit->s2();
        const Vertex p3 = tit->s3();
        polyhedron.make_triangle(Point_3(p1.x(), p1.y(), p1.z()), Point_3(p2.x(), p2.y(), p2.z()), Point_3(p3.x(), p3.y(), p3.z()));
    }
    // Create domain
    Mesh_domain domain(polyhedron);

    C3t3 c3t3;
    unsigned i = 0;
    unsigned nb_initial_points = 5;
    for ( Polyhedron::Vertex_iterator it = polyhedron.vertices_begin(); it != polyhedron.vertices_end(); it++, i++) {
        if ( i% (m_in.nb_vertices() / (1+nb_initial_points)) == 0 ) {
            c3t3.triangulation().insert(it->point());
        }
    }
    // Meshing
    CGAL::refine_mesh_3(c3t3, domain, criteria, no_exude(), no_perturb());

    Mesh m_out = CGAL_to_OM(c3t3);
    m_out.save(output_filename);
    m_out.info();
    return 0;
}

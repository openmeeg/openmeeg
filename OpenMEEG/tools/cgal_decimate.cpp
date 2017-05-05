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

// Mesh simplification adapted from Kai Dang 2015 @ Inria

#include <mesh.h>
#include <options.h>
#include "cgal_mesh.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

namespace SMS = CGAL::Surface_mesh_simplification ;

using namespace OpenMEEG;

int main(int argc, char **argv)
{
    command_usage("Decimate a mesh:");
    const char * input_filename  = command_option("-i", (const char *) NULL,"Input image or mesh");
    const int    nb_points       = command_option("-n", 1000,"desired output number of vertices");
    const char * output_filename = command_option("-o", (const char *) NULL,"Output Mesh");

    Mesh m_in(input_filename, false);
    std::cout << "Input surface:\n nb of points: " << m_in.nb_vertices() << "\t nb of triangles:\t" << m_in.nb_triangles() << std::endl;
    Polyhedron polyhedron = Mesh2Polyhedron(m_in);

    int nb_edges = 3*nb_points-6;

    std::cout << "Target: less than " << nb_points << " vertices <=> " << nb_edges << " edges.\n";

    SMS::Count_stop_predicate<Polyhedron> stop(nb_edges);

    int r = SMS::edge_collapse(
            polyhedron,
            stop,
            CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, polyhedron)).
            halfedge_index_map(get(CGAL::halfedge_external_index, polyhedron)).
            get_cost(SMS::Edge_length_cost<Polyhedron>()).
            get_placement(SMS::Midpoint_placement<Polyhedron>()));

    std::cout << "polyhedron.size_of_vertices() = " <<  polyhedron.size_of_vertices() << std::endl;
    std::cout << "polyhedron.size_of_facets() = " <<  polyhedron.size_of_facets() << std::endl;

    unsigned inum = 0;
    // write the output mesh
    Mesh m(polyhedron.size_of_vertices(), polyhedron.size_of_facets());
    std::map< Polyhedron::Vertex_const_handle, unsigned> V;
    for ( Polyhedron::Vertex_iterator vit = polyhedron.vertices_begin(); vit != polyhedron.vertices_end(); ++vit)
    {
        const Polyhedron::Vertex::Point& p = vit->point();
        m.add_vertex(Vertex(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())));
        V[vit] = inum++;
    }

    for ( Polyhedron::Facet_iterator fit = polyhedron.facets_begin(); fit != polyhedron.facets_end(); ++fit) {
        Polyhedron::Facet::Halfedge_around_facet_circulator j = fit->facet_begin();
        // Facets in polyhedral surfaces are triangles.
        CGAL_assertion( CGAL::circulator_size(j) == 3);
        const int index1 = V[j->vertex()];
        j++;
        const int index2 = V[j->vertex()];
        j++;
        const int index3 = V[j->vertex()];
        Triangle t(m.vertices()[index1], m.vertices()[index2], m.vertices()[index3] );
        m.push_back(t);
    }

    m.update();
    m.correct_global_orientation();

    std::cout << "Finished...\n" << r << " edges removed.\n" << m.nb_vertices() << " final nb_points.\n";
    m.info();
    m.save(output_filename);

    return 0 ; 
}

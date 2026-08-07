// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the GPL open source license.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include "OM_CGAL.h"
#include "Conversions.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>

// Simplification function

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policies

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

namespace OpenMEEG::OMCGAL {

    /// Safely decimate a mesh.

    Mesh decimate_mesh(const Mesh& mesh,const unsigned nb_points) {

        namespace params = CGAL::parameters;
        namespace SMS    = CGAL::Surface_mesh_simplification;

        // Create a working copy of the input mesh.

        SurfaceMesh cgal_mesh = OM_to_CGAL_SurfaceMesh(mesh);

        const unsigned nb_edges = 3*nb_points-6; // Assumes closed surfaces.

        std::cout << "Target: less than " << nb_points << " vertices <=> " << nb_edges << " edges.\n";

        SMS::Edge_count_stop_predicate<SurfaceMesh> stop(nb_edges);

        const int nb_edge_removed = SMS::edge_collapse(cgal_mesh,stop,params::get_cost(SMS::Edge_length_cost<SurfaceMesh>()).get_placement(SMS::Midpoint_placement<SurfaceMesh>()));

        std::cout << "Finished..." << std::endl
                  << "Number of vertices  = " << cgal_mesh.number_of_vertices() << std::endl
                  << "Number of triangles = " << cgal_mesh.number_of_faces()    << std::endl
                  << "Edges removed: " << nb_edge_removed << std::endl;

        return CGAL_SurfaceMesh_to_OM(cgal_mesh);
    }
}

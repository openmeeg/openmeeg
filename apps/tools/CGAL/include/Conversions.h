// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the GPL open source license.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <map>

#include <mesh.h>
#include <geometry.h>

#include <CGAL/Surface_mesh.h>

#include "Defs.h"

namespace OpenMEEG::OMCGAL {

    using SurfaceMesh = CGAL::Surface_mesh<Point3>;

    // Convert an OpenMEEG mesh to a CGAL SurfaceMesh.

    inline
    SurfaceMesh OM_to_CGAL_SurfaceMesh(const Mesh& mesh) {

        SurfaceMesh cgal_mesh;

        std::map<const Vertex*,SurfaceMesh::Vertex_index> vmap;
        for (auto vit=mesh.vertices().begin(); vit!=mesh.vertices().end(); ++vit)
            vmap[*vit] = cgal_mesh.add_vertex(Point3((*vit)->x(),(*vit)->y(),(*vit)->z()));

        for (auto tit=mesh.triangles().begin(); tit!=mesh.triangles().end(); ++tit)
            cgal_mesh.add_face(vmap[&tit->vertex(0)],vmap[&tit->vertex(1)],vmap[&tit->vertex(2)]);

        return cgal_mesh;
    }

    // Convert a CGAL SurfaceMesh to an OpenMEEG mesh.

    inline
    Mesh CGAL_SurfaceMesh_to_OM(const SurfaceMesh& cgal_mesh) {

        Mesh mesh(cgal_mesh.number_of_vertices(),cgal_mesh.number_of_faces());
        Geometry& geom = mesh.geometry();

        // Insert points in the geometry (assumes that all points are used).

        std::map<SurfaceMesh::Vertex_index,unsigned> vertex_index_mapping;
        for (SurfaceMesh::Vertex_index v : cgal_mesh.vertices()) {

            using CGAL::to_double;

            const Point3&     p     = cgal_mesh.point(v);
            const unsigned    index = geom.add_vertex(Vertex(to_double(p.x()),to_double(p.y()),to_double(p.z())));
            vertex_index_mapping[v] = index;
            mesh.vertices().push_back(&geom.vertices()[index]);
        }

        // Insert triangles in the mesh.

        for (SurfaceMesh::Face_index f : cgal_mesh.faces()) {

            unsigned k = 0;
            unsigned inds[3];
            for (auto v : CGAL::vertices_around_face(cgal_mesh.halfedge(f),cgal_mesh))
                inds[k++] = vertex_index_mapping[v];

            mesh.add_triangle(TriangleIndices(inds));
        }

        mesh.update();
        mesh.correct_global_orientation();

        return mesh;
    }

    // Convert CGAL meshes to OpenMEEG meshes.

    template <typename C3t3>
    Mesh CGAL_to_OM(const C3t3& c3t3) {

        Mesh mesh(c3t3.triangulation().number_of_vertices(),c3t3.number_of_facets());
        Geometry& geom = mesh.geometry();

        using Triangulation = C3t3::Triangulation;

        std::map<typename C3t3::Vertex_handle,unsigned> vertex_index_mapping;
        for (typename Triangulation::Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin(); vit!=c3t3.triangulation().finite_vertices_end(); ++vit) {

            using CGAL::to_double;
            
            const typename Triangulation::Point& p = vit->point();
            const unsigned index = geom.add_vertex(Vertex(to_double(p.x()),to_double(p.y()),to_double(p.z())));
            vertex_index_mapping[vit] = index;
            mesh.vertices().push_back(&geom.vertices()[index]);
        }

        for (typename C3t3::Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin(); fit!=c3t3.facets_in_complex_end(); ++fit) {
            const typename Triangulation::Cell_handle cell = fit->first;
            const int& index = fit->second;
            const unsigned index1 = vertex_index_mapping[cell->vertex(c3t3.triangulation().vertex_triple_index(index,0))];
            const unsigned index2 = vertex_index_mapping[cell->vertex(c3t3.triangulation().vertex_triple_index(index,1))];
            const unsigned index3 = vertex_index_mapping[cell->vertex(c3t3.triangulation().vertex_triple_index(index,2))];
            mesh.add_triangle({ index1, index2, index3 });
        }

        mesh.update();
        mesh.info();

        mesh.correct_global_orientation();

        return mesh;
    }
}

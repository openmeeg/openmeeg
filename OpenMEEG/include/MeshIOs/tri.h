// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <iostream>
#include <string>

#include <om_utils.h>
#include <MeshIO.h>

namespace OpenMEEG::MeshIOs {

    /// \brief Mesh io for TRI file format.

    class OPENMEEG_EXPORT Tri: public MeshIO {

        typedef MeshIO base;

    public:

        void load_points(Geometry& geom) override {
            char ch;
            unsigned npts;
            fs >> ch >> npts;

            Vertices vertices;
            for (unsigned i=0; i<npts; ++i) {
                Vertex v;
                Normal n;
                fs >> v >> n;
                vertices.push_back(v);
            }
            indmap = geom.add_vertices(vertices);
        }

        void load_triangles(OpenMEEG::Mesh& mesh) override {
            reference_vertices(mesh);

            char ch;
            unsigned ntrgs;
            fs >> ch >> ntrgs >> ntrgs >> ntrgs; // This number is repeated 3 times

            mesh.triangles().reserve(ntrgs);
            for (unsigned i=0; i<ntrgs; ++i) {
                TriangleIndices t;
                fs >> t[0] >> t[1] >> t[2];
                mesh.add_triangle(t,indmap);
            }
        }

        void save(const OpenMEEG::Mesh& mesh,std::ostream& os) const override {
            os << "- " << mesh.vertices().size() << std::endl;

            const VertexIndices& vertex_index(mesh);
            for (const auto& vertex : mesh.vertices())
                os << *vertex << " " << mesh.normal(*vertex) << std::endl;

            const unsigned ntriangles = mesh.triangles().size();
            os << "- " << ntriangles << ' ' << ntriangles << ' ' << ntriangles << std::endl;
            for (const auto& triangle : mesh.triangles())
                os << vertex_index(triangle,0) << ' '
                   << vertex_index(triangle,1) << ' '
                   << vertex_index(triangle,2) << std::endl;
        }

    private:

        MeshIO* clone(const std::string& filename) const override { return new Tri(filename); }

        Tri(const std::string& filename=""): base(filename,"tri") { }

        static const Tri prototype;

        const char* name() const override { return "TRI"; }
    };
}

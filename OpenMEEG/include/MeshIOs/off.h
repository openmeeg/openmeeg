// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <iostream>
#include <fstream>

#include <map>
#include <string>

#include <om_utils.h>
#include <MeshIO.h>
#include <GeometryExceptions.H>

namespace OpenMEEG::MeshIOs {

    /// \brief Mesh io for TRI file format.

    class OPENMEEG_EXPORT Off: public MeshIO {

        typedef MeshIO base;

    public:

        void load_points(Geometry& geom) override {
            std::string magic;
            fs >> magic;
            if (magic!="OFF")
                OpenMEEG::WrongFileFormat("File is not in OFF format.");

            unsigned npts;
            fs >> io_utils::skip_comments("#") >> npts;
            
            unsigned trash;
            fs >> ntriangles >> trash;
            
            Vertices vertices;
            for (unsigned i=0; i<npts; ++i) {
                Vertex v;
                fs >> v;
                vertices.push_back(v);
            }
            indmap = geom.add_vertices(vertices);
        }

        void load_triangles(OpenMEEG::Mesh& mesh) override {
            reference_vertices(mesh);

            mesh.triangles().reserve(ntriangles);
            for (unsigned i=0; i<ntriangles; ++i) {
                unsigned trash;
                TriangleIndices t;
                fs >> trash >> t[0] >> t[1] >> t[2];
                mesh.add_triangle(t,indmap);
            }
        }

        void save(const OpenMEEG::Mesh& mesh,std::ostream& os) const override {
            os << "OFF" << std::endl;
            os << mesh.vertices().size() << " " << mesh.triangles().size() << " 0" << std::endl;
            const VertexIndices& vertex_index(mesh);

            for (const auto& vertex : mesh.vertices())
                os << *vertex << std::endl;

            for (const auto& triangle : mesh.triangles())
                os << "3 " << vertex_index(triangle,0) << ' '
                           << vertex_index(triangle,1) << ' '
                           << vertex_index(triangle,2) << std::endl;
        }

        MeshIO* clone(const std::string& filename) const override { return new Off(filename); }

    private:

        Off(const std::string& filename=""): base(filename,"off") { }

        static const Off prototype;

        const char* name() const override { return "OFF"; }

        unsigned ntriangles;
    };
}

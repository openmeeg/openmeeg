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

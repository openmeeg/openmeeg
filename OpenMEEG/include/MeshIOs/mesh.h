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

namespace OpenMEEG::MeshIOs {

    /// \brief Mesh io for TRI file format.

    class OPENMEEG_EXPORT Mesh: public MeshIO {

        typedef MeshIO base;

    public:

        void load_points(Geometry& geom) override {

            unsigned char uc[5];
            fs.read(reinterpret_cast<char*>(uc),5); // File format
            fs.read(reinterpret_cast<char*>(uc),4); // lbindian
            //  TODO: we should check that these values are correct.

            unsigned arg_size;
            fs.read(reinterpret_cast<char*>(&arg_size),sizeof(unsigned));
            fs.ignore(arg_size);

            unsigned vertex_per_face;
            fs.read(reinterpret_cast<char*>(&vertex_per_face),sizeof(unsigned));

            unsigned mesh_time;
            fs.read(reinterpret_cast<char*>(&mesh_time),sizeof(unsigned));
            fs.ignore(sizeof(unsigned)); // mesh_step
            unsigned npts;
            fs.read(reinterpret_cast<char*>(&npts),sizeof(unsigned));
            
            // Support only for triangulations and one time frame.

            if (vertex_per_face!=3)
                throw std::invalid_argument("OpenMEEG only handles 3D surfacic meshes.");

            if (mesh_time!=1)
                throw std::invalid_argument("OpenMEEG only handles 3D surfacic meshes with one time frame.");

            float* coords = new float[3*npts]; // Point coordinates
            fs.read(reinterpret_cast<char*>(coords),3*npts*sizeof(float));
            Vertices vertices;
            for (unsigned i=0,j=0; i<npts; ++i,j+=3)
                vertices.push_back(Vertex(coords[j],coords[j+1],coords[j+2]));
            indmap = geom.add_vertices(vertices);
            delete[] coords;

            fs.ignore(sizeof(unsigned));
            fs.ignore(3*npts*sizeof(float)); // Ignore normals.
            fs.ignore(sizeof(unsigned));
        }

        void load_triangles(OpenMEEG::Mesh& mesh) override {
            reference_vertices(mesh);

            unsigned ntrgs; // Number of faces
            fs.read(reinterpret_cast<char*>(&ntrgs),sizeof(unsigned));

            unsigned* pts_inds = new unsigned[3*ntrgs]; // Faces
            fs.read(reinterpret_cast<char*>(pts_inds),3*ntrgs*sizeof(unsigned));
            mesh.triangles().reserve(ntrgs);
            for (unsigned i=0,j=0; i<ntrgs; ++i,j+=3) {
                const TriangleIndices t = { pts_inds[j], pts_inds[j+1], pts_inds[j+2] };
                mesh.add_triangle(t,indmap);
            }
            delete[] pts_inds;
        }

        void save(const OpenMEEG::Mesh& mesh,std::ostream& os) const override {
            unsigned char format[5] = {'b', 'i', 'n', 'a', 'r'}; // File format
            os.write(reinterpret_cast<char*>(format),5);

            unsigned char lbindian[4] = {'D', 'C', 'B', 'A'}; // lbindian
            os.write(reinterpret_cast<char*>(lbindian),4);

            unsigned arg_size = 4;
            os.write(reinterpret_cast<char*>(&arg_size),sizeof(unsigned));

            unsigned char VOID[4] = {'V', 'O', 'I', 'D'}; // Trash
            os.write(reinterpret_cast<char*>(VOID),4);

            unsigned vertex_per_face = 3;
            os.write(reinterpret_cast<char*>(&vertex_per_face),sizeof(unsigned));

            unsigned mesh_time = 1;
            os.write(reinterpret_cast<char*>(&mesh_time),sizeof(unsigned));

            unsigned mesh_step = 0;
            os.write(reinterpret_cast<char*>(&mesh_step),sizeof(unsigned));

            //  Vertices

            float* pts_raw     = new float[mesh.vertices().size()*3]; // Points
            float* normals_raw = new float[mesh.vertices().size()*3]; // Normals

            const VertexIndices& vertex_index(mesh);

            unsigned i = 0;
            for (const auto& vertex : mesh.vertices()) {
                pts_raw[i*3+0]     = static_cast<float>(vertex->x());
                pts_raw[i*3+1]     = static_cast<float>(vertex->y());
                pts_raw[i*3+2]     = static_cast<float>(vertex->z());
                const Normal& n = mesh.normal(*vertex);
                normals_raw[i*3+0] = static_cast<float>(n.x());
                normals_raw[i*3+1] = static_cast<float>(n.y());
                normals_raw[i*3+2] = static_cast<float>(n.z());
                ++i;
            }
            unsigned vertex_number = mesh.vertices().size();
            os.write(reinterpret_cast<char*>(&vertex_number),sizeof(unsigned));
            os.write(reinterpret_cast<char*>(pts_raw),sizeof(float)*vertex_number*3);
            os.write(reinterpret_cast<char*>(&vertex_number),sizeof(unsigned));
            os.write(reinterpret_cast<char*>(normals_raw),sizeof(float)*vertex_number*3);

            delete[] normals_raw;
            delete[] pts_raw;

            //  Triangles

            unsigned* faces_raw = new unsigned[mesh.triangles().size()*3]; // Faces
            i = 0;
            for (const auto& triangle : mesh.triangles()) {
                faces_raw[i*3+0] = vertex_index(triangle,0);
                faces_raw[i*3+1] = vertex_index(triangle,1);
                faces_raw[i*3+2] = vertex_index(triangle,2);
                ++i;
            }

            unsigned char zero = 0;
            os.write(reinterpret_cast<char*>(&zero),1);
            unsigned ntrgs = mesh.triangles().size();
            os.write(reinterpret_cast<char*>(&ntrgs),sizeof(unsigned));
            os.write(reinterpret_cast<char*>(faces_raw),sizeof(unsigned)*ntrgs*3);

            delete[] faces_raw;
        }

        MeshIO* clone(const std::string& filename) const override { return new Mesh(filename); }

    private:

        bool binary() const override { return true; }

        Mesh(const std::string& filename=""): base(filename,"mesh") { }

        static const Mesh prototype;

        const char* name() const override { return "MESH"; }
    };
}

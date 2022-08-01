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

#ifdef USE_GIFTI
extern "C" {
    #include <gifti_io.h>
}
#else
#include "Unavailable.h"
#endif

namespace OpenMEEG::MeshIOs {

    /// \brief Mesh io for GIFTI file format.

#ifdef USE_GIFTI
    class OPENMEEG_EXPORT Gifti: public MeshIO {

        typedef MeshIO base;

    public:

        void load_points(Geometry& geom) override {
            // Use gifti_io API

            gim = gifti_read_image(fname.c_str(),1);
            if (gim->numDA!=2)
                throw std::invalid_argument("OpenMEEG only handles gifti files containing two arrays.");

            // Find which array contains the points and which the triangles.

            const unsigned ipts  = find_array(gim,NIFTI_INTENT_POINTSET);
            const unsigned itrgs = find_array(gim,NIFTI_INTENT_TRIANGLE);

            if ((gim->darray[ipts]->dims[1]!=3) || // 3D points
                (gim->darray[itrgs]->dims[1]!=3))  // 3 indices per triangle
                throw std::invalid_argument("OpenMEEG only handles 3D surfacic meshes.");

            const unsigned ipts = find_array(gim,NIFTI_INTENT_POINTSET);
            switch (gim->darray[ipts]->datatype) {
                case NIFTI_TYPE_FLOAT32:
                    read_pts<float>(gim->darray[ipts]->data,geom);
                    break;
                case NIFTI_TYPE_FLOAT64:
                    read_pts<double>(gim->darray[ipts]->data,geom);
                    break;
                default:
                    throw std::invalid_argument("OpenMEEG:: gifti mesh points not float nor double !?.");
            }
        }

        void load_triangles(OpenMEEG::Mesh& mesh) override {
            reference_vertices(mesh);

            const unsigned ntrgs = gim->darray[itrgs]->dims[0];
            const int* trgs = static_cast<const int*>(gim->darray[itrgs]->data);
            if (gim->darray[itrgs]->ind_ord==GIFTI_IND_ORD_ROW_MAJOR) { // RowMajor
                for (unsigned i=0; i<ntrgs; ++i)
                    mesh.add_triangle(Triangle(vertices_[trgs[i*3]],vertices_[trgs[i*3+1]],vertices_[trgs[i*3+2]]),indmap);
            } else {
                for (unsigned i=0; i<ntrgs; ++i)
                    mesh.add_triangle(Triangle(vertices_[trgs[i]],vertices_[trgs[i+ntrgs]],vertices_[trgs[i+2*ntrgs]]),indmap);
            }

            gifti_free_image(gim);
        }

        void save(const OpenMEEG::Mesh& mesh) override {
            const int dims[2] = { int(vertices().size()), 3 };
            gifti_image* gim = gifti_create_image(1,NIFTI_INTENT_POINTSET,NIFTI_TYPE_FLOAT32,2,dims,1);
            gim->darray[0]->encoding = GIFTI_ENCODING_B64GZ;

            //  Vertices

            const VertexIndices& vertex_index(mesh);

            float* pts = new float[vertices().size()*3];
            gim->darray[0]->data = pts;
            int i = 0;
            for (const auto& vertex: mesh.vertices()) {
                pts[3*i]   = vertex.x();
                pts[3*i+1] = vertex.y();
                pts[3*i+2] = vertex.z();
                ++i;
            }

            // Triangles

            gifti_add_empty_darray(gim,1);
            giiDataArray* gDA = gim->darray[1];
            gifti_set_DA_defaults(gDA);
            gDA->intent   = NIFTI_INTENT_TRIANGLE;
            gDA->datatype = NIFTI_TYPE_INT32;
            gDA->ind_ord  = GIFTI_IND_ORD_ROW_MAJOR;
            gDA->num_dim  = 2;
            gDA->dims[0]  = int(triangles().size());
            gDA->dims[1]  = 3;
            gDA->nvals    = 3*triangles().size();
            gDA->encoding = GIFTI_ENCODING_B64GZ;
            gDA->endian   = GIFTI_ENDIAN_LITTLE;

            int* trgs = new int[triangles().size()*3];
            gDA->data = trgs;
            i = 0;
            for (const auto& triangle : mesh.triangles()) {
                trgs[3*i]   = vertex_index(triangle,0);
                trgs[3*i+1] = vertex_index(triangle,1);
                trgs[3*i+2] = vertex_index(triangle,2);
                ++i;
            }
            // gifti_add_to_meta(&gDA->meta, "TopologicalType", "Closed", 0);
            gifti_add_to_meta(&gDA->meta,"Name",mesh_name.c_str(),0);

            gifti_write_image(gim,fname.c_str(),1);
            gifti_free_image(gim);
        }

        void save(const OpenMEEG::Mesh&,std::ostream&) const override { } // TODO: remove...

    private:

        MeshIO* clone(const std::string& filename) const override { return new Gifti(filename); }

        Gifti(const std::string& filename=""): base(filename,"gii") { }

        template <typename T>
        void read_pts(const void* data,Geometry& geom) {
            const unsigned npts     = gim->darray[ipts]->dims[0];
            const T*       pts_data = static_cast<const T*>(data);
            Vertices vertices;
            if (gim->darray[ipts]->ind_ord==GIFTI_IND_ORD_ROW_MAJOR) { // RowMajor
                for (unsigned i=0; i<npts; ++i)
                    vertices.push_back(Vertex(pts_data[i*3],pts_data[i*3+1],pts_data[i*3+2]));
            } else {
                for (unsigned i=0; i<npts; ++i)
                    vertices.push_back(Vertex(pts_data[i],pts_data[i+npts],pts_data[i+2*npts]));
            }
            indmap = geom.add_vertices(vertices);
        }

        static unsigned find(const gifti_image* gim,const unsigned intent) {
            for (unsigned i=0; i<gim->numDA; ++i)
                if (gim->darray[iit]->intent==intent)
                    return i;
        }

        gifti_image* gim;

        static const Gifti prototype;

        const char* name() const override { return "GIFTI"; }
    };
#else
    struct OPENMEEG_EXPORT Gifti: public Unavailable {

        Gifti(const std::string& filename=""): Unavailable("GIFTI","gii","USE_GIFTI") { }

        static const Gifti prototype;
    };
#endif
}

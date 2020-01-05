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

#include <fstream>
#include <map>
#include <string>

#include <triangle.h>
#include <mesh.h>
#include <geometry.h>
#include <filenames.h>

namespace OpenMEEG {

    //  Mesh class
    //  \brief Mesh is a collection of triangles associated to a geometry containing the points
    //  on which triangles are based.

    class OPENMEEG_EXPORT MeshIO {
    public:

        virtual ~MeshIO() { }

        static MeshIO* create(const std::string& filename) {
            const std::string& extension = tolower(getFilenameExtension(filename));
            return registery.at(extension)->clone(filename);
        }

        virtual const char* name() const = 0;

        void open(const std::ios_base::openmode mode=std::ios_base::in) {
            fs.open(fname,(binary()) ? mode|std::ios_base::binary : mode);
            if (!fs.is_open()) {
                std::ostringstream ost;
                ost << "Error opening " << name() << " file: " << fname << " for reading." << std::endl;
                throw std::invalid_argument(ost.str());
            }
        }

        virtual void load_points(Geometry& geom) = 0;
        virtual void load_triangles(Mesh& m)     = 0;

        virtual void load(Mesh& m) {
            open(std::ios_base::in);
            load_points(m.geometry());
            // TODO
            load_triangles(m);
            fs.close();
            m.update(true);
        }

        virtual void save(const Mesh& mesh,std::ostream& os) const = 0;

        virtual void save(const Mesh& mesh) {
            open(std::ios_base::out);
            save(mesh,fs);
            fs.close();
        }

    protected:

        typedef std::map<std::string,MeshIO*> Registery;

        class VertexIndices {
        public:

            VertexIndices(const Mesh& mesh) {
                unsigned i = 0;
                for (const auto& vertex : mesh.vertices())
                    vmap[vertex] = i++;
            }

            unsigned operator()(const Triangle& triangle,const unsigned ind) const {
                return vmap.at(&(triangle.vertex(ind)));
            }

        private:

            std::map<const Vertex*,unsigned> vmap;
        };

        virtual MeshIO* clone(const std::string& filename) const = 0;
        virtual bool binary() const { return false; }

        void reference_vertices(Mesh& mesh) const { mesh.reference_vertices(indmap); }

        static Registery registery;

        MeshIO(const std::string& filename,const char* name): fname(filename) { registery.insert({ name, this }); }

        std::string  fname;
        std::fstream fs;
        Mesh*        mesh;
        IndexMap     indmap;
    };
}

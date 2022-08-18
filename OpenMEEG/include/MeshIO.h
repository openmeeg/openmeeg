// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <fstream>
#include <map>
#include <string>

#include <triangle.h>
#include <mesh.h>
#include <geometry.h>
#include <filenames.h>
#include <OMExceptions.H>

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
            if (!fs.is_open())
                throw OpenMEEG::OpenError(name());
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
        IndexMap     indmap;
    };
}

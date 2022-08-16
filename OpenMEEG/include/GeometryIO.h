// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <map>
#include <string>

#include <geometry.h>
#include <matrix.h>
#include <filenames.h>

namespace OpenMEEG {

    //  Geometry IO class

    class OPENMEEG_EXPORT GeometryIO {
    public:

        static GeometryIO* create(const std::string& filename) {
            const std::string& extension = tolower(getFilenameExtension(filename));
            try {
                return registery.at(extension)->clone(filename);
            } catch(std::out_of_range&) {
                 throw UnknownFileSuffix(extension);
            }
        }

        virtual const char* name() const = 0;

        void load(Geometry& geometry) {
            load_meshes(geometry);
            load_domains(geometry);
        }

        void load(Geometry& geometry,Matrix& matrix) {
            load(geometry);
            matrix = load_data();
        }

        void save(const Geometry& geometry) {
            save_geom(geometry);
            write();
        }

        void save(const Geometry& geometry,const Matrix& matrix) {
            save_geom(geometry);
            save_data(geometry,matrix);
            write();
        }

        virtual ~GeometryIO() = default;

    protected:

        virtual void   load_meshes(Geometry& geometry) = 0;
        virtual void   load_domains(Geometry&) { }
        virtual Matrix load_data() const = 0;

        virtual void save_geom(const Geometry& geometry) = 0;
        virtual void save_data(const Geometry& geometry,const Matrix& matrix) const = 0;
        virtual void write() const = 0;

        virtual GeometryIO* clone(const std::string& filename) const = 0;

        typedef std::map<std::string,GeometryIO*> Registery;

        static Registery registery;

        GeometryIO(const std::string& filename,const char* name): fname(filename) { registery.insert({ name, this }); }

        std::string fname;
    };
}

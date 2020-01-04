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
            return registery.at(extension)->clone(filename);
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

    protected:

        virtual void   load_meshes(Geometry& geometry) = 0;
        virtual void   load_domains(Geometry& geometry) { }
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

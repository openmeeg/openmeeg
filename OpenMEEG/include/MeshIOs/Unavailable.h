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

#ifdef USE_GIFTI
extern "C" {
    #include <gifti_io.h>
}
#endif

namespace OpenMEEG::MeshIOs {

    /// \brief Not existent IO.

    class OPENMEEG_EXPORT Unavailable: public MeshIO {

        typedef MeshIO base;

    public:

        void load_points(Geometry&)          override { }
        void load_triangles(OpenMEEG::Mesh&) override { }
        void load(OpenMEEG::Mesh&)           override { }

        void save(const OpenMEEG::Mesh&) override { }
        void save(const OpenMEEG::Mesh&,std::ostream&) const override { } // TODO: remove...

        const char* name() const override { return ioname; }

    protected:

        Unavailable(const char* name,const char* extension,const char* cmake):
            base("",extension),ioname(name),cmakevar(cmake)
        { }

    private:

        MeshIO* clone(const std::string&) const override {
            std::cerr << "OpenMEEG not compiled with " << name() << " support. Specify " << cmakevar << " in cmake." << std::endl;
            return const_cast<MeshIO*>(static_cast<const base*>(this)); 
        }

        const char* ioname;
        const char* cmakevar;
    };
}

// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <vector>
#include <vect3.h>

namespace OpenMEEG {

    /// \brief  Vertex
    ///
    ///   Vertex Class
    ///   derived from a Vect3 Class, has an index

    class OPENMEEG_EXPORT Vertex: public Vect3 {
    public:
        #ifdef _MSC_VER
        #pragma warning( suppress : 4245 )  // 'initializing': conversion from 'int' to 'unsigned int', signed/unsigned mismatch
        Vertex(): ind(-1) { };
        #else
        Vertex(): ind(-1) { };
        #endif

        Vertex(const Vect3& V,const unsigned id=-1): Vect3(V),ind(id) { }
        Vertex(const double V[3],const unsigned id=-1): Vect3(V[0],V[1],V[2]),ind(id) { }
        Vertex(const double& x,const double& y,const double& z,const unsigned id=-1): Vect3(x,y,z),ind(id) { }

              unsigned& index()       { return ind; }
        const unsigned& index() const { return ind; }

    private:

        unsigned ind; ///< Index of the vertex
    };

    typedef std::vector<Vertex>  Vertices;
    typedef std::vector<Vertex*> VerticesRefs;
}

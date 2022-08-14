// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <vector>
#include <vertex.h>

namespace OpenMEEG {

    /// \brief Edge
    /// Edge class

    class OPENMEEG_EXPORT Edge {
    public:

        Edge() {
            vertices[0] = nullptr;
            vertices[1] = nullptr;
        }

        Edge(const Vertex& V1,const Vertex& V2) {
            vertices[0] = &V1;
            vertices[1] = &V2;
        }

        const Vertex& vertex(const unsigned i) const { return *vertices[i]; }

        bool operator==(const Edge& e) const { return (e.vertex(0)==vertex(0)) && (e.vertex(1)==vertex(1)); }
        bool operator!=(const Edge& e) const { return (e.vertex(0)!=vertex(0)) || (e.vertex(1)!=vertex(1)); }

    private:

        const Vertex* vertices[2];
    };

    typedef std::vector<Edge> Edges;
}

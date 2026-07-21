// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <vector>
#include <map>

#include <vect3.h>
#include <vertex.h>
#include <edge.h>
#include <OMExceptions.H>

namespace OpenMEEG {

    struct TriangleIndices {

        TriangleIndices() { }

        TriangleIndices(const unsigned i,const unsigned j,const unsigned k) {
            indices[0] = i;
            indices[1] = j;
            indices[2] = k;
        }

        TriangleIndices(const unsigned ind[3]): TriangleIndices(ind[0],ind[1],ind[2]) { }

        TriangleIndices(const TriangleIndices& ind) { std::copy(&ind[0],&ind[3],&indices[0]); }

        TriangleIndices& operator=(const TriangleIndices&) = default;

              unsigned& operator[](const unsigned i)       { return indices[i]; }
        const unsigned& operator[](const unsigned i) const { return indices[i]; }

        unsigned indices[3];
    };

    /// \brief  Triangle
    /// Triangle class

    class OPENMEEG_EXPORT Triangle {
    public:

        typedef       Vertex**       iterator;
        typedef const Vertex** const_iterator;

        /// Constructors

        Triangle() { } // Is this needed besides for SWIG ?

        /// Create a new triangle from 3 vertex adresses.

        Triangle(const Vertex* p1,const Vertex* p2,const Vertex* p3,const unsigned index=-1): vertices{p1,p2,p3},ind(index) {
            out_normal = crossprod(vertex(1)-vertex(0),vertex(2)-vertex(0));
            const double nrm = out_normal.norm();
            surface = 0.5*nrm;
            out_normal /= nrm;
        }

        /// Create a new triangle from a list of 3 vertices.

        Triangle(const Vertex* pts[3],const unsigned index=-1): Triangle(pts[0],pts[1],pts[2],index) { }

        /// Create a new triangle from 3 vertices.

        Triangle(const Vertex& p1,const Vertex& p2,const Vertex& p3,const unsigned index=-1): Triangle(&p1,&p2,&p3,index) { }

        /// Iterators.

        const_iterator begin() const { return const_iterator(vertices);   }
        const_iterator end()   const { return const_iterator(vertices+3); }
        iterator       begin()       { return iterator(vertices);         }
        iterator       end()         { return iterator(vertices+3);       }

        /// Operators

        bool operator==(const Triangle& T) const {
            return (&T.vertex(0)==&vertex(0)) && (&T.vertex(1)==&vertex(1)) && (&T.vertex(2)==&vertex(2));
        }

        const Vertex& vertex(const unsigned& vindex) const { return *vertices[vindex]; }

        Edge edge(const Vertex& V) const {
            const unsigned indx = vertex_index(V);
            return Edge(vertex(indices[indx][0]),vertex(indices[indx][1]));
        }

        Edges edges() const {
            return { Edge(vertex(1),vertex(2)), Edge(vertex(2),vertex(0)), Edge(vertex(0),vertex(1)) };
        }

        const Normal& normal() const { return out_normal; }

        double  area() const { return surface; }

        unsigned& index()       { return ind; }
        unsigned  index() const { return ind; }

        Vect3 center() const { return (vertex(0)+vertex(1)+vertex(2))/3; }

        bool contains(const Vertex& p) const {
            for (unsigned i=0;i<3;++i)
                if (&vertex(i)==&p )
                    return true;
            return false;
        }

        /// Change triangle orientation by flipping two of the vertices.

        void change_orientation() {
            std::swap(vertices[0],vertices[1]);
            out_normal = -out_normal;
        }

        /// Check for intersection with another triangle.

        bool intersects(const Triangle& triangle) const;

    private:

        unsigned vertex_index(const Vertex& V) const {
            for (unsigned i=0;i<3;++i)
                if (&vertex(i)==&V)
                    return i;
            std::ostringstream oss;
            oss << "The vertex with address " << static_cast<const void*>(&V) << " with coordinates " << V
                << " does not belong to the triangle ";
            if (ind!=static_cast<unsigned>(-1))
                oss << ind;
            else
                oss << static_cast<const void*>(this);
            throw UnknownVertex(oss.str());
        }

        static constexpr unsigned indices[3][2] = {{1,2},{2,0},{0,1}};

        const Vertex* vertices[3];   ///< &Vertex-triplet defining the triangle
        double        surface = 0.0; ///< Area
        Normal        out_normal;    ///< Normal
        unsigned      ind;           ///< Index of the triangle
    };

    inline double
    Vect3::solid_angle(const Triangle& T) const { return solid_angle(T.vertex(0),T.vertex(1),T.vertex(2)); }

    typedef std::vector<Triangle>  Triangles;
    typedef std::vector<Triangle*> TrianglesRefs;

    typedef std::map<unsigned,unsigned> IndexMap;
}

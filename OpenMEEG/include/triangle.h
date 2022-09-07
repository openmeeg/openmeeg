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

        /// Create a new triangle from a set of vertices.

        Triangle(Vertex* pts[3],const unsigned index=-1): ind(index) {
            for (unsigned i=0;i<3;++i)
                vertices_[i] = pts[i];
        }

        /// Create a new triangle from a 3 vertex adresses.

        Triangle(Vertex* p1,Vertex* p2,Vertex* p3,const unsigned index=-1): vertices_{p1,p2,p3},ind(index) { }

        /// Create a new triangle from a 3 vertices.

        Triangle(Vertex& p1,Vertex& p2,Vertex& p3,const unsigned index=-1): Triangle(&p1,&p2,&p3,index) { }

        /// Iterators.

        const_iterator begin() const { return const_iterator(vertices_);   }
        const_iterator end()   const { return const_iterator(vertices_+3); }
        iterator       begin()       { return iterator(vertices_);         }
        iterator       end()         { return iterator(vertices_+3);       }

        /// Operators

        bool operator==(const Triangle& T) const {
            return (&T.vertex(0)==&vertex(0)) && (&T.vertex(1)==&vertex(1)) && (&T.vertex(2)==&vertex(2));
        }

              Vertex& vertex(const unsigned& vindex)       { return *vertices_[vindex]; }
        const Vertex& vertex(const unsigned& vindex) const { return *vertices_[vindex]; }

        Edge edge(const Vertex& V) const {
            const unsigned indx = vertex_index(V);
            return Edge(vertex(indices[indx][0]),vertex(indices[indx][1]));
        }

        Edges edges() const {
            return { Edge(vertex(1),vertex(2)), Edge(vertex(2),vertex(0)), Edge(vertex(0),vertex(1)) };
        }

              Normal&  normal()       { return normal_; }
        const Normal&  normal() const { return normal_; }

        double  area() const { return area_; }
        double& area()       { return area_; }

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

        void change_orientation() { std::swap(vertices_[0],vertices_[1]); }

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
            if (ind!=-1U)
                oss << ind;
            else
                oss << static_cast<const void*>(this);
            throw UnknownVertex(oss.str());
        }

        static constexpr unsigned indices[3][2] = {{1,2},{2,0},{0,1}};

        Vertex*  vertices_[3]; ///< &Vertex-triplet defining the triangle
        double   area_ = 0.0;  ///< Area
        Normal   normal_;      ///< Normal
        unsigned ind;          ///< Index of the triangle
    };

    typedef std::vector<Triangle>  Triangles;
    typedef std::vector<Triangle*> TrianglesRefs;

    typedef std::map<unsigned,unsigned> IndexMap;
}

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

#include <vector>
#include <map>

#include <vect3.h>
#include <vertex.h>
#include <edge.h>
#include <GeometryExceptions.H>

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

        Triangle(): ind(-1) { }

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
            const unsigned ind = vertex_index(V);
            return Edge(vertex(indices[ind][0]),vertex(indices[ind][1]));
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
            throw UnknownVertex();
        }

        static constexpr unsigned indices[3][2] = {{1,2},{2,0},{0,1}};

        Vertex*  vertices_[3]; ///< &Vertex-triplet defining the triangle
        double   area_;        ///< Area
        Normal   normal_;      ///< Normal
        unsigned ind;          ///< Index of the triangle
    };

    typedef std::vector<Triangle>  Triangles;
    typedef std::vector<Triangle*> TrianglesRefs;

    typedef std::map<unsigned,unsigned> IndexMap;
}

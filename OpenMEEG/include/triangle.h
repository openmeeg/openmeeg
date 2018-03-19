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

#include <cstdlib>
#include <vector>
#include <set>
#include <vect3.h>
#include <vertex.h>

namespace OpenMEEG {

    /// \brief  Triangle
    ///
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

        /// Create a new triangle from a 3 vertices.

        Triangle(Vertex& p1,Vertex& p2,Vertex& p3,const unsigned index=-1): ind(index) {
            vertices_[0] = &p1;
            vertices_[1] = &p2;
            vertices_[2] = &p3;
        }

        /// Create a new triangle from a 3 vertex adresses.
        
        Triangle(Vertex* p1,Vertex* p2,Vertex* p3,const unsigned index=-1): ind(index) {
            vertices_[0] = p1;
            vertices_[1] = p2;
            vertices_[2] = p3;
        }
        
        /// Operators
        // Having both [] and () doing different things is prone to errors. Also the %3 in indexing seems vety costly.

              Vertex*   operator[](const unsigned& vindex)       { return vertices_[vindex%3];  } // 0 <= 'index' <= '2'
        const Vertex*   operator[](const unsigned& vindex) const { return vertices_[vindex%3];  }
              Vertex&   operator()(const unsigned& vindex)       { return *vertices_[vindex%3]; } // 0 <= 'index' <= '2'
        const Vertex&   operator()(const unsigned& vindex) const { return *vertices_[vindex%3]; }

              bool      operator==(const Triangle& T)     const { return (T[0]==vertices_[0])&(T[1]==vertices_[1])&(T[2]==vertices_[2]); }
                                                 
              Vertex&        vertex(const unsigned& vindex)       { return operator()(vindex); }
        const Vertex&        vertex(const unsigned& vindex) const { return operator()(vindex); }

        /// Iterators.

        const_iterator begin() const { return const_iterator(vertices_);   }
        const_iterator end()   const { return const_iterator(vertices_+3); }
        iterator       begin()       { return iterator(vertices_);         }
        iterator       end()         { return iterator(vertices_+3);       }

        //  These (s[1-3]) are ugly.. Remove ?

        const Vertex&  s1()    const { return *vertices_[0]; }
        const Vertex&  s2()    const { return *vertices_[1]; }
        const Vertex&  s3()    const { return *vertices_[2]; }

              Vertex&  s1()          { return *vertices_[0]; }
              Vertex&  s2()          { return *vertices_[1]; }
              Vertex&  s3()          { return *vertices_[2]; }

              Normal&  normal()       { return normal_; }
        const Normal&  normal() const { return normal_; }
                                
              double&  area()         { return area_; }
        const double&  area()   const { return area_; }
                                
              unsigned& index()        { return ind; }
        const unsigned& index()  const { return ind; }

        // These two methods are ugly and used exactly once. There is probably a better way.

        const Vertex& prev(const Vertex& V) const { 
            if ( &V == vertices_[0]) {
                return *vertices_[2];
            } else if ( &V == vertices_[1] ) {
                return *vertices_[0];
            } else if ( &V == vertices_[2] ) {
                return *vertices_[1];
            } else {
                static Vertex v;
                return v;
            }
        }
        const Vertex& next(const Vertex& V) const { 
            if ( &V == vertices_[0]) {
                return *vertices_[1];
            } else if ( &V == vertices_[1] ) {
                return *vertices_[2];
            } else if ( &V == vertices_[2] ) {
                return *vertices_[0];
            } else {
                static Vertex v;
                return v;
            }
        }

        Vect3 center() const {
            return (*vertices_[0]+*vertices_[1]+*vertices_[2])/3;
        }

        bool contains(const Vertex& p) const {
            for (unsigned i=0;i<3;++i)
                if (&vertex(i)==&p )
                    return true;
            return false;
        }

        /// flip two of the three vertex address

        void flip() { std::swap(vertices_[0],vertices_[1]); }

    private:

        Vertex*  vertices_[3]; ///< &Vertex-triplet defining the triangle
        double   area_;        ///< Area
        Normal   normal_;      ///< Normal
        unsigned ind;          ///< Index of the triangle
    };

    typedef std::vector<Triangle> Triangles;
}

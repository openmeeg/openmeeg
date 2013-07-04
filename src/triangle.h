/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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

#ifndef OPENMEEG_TRIANGLE_H
#define OPENMEEG_TRIANGLE_H

#include <cstdlib>
#include <vector>
#include <set>
#include "vect3.h"
#include "vertex.h"

namespace OpenMEEG {

    /** \brief  Triangle

        Triangle class

    **/

    //  The default class to handle pointers to a vertex.
    //  This is simply a pointer here, but sometimes we want to attach
    //  some information to this pointer.

    class Reference {

    public:

        Reference() { }

        Reference&    operator=(Reference v) { vref = v.vref; return * this; }
        Reference&    operator=(Vertex* v)   { vref = v;      return * this; }

              Vertex& vertex()               { return *vref; }
        const Vertex& vertex() const         { return *vref; }

              size_t& index()                { return vref->index(); }
        const size_t& index()  const         { return vref->index(); }

    private:

        Vertex* vref;
    };

    class Triangle {

    public:

        typedef       Reference*       iterator;
        typedef const Reference* const_iterator;

        // Constructors
        Triangle(): index_(-1) {}
        Triangle(const Triangle& t); // copy constructor
        Triangle(Vertex *pts[3]); // Create a new face from a set of vertices.
        Triangle(Vertex& p1, Vertex& p2, Vertex& p3);
        Triangle(Vertex * p1, Vertex * p2, Vertex * p3);
        
        // Destructor
        ~Triangle() { destroy(); }

              Triangle&  operator=(const Triangle& t);

        // 0 <= 'index' <= '2'
              Reference& operator()(const size_t& vindex)       { return vertices[vindex]; }
        const Reference& operator()(const size_t& vindex) const { return vertices[vindex]; }
                                                 
              Vertex&        vertex(const size_t& vindex)       { return operator()(vindex).vertex(); }
        const Vertex&        vertex(const size_t& vindex) const { return operator()(vindex).vertex(); }

        // Iterators.
              iterator       begin()                     { return iterator(vertices);       }
              const_iterator begin()               const { return const_iterator(vertices); }
              iterator       end()                       { return iterator(vertices+3);       }
              const_iterator end()                 const { return const_iterator(vertices+3); }

        const Vertex&        prev(const Vertex& V) const { return (s1().vertex() == V)?s3().vertex():(s2().vertex() == V)?s1().vertex():s2().vertex(); }
        const Vertex&        next(const Vertex& V) const { return (s1().vertex() == V)?s2().vertex():(s2().vertex() == V)?s3().vertex():s1().vertex(); }

              Reference&     s1()                        { return vertices[0]; }
              Reference&     s2()                        { return vertices[1]; }
              Reference&     s3()                        { return vertices[2]; }
                                 
        const Reference&     s1()                  const { return vertices[0]; }
        const Reference&     s2()                  const { return vertices[1]; }
        const Reference&     s3()                  const { return vertices[2]; }

              Normal&        normal()                    { return normal_; }
        const Normal&        normal()              const { return normal_; }
                                     
              double&        area()                      { return area_; }
        const double&        area()                const { return area_; }
                                     
              size_t&        index()                     { return index_; }
        const size_t&        index()               const { return index_; }

        bool contains(const Vertex& p) const {
            for (size_t i = 0; i < 3; i++) {
                if (&vertex(i) == &p) {
                    return true;
                }
            }
            return false;
        }

        bool operator==(const Triangle& T) const;

    private:

        void copy(const Triangle &t);
        void destroy();

        Reference vertices[3]; // &Vertex-triplet defining the triangle
        double    area_;       // Area
        Normal    normal_;     // Normal
        size_t    index_;      // Index of the triangle
    };

    typedef std::vector<Triangle>       Triangles;
    typedef std::set<const Triangle *>  SetPTriangle;
}

#endif  //! OPENMEEG_TRIANGLE_H

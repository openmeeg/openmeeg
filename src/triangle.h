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

        Triangle() {}

        //  Create a new face from a set of vertices.
        Triangle(Vertex *pts[3]): _index(-1) {
            Triangle& f = *this;
            for (unsigned i = 0; i < 3; i++) {
                f(i) = pts[i];
            }
        }

        Triangle(Vertex& p1, Vertex& p2, Vertex& p3): _index(-1) {
            Triangle& f = *this;
            f(0) = &p1;
            f(1) = &p2;
            f(2) = &p3;
        }

        //  0 <= 'index' <= '2'
              Reference& operator()(const size_t vindex)       { return vertices[vindex]; }
        const Reference& operator()(const size_t vindex) const { return vertices[vindex]; }

              Vertex&        vertex(const size_t vindex)       { return operator()(vindex).vertex(); }
        const Vertex&        vertex(const size_t vindex) const { return operator()(vindex).vertex(); }

        //  Iterators.
              iterator       begin()                           { return iterator(vertices);       }
              const_iterator begin()                     const { return const_iterator(vertices); }
              iterator       end()                             { return iterator(vertices+3);       }
              const_iterator end()                       const { return const_iterator(vertices+3); }

        const Vertex&        next(const size_t i)        const { return operator()((1+i)%3).vertex(); }
        const Vertex&        prev(const size_t i)        const { return operator()((i-1)%3).vertex(); }
        const Vertex&        next(const Vertex V)        const { return next(vertex_index(V)); }
        const Vertex&        prev(const Vertex V)        const { return prev(vertex_index(V)); }

              Reference&     s1()                              { return vertices[0]; }
              Reference&     s2()                              { return vertices[1]; }
              Reference&     s3()                              { return vertices[2]; }
                                 
              Reference      s1()                        const { return vertices[0]; }
              Reference      s2()                        const { return vertices[1]; }
              Reference      s3()                        const { return vertices[2]; }

              Vect3&         normal()                          { return _normal; }
        const Vect3&         normal()                    const { return _normal; }
                                     
              double&        area()                            { return _area; }
        const double&        area()                      const { return _area; }
                                     
              size_t&        index()                           { return _index; }
        const size_t&        index()                     const { return _index; }

        bool contains(const Vertex& p) const {
            for (Triangle::const_iterator tit = this->begin(); tit != this->end(); tit++) {
                if (tit->vertex() == p) {
                    return true;
                }
            }
            return false;
        }

    private:

        Reference vertices[3]; // Vertex-triplet defining the triangle
        double    _area;       // Area
        Vect3     _normal;     // Normal
        size_t    _index;      // Index of the triangle

        size_t vertex_index(const Vertex& p) const { // returns 0, 1, or 2 (or the size_t MAX VALUE)
            size_t vindex = 0;
            for (Triangle::const_iterator tit = this->begin(); tit != this->end(); tit++, vindex++) {
                if (tit->vertex() == p) {
                    return vindex;
                }
            }
            return -1;
        }
    };

    bool operator==(const Triangle& T1, const Triangle& T2) {
        for (Triangle::const_iterator i1 = T1.begin(), i2 = T2.begin(); i1 != T1.end(); ++i1, ++i2) {
            if (&*i1 != &*i2) {
                return false;
            }
        }
        return true;
    }

    bool operator<(const Triangle& T1, const Triangle& T2) {
        for (Triangle::const_iterator i1 = T1.begin(), i2 = T2.begin(); i1 != T1.end(); ++i1, ++i2) {
            if (i1->vertex() < i2->vertex()) {
                return true;
            }
        }
        return false;
    }

    typedef std::vector<Triangle> Triangles;
    typedef std::set<Triangle>    SetTriangle;
}

#endif  //! OPENMEEG_TRIANGLE_H

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
#include "vect3.h"
#include "point.h"


namespace OpenMEEG {

    /** \brief  Triangle

        Triangle class

    **/

    //  The default class to handle pointers to a point.
    //  This is simply a pointer here, but sometimes we want to attach
    //  some information to this pointer, this is why this class can be
    //  user defined in the Triangle structure.
    class Reference {
    public:

        Reference() { }

        Reference& operator=(Reference v) { vref = v.vref; return * this; }
        Reference& operator=(Point* v)    { vref = v;      return * this; }

              Point& point()       { return *vref; }
        const Point& point() const { return *vref; }

    private:

        Point* vref;
    };

    class Triangle {

    public:

        typedef unsigned         Index;
        typedef Reference*       iterator;
        typedef const Reference* const_iterator;

        //  Create a new face from a set of points.
        Triangle(Point *pts[3]) {
            Triangle& f = *this;
            for (unsigned i = 0; i < 3; i++) {
                f(i) = pts[i];
            }
        }

        Triangle(Point& p1, Point& p2, Point& p3) {
            Triangle& f = *this;
            f(0) = &p1;
            f(1) = &p2;
            f(2) = &p3;
        }

        //  0 <= 'index' <= '2'
              Reference& operator()(const Index index)       { return points[index]; }
        const Reference& operator()(const Index index) const { return points[index]; }

              Point& point(const Index index)       { return operator()(index).point(); }
        const Point& point(const Index index) const { return operator()(index).point(); }

        //  Iterator.
              iterator begin()       { return iterator(points);       }
        const_iterator begin() const { return const_iterator(points); }

              iterator end()       { return iterator(points+3);       }
        const_iterator end() const { return const_iterator(points+3); }

        Index index(const iterator i)       const { return i-begin(); }
        Index index(const const_iterator i) const { return i-begin(); }

        Point center() const {
            Point sum(0., 0., 0.);
            for (Triangle::const_iterator v = begin(); v != end(); ++v) {
                sum += v->point();
            }
            return static_cast<Point>(sum/static_cast<double>(3.));
        }

        const Point& next(const Index i) const {
            return operator()((1+i)%3).point();
        }

        const Point& prev(const Index i) const {
            return operator()((i-1)%3).point();
        }

        Reference&   s1()           { return points[0]; }
        Reference&   s2()           { return points[1]; }
        Reference&   s3()           { return points[2]; }
                         
        Reference    s1()     const { return points[0]; }
        Reference    s2()     const { return points[1]; }
        Reference    s3()     const { return points[2]; }

        const Vect3& normal() const { return n; }
        Vect3&       normal()       { return n; }

        bool contains(const Point& p) const {
            for (Triangle::const_iterator tit = this->begin(); tit != this->end(); tit++) {
                if (tit->point() == p) {
                    return true;
                }
            }
            return false;
        }

        inline double  getArea() const   { return _area; };
        inline void    setArea(double a) { _area = a; };
        inline double& area()            { return _area; }
    private:

        Reference points[3];
        double    _area;    // Area
        Vect3     n;        // Normal
    };

    bool operator==(const Triangle& F1, const Triangle& F2) {
        for (Triangle::const_iterator i1 = F1.begin(), i2 = F2.begin(); i1 != F1.end(); ++i1, ++i2) {
            if (&*i1 != &*i2) {
                return false;
            }
        }
        return true;
    }

    // template <typename MESH,typename CONTEXT>
    // std::istream& operator>>(std::istream& is,TriangleReader<MESH,CONTEXT>& reader) {
        // typedef typename MESH::Point Point;
        // typedef typename MESH::Triangle   Triangle;
        // Point *facePoints[Triangle::Dim+1];
        // for (unsigned j=0;j<=Triangle::Dim;++j) {
            // Index pointIndex;
            // is >> pointIndex;
            // Maths::minmax(--pointIndex,reader.min_point,reader.max_point);
            // facePoints[j] = &reader.mesh.point(pointIndex);
        // }
        // reader.face_ptr = new (&*(reader.mesh.faces().end())) Triangle(is,facePoints,reader.context);
        // return is;
    // }


        // inline bool operator!= (const Triangle &t ) const {return (m_s1!=t[0] || m_s2!=t[1] || m_s3!=t[2]);}

        // friend std::istream& operator>>(std::istream &is, Triangle &t);

    // inline std::istream& operator>>(std::istream &is, Triangle &t) {
        // return is >> t.m_s1 >> t.m_s2 >> t.m_s3;
    // }

    // inline std::ostream& operator<<(std::ostream &os, const Triangle &t) {
        // return os << t[0] << " " << t[1] << " " << t[2];
    // }

    // assigment operator
    // inline Triangle& Triangle::operator= (const Triangle &t) {
        // check for self-assignment
        // if (this == &t) {
            // return *this;
        // }

        // do the copy
        // m_s1 = t.m_s1;
        // m_s2 = t.m_s2;
        // m_s3 = t.m_s3;
        // m_area = t.m_area;
        // n = t.n;
        // return the existing object
        // return *this;
    // }

}
#endif  //! OPENMEEG_TRIANGLE_H

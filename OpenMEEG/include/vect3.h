// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#ifdef WIN32
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_DEPRECATE 1
#endif
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <om_common.h>
#include <vector.h>
#include <vector>

#include <OpenMEEG_Export.h>

namespace OpenMEEG {

    inline double sqr(const double x) { return x*x; }

    class Triangle;

    /// \brief  Vect3

    class OPENMEEG_EXPORT Vect3 {

        double m[3]; //!< Coordinates of the vector

    public:

        Vect3(const double x1,const double x2,const double x3) { m[0] = x1; m[1] = x2; m[2] = x3; }
        Vect3(const double a=0.0) { std::fill(&m[0],&m[3],a); }

        Vect3(const Vect3& v) { std::copy(&v.m[0],&v.m[3],&m[0]); }

        operator Vector() const {
            Vector V(3);
            std::copy(&m[0],&m[3],V.data());
            return V;
        }

        operator const double*() const { return m; }

        Vect3& operator=(const double v) {
            for (unsigned i=0;i<3;++i)
                m[i] = v;
            return *this;
        }

        Vect3& operator=(const Vect3& v) {
            std::copy(&v.m[0],&v.m[3],&m[0]);
            return *this;
        }

        double& x()       { return m[0]; }
        double  x() const { return m[0]; }

        double& y()       { return m[1]; }
        double  y() const { return m[1]; }

        double& z()       { return m[2]; }
        double  z() const { return m[2]; }

        double operator<(const Vect3& v) const { return ((m[0]!=v.x()) ? (m[0]<v.x()) : ((m[1]!=v.y()) ? (m[1]<v.y()) : (m[2]<v.z()))); }

        double norm()  const { return sqrt(norm2());                 }
        double norm2() const { return sqr(m[0])+sqr(m[1])+sqr(m[2]); }

        Vect3  unit_vector() const { return *this/norm(); }

        Vect3& normalize() { return *this /= norm(); }

        bool operator==(const Vect3& v) const { return (m[0]==v.x() && m[1]==v.y() && m[2]==v.z()); }
        bool operator!=(const Vect3& v) const { return (m[0]!=v.x() || m[1]!=v.y() || m[2]!=v.z()); }

        Vect3& operator+=(const Vect3& v) { m[0] += v.x(); m[1] += v.y(); m[2] += v.z(); return *this; }
        Vect3& operator-=(const Vect3& v) { m[0] -= v.x(); m[1] -= v.y(); m[2] -= v.z(); return *this; }
        Vect3& operator*=(const double d) { m[0] *= d; m[1] *= d; m[2] *= d; return *this; }
        Vect3& operator/=(const double d) { return operator*=(1.0/d); }

        void multadd(const double d,const Vect3& v) {m[0] += d*v.x(); m[1] += d*v.y(); m[2] += d*v.z();}

        Vect3 operator+(const Vect3& v) const { return Vect3(m[0]+v.x(),m[1]+v.y(),m[2]+v.z()); }
        Vect3 operator-(const Vect3& v) const { return Vect3(m[0]-v.x(),m[1]-v.y(),m[2]-v.z()); }
        Vect3 operator*(const double d) const { return Vect3(d*m[0],d*m[1],d*m[2]); }
        Vect3 operator/(const double d) const { return Vect3(m[0]/d,m[1]/d,m[2]/d); }

        double operator()(const int i) const {
            om_assert(i>=0 && i<3);
            return m[i];
        }

        double& operator()(const int i) {
            om_assert(i>=0 && i<3);
            return m[i];
        }

        Vect3 operator-() const { return Vect3(-m[0],-m[1],-m[2]); }

        friend std::ostream& operator<<(std::ostream& os,const Vect3& v);
        friend std::istream& operator>>(std::istream& is,Vect3& v);
    };

    inline Vect3  operator*(const double d,const Vect3& V)   { return V*d;   }
    inline double dotprod(const Vect3& V1,const Vect3& V2)   { return V1.x()*V2.x()+V1.y()*V2.y()+V1.z()*V2.z(); }
    inline Vect3  crossprod(const Vect3& V1,const Vect3& V2) { return Vect3(V1.y()*V2.z()-V1.z()*V2.y(),V1.z()*V2.x()-V1.x()*V2.z(),V1.x()*V2.y()-V1.y()*V2.x()); }

    inline double det(const Vect3& V1,const Vect3& V2,const Vect3& V3) { return dotprod(V1,crossprod(V2,V3)); }

    inline std::istream& operator>>(std::istream& is,Vect3& v) {
        return is >> v.x() >> v.y() >> v.z();
    }

    inline std::ostream& operator<<(std::ostream& os,const Vect3& v) {
        return os << v.x() << ' ' << v.y() << ' ' << v.z() ;
    }

    typedef Vect3                Normal;
    typedef std::vector<Normal>  Normals;
}

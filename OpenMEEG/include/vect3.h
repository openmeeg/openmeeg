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

#if WIN32
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_DEPRECATE 1
#endif
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <om_common.h>
#include <vector>

#include <OpenMEEG_Export.h>

namespace OpenMEEG {

    inline double sqr(const double& x) { return x*x; }

    /// \brief  Vect3

    class OPENMEEG_EXPORT Vect3 {

        double m[3]; //!< Coordinates of the vector

    public:

        Vect3() { }
        Vect3(const double& x1,const double& x2,const double& x3) { m[0] = x1; m[1] = x2; m[2] = x3; }
        Vect3(const double& a) { std::fill(&m[0],&m[3],a); }
        ~Vect3() { }

        Vect3(const Vect3& v) {
            for (unsigned i=0;i<3;++i)
                m[i] = v.m[i];
        }

        operator const double*() const { return m; }

        Vect3& operator=(const double& v) {
            for (unsigned i=0;i<3;++i)
                m[i] = v;
            return *this;
        }

        Vect3& operator=(const Vect3& v) {
            std::copy(&v.m[0],&v.m[3],&m[0]);
            return *this;
        }

              double& x()       { return m[0]; }
        const double& x() const { return m[0]; }

              double& y()       { return m[1]; }
        const double& y() const { return m[1]; }

              double& z()       { return m[2]; }
        const double& z() const { return m[2]; }

        double operator<(const Vect3& v) const { return ((m[0]!=v.x()) ? (m[0]<v.x()) : ((m[1]!=v.y()) ? (m[1]<v.y()) : (m[2]<v.z()))); }

        double norm()  const { return sqrt(norm2());                 }
        double norm2() const { return sqr(m[0])+sqr(m[1])+sqr(m[2]); }

        bool operator==(const Vect3& v ) const { return (m[0]==v.x() && m[1]==v.y() && m[2]==v.z()); }
        bool operator!=(const Vect3& v ) const { return (m[0]!=v.x() || m[1]!=v.y() || m[2]!=v.z()); }

        void operator+=(const Vect3& v)  { m[0] += v.x(); m[1] += v.y(); m[2] += v.z(); }
        void operator-=(const Vect3& v)  { m[0] -= v.x(); m[1] -= v.y(); m[2] -= v.z(); }
        void operator*=(const double& d) { m[0] *= d; m[1] *= d; m[2] *= d; }
        void operator/=(const double& d) { operator*=(1.0/d); }

        void multadd(const double& d, const Vect3& v) {m[0] += d*v.x(); m[1] += d*v.y(); m[2] += d*v.z();}

        Vect3 operator+(const Vect3& v)  const { return Vect3(m[0]+v.x(),m[1]+v.y(),m[2]+v.z()); }
        Vect3 operator-(const Vect3& v)  const { return Vect3(m[0]-v.x(),m[1]-v.y(),m[2]-v.z()); }
        Vect3 operator^(const Vect3& v)  const { return Vect3(m[1]*v.z()-m[2]*v.y(),m[2]*v.x()-m[0]*v.z(),m[0]*v.y()-m[1]*v.x()); }
        Vect3 operator*(const double& d) const { return Vect3(d*m[0],d*m[1],d*m[2]); }
        Vect3 operator/(const double& d) const { return Vect3(m[0]/d,m[1]/d,m[2]/d); }

        double operator()(const int i) const {
            om_assert(i>=0 && i<3);
            return m[i];
        }

        double& operator()(const int i) {
            om_assert(i>=0 && i<3);
            return m[i];
        }

        Vect3 operator-() const { return Vect3(-m[0],-m[1],-m[2]); }

        inline double solid_angle(const Vect3& v1,const Vect3& v2,const Vect3& v3) const;

        Vect3& normalize() {
            *this /= (*this).norm();
            return *this;
        }

        friend std::ostream& operator<<(std::ostream& os, const Vect3& v);
        friend std::istream& operator>>(std::istream& is, Vect3& v);
    };

    inline Vect3  operator*(const double& d,const Vect3& V)  { return V*d;   }
    inline double dotprod(const Vect3& V1,const Vect3& V2)   { return V1.x()*V2.x()+V1.y()*V2.y()+V1.z()*V2.z(); }
    inline Vect3  crossprod(const Vect3& V1,const Vect3& V2) { return V1^V2; }
    inline double det(const Vect3& V1,const Vect3& V2,const Vect3& V3) { return dotprod(V1,crossprod(V2,V3)); }

    inline double Vect3::solid_angle(const Vect3& V1,const Vect3& V2,const Vect3& V3) const {
        // De Munck : Good sign directly
        const Vect3& V0 = *this;
        const Vect3& Y1 = V1-V0;
        const Vect3& Y2 = V2-V0;
        const Vect3& Y3 = V3-V0;
        const double y1 = Y1.norm();
        const double y2 = Y2.norm();
        const double y3 = Y3.norm();
        const double d = det(Y1,Y2,Y3);
        return (fabs(d)<1e-10) ? 0.0 : 2*atan2(d,(y1*y2*y3+y1*dotprod(Y2,Y3)+y2*dotprod(Y3,Y1)+y3*dotprod(Y1,Y2)));
    }

    inline std::istream& operator>>(std::istream& is, Vect3& v) {
        return is >> v.x() >> v.y() >> v.z();
    }

    inline std::ostream& operator<<(std::ostream& os, const Vect3& v) {
        return os << v.x() << ' ' << v.y() << ' ' << v.z() ;
    }

    typedef Vect3                Normal;
    typedef std::vector<Normal>  Normals;
}

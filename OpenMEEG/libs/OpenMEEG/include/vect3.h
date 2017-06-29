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

    /** \brief  Vect3

        Mesh Class

    **/

    class OPENMEEG_EXPORT Vect3 {

        double m[3]; //!< Coordinates of the vector

    public:

        inline  Vect3() { }
        inline  Vect3(const double& xx, const double& yy, const double& zz) { m[0] = xx; m[1] = yy; m[2] = zz; }
        inline  Vect3(const double& a) { std::fill(&m[0], &m[3], a); }
        inline ~Vect3() { }

        Vect3& operator=(const Vect3& v) {
            std::copy(&v.m[0], &v.m[3], &m[0]);
            return *this;
        }

        Vect3(const Vect3& v) {
            m[0] = v.x();
            m[1] = v.y();
            m[2] = v.z();
        }

        inline       double& x()       { return m[0]; }
        inline const double& x() const { return m[0]; }

        inline       double& y()       { return m[1]; }
        inline const double& y() const { return m[1]; }

        inline       double& z()       { return m[2]; }
        inline const double& z() const { return m[2]; }

        inline double operator*(const Vect3& v) const { return m[0]*v.x()+m[1]*v.y()+m[2]*v.z(); }
        inline double operator<(const Vect3& v) const { return ((m[0] != v.x())?(m[0] < v.x()):((m[1] != v.y())? (m[1] < v.y()):(m[2] < v.z()))); }

        inline double norm()  const { return sqrt(norm2()); }
        inline double norm2() const { return m[0]*m[0]+m[1]*m[1]+m[2]*m[2]; }

        inline bool operator==(const Vect3& v ) const { return (m[0]==v.x() && m[1]==v.y() && m[2]==v.z()); }
        inline bool operator!=(const Vect3& v ) const { return (m[0]!=v.x() || m[1]!=v.y() || m[2]!=v.z()); }

        inline void operator+=(const Vect3& v)  { m[0] += v.x(); m[1] += v.y(); m[2] += v.z(); }
        inline void operator-=(const Vect3& v)  { m[0] -= v.x(); m[1] -= v.y(); m[2] -= v.z(); }
        inline void operator*=(const double& d) { m[0] *= d; m[1] *= d; m[2] *= d; }
        inline void operator/=(const double& d) { operator*=(1.0/d); }

        inline void multadd(const double& d, const Vect3& v) {m[0] += d*v.x(); m[1] += d*v.y(); m[2] += d*v.z();}

        inline Vect3 operator+(const Vect3& v)  const { return Vect3(m[0]+v.x(), m[1]+v.y(), m[2]+v.z()); }
        inline Vect3 operator-(const Vect3& v)  const { return Vect3(m[0]-v.x(), m[1]-v.y(), m[2]-v.z()); }
        inline Vect3 operator^(const Vect3& v)  const { return Vect3(m[1]*v.z()-m[2]*v.y(), m[2]*v.x()-m[0]*v.z(), m[0]*v.y()-m[1]*v.x()); }
        inline Vect3 operator*(const double& d) const { return Vect3(d*m[0], d*m[1], d*m[2]); }
        inline Vect3 operator/(const double& d) const { return Vect3(m[0]/d, m[1]/d, m[2]/d); }

        inline double operator() (const int i) const {
            om_assert(i>=0 && i<3);
            return m[i];
        }

        inline double& operator()(const int i) {
            om_assert(i>=0 && i<3);
            return m[i];
        }

        inline Vect3 operator-() { return Vect3(-m[0], -m[1], -m[2]); }

        inline double det(const Vect3& y2, const Vect3& y3) const {
            return (*this)*(y2^y3); // y1.det(y2, y3):= y1/(y2^y3)
        }

        inline double solangl(const Vect3& v1, const Vect3& v2, const Vect3& v3) const {
            // De Munck : Good sign directly
            const Vect3 Y1 = v1 - *this;
            const Vect3 Y2 = v2 - *this;
            const Vect3 Y3 = v3 - *this;
            const double y1 = Y1.norm();
            const double y2 = Y2.norm();
            const double y3 = Y3.norm();
            const double d = Y1*(Y2^Y3);
            return 2.*atan2(d, (y1*y2*y3+y1*(Y2*Y3)+y2*(Y3*Y1)+y3*(Y1*Y2)));
        }

        inline void normalize() {
            *this /= (*this).norm();
        }

        friend std::ostream& operator<<(std::ostream& os, const Vect3& v);
        friend std::istream& operator>>(std::istream& is, Vect3& v);
    };

    inline Vect3 operator*(const double& d, const Vect3& v) { return v*d; }

    inline std::istream& operator>>(std::istream& is, Vect3& v) {
        return is >> v.x() >> v.y() >> v.z();
    }

    inline std::ostream& operator<<(std::ostream& os, const Vect3& v) {
        return os << v.x() << ' ' << v.y() << ' ' << v.z() ;
    }

    typedef Vect3                Normal;
    typedef std::vector<Normal>  Normals;
}

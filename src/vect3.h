/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

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

#ifndef OPENMEEG_VECT3_H
#define OPENMEEG_VECT3_H

#if WIN32
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_DEPRECATE 1
#endif
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <assert.h>

#include "DLLDefinesOpenMEEG.h"

namespace OpenMEEG {

    /** \brief  Vect3

        Mesh Class

    **/

    class OPENMEEG_EXPORT Vect3 {

    private:
        double m_x,m_y,m_z; //!< Coordinates of the vector

    public:
        inline Vect3 (const double &x, const double &y, const double &z) : m_x(x),m_y(y),m_z(z) {}
        inline Vect3 (const double &a) : m_x(a),m_y(a),m_z(a) {}
        inline Vect3() {}
        inline ~Vect3() {}
        inline double & x(){return m_x;}
        inline double & y(){return m_y;}
        inline double & z(){return m_z;}
        inline const double& x() const {return m_x;}
        inline const double& y() const {return m_y;}
        inline const double& z() const {return m_z;}

        inline double operator* (const Vect3 &v) const {return m_x*v.x() + m_y*v.y() + m_z*v.z();}
        inline double norm() const {return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);}
        inline double norm2() const {return m_x*m_x+m_y*m_y+m_z*m_z;}
        inline bool operator== (const Vect3 &v ) const {return (m_x==v.x() &&  m_y==v.y() && m_z==v.z());}
        inline void operator+= (const Vect3 &v) {m_x+=v.x(); m_y+=v.y(); m_z+=v.z();}
        inline void operator*= (const double &d) {m_x*=d; m_y*=d; m_z*=d;}
        inline void multadd(const double &d, const Vect3&v) {m_x = m_x+d*v.x(); m_y = m_y+d*v.y(); m_z = m_z+d*v.z();}
        inline Vect3 operator+ (const Vect3&v) const {return Vect3(m_x+v.x(), m_y+v.y(), m_z+v.z());}
        inline Vect3 operator- (const Vect3 &v) const {return Vect3(m_x-v.x(), m_y-v.y(), m_z-v.z());}
        inline Vect3 operator^ (const Vect3 &v) const {return Vect3( m_y*v.z()-m_z*v.y(), m_z*v.x()-m_x*v.z(), m_x*v.y()-m_y*v.x());}
        inline Vect3 operator* (const double &d) const {return Vect3(d*m_x, d*m_y, d*m_z);}
        inline Vect3 operator/ (const double &d) const {double d2 = 1.0/d; return Vect3(d2*m_x, d2*m_y, d2*m_z);}

        inline double operator() (const int i) const
        {
            assert(i>=0 && i<3);
            switch(i)
            {
                case 0 : return m_x;
                case 1 : return m_y;
                case 2 : return m_z;
                default : exit(1);
            }
        }

        inline double& operator() (const int i)
        {
            assert(i>=0 && i<3);
            switch(i)
            {
                case 0 : return m_x;
                case 1 : return m_y;
                case 2 : return m_z;
                default : exit(1);
            }
        }

        inline Vect3 operator- () {return Vect3(-m_x,-m_y,-m_z);}

        inline double det(const Vect3 &y2, const Vect3 &y3) const
        {
            return (*this)*(y2^y3); // y1.det(y2,y3):= y1/(y2^y3)
        }

        inline double solangl(const Vect3 &v1,const Vect3 &v2,const Vect3 &v3,double *omega_i = NULL) const
        {
            double omega;

            // De Munck : Good sign directly
            Vect3 Y1 = v1-*this;
            Vect3 Y2 = v2-*this;
            Vect3 Y3 = v3-*this;
            double y1 = Y1.norm();
            double y2 = Y2.norm();
            double y3 = Y3.norm();
            double d = Y1*(Y2^Y3);
            omega = 2*atan2(d,(y1*y2*y3+y1*(Y2*Y3)+y2*(Y3*Y1)+y3*(Y1*Y2)));

            if (omega_i==NULL)
                return omega;

            Vect3 Z1 = Y2^Y3;
            Vect3 Z2 = Y3^Y1;
            Vect3 Z3 = Y1^Y2;
            Vect3 D1 = Y2-Y1;
            Vect3 D2 = Y3-Y2;
            Vect3 D3 = Y1-Y3;
            double d1 = D1.norm();
            double d2 = D2.norm();
            double d3 = D3.norm();
            double g1 = -1/d1*log((y1*d1+Y1*D1)/(y2*d1+Y2*D1));
            double g2 = -1/d2*log((y2*d2+Y2*D2)/(y3*d2+Y3*D2));
            double g3 = -1/d3*log((y3*d3+Y3*D3)/(y1*d3+Y1*D3));
            Vect3 N = Z1+Z2+Z3;
            double A = N.norm();
            Vect3 S = D1*g1+D2*g2+D3*g3;
            omega_i[0] = 1/A/A*(Z1*N*omega+d*(D2*S));
            omega_i[1] = 1/A/A*(Z2*N*omega+d*(D3*S));
            omega_i[2] = 1/A/A*(Z3*N*omega+d*(D1*S));

            return omega;
        }

        inline Vect3 normal(const Vect3 &v2, const Vect3 &v3) const
        {
            Vect3 v=*this;
            return ( (v-v2)^(v-v3) ) ;
        }

        inline Vect3& normalize() {*this=*this*(1/(*this).norm()); return *this;}

        friend std::ostream& operator<<(std::ostream &os,const Vect3 &v);
        friend std::istream& operator>>(std::istream &is,Vect3 &v);
    };

    inline Vect3 operator * (const double &d, const Vect3 &v) {return v*d;}

    inline std::istream& operator>>(std::istream &is,Vect3 &v)
    {
        return is >> v.x() >> v.y() >> v.z();
    }

    inline std::ostream& operator<<(std::ostream &os,const Vect3 &v)
    {
        return os << v.x() << " " << v.y() << " " << v.z() ;
    }
}

#endif  //! OPENMEEG_VECT3#_H

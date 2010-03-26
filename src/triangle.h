/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
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

namespace OpenMEEG {

    /** \brief  Triangle
        
        Triangle class
        
    **/

    class OPENMEEG_EXPORT Triangle {

    private:
        int m_s1,m_s2,m_s3; //!< index of vertices of the triangle
        double m_area; //!< area of the triangle
        Vect3 n; // Normale

    public:
        inline Triangle(int a, int b, int c, Vect3 m) {
            m_s1=a; m_s2=b; m_s3=c; n=m;
        }

        inline Triangle() {}

        inline ~Triangle() {}

        inline int som(int i) const {
            switch (i){
                case 1:
                    return m_s1;
                case 2:
                    return m_s2;
                case 3:
                    return m_s3;
                default:
                    static int foo;
                    std::cerr << "bad idx in som\n";
                    return foo;
            }
        }

        inline int next(int i) const {
            return som(1+(i%3));
        }

        inline int prev(int i) const {
            return som(1+((1+i)%3));
        }

        inline int& s1() { return m_s1; }
        inline int& s2() { return m_s2; }
        inline int& s3() { return m_s3; }

        inline int s1() const { return m_s1; }
        inline int s2() const { return m_s2; }
        inline int s3() const { return m_s3; }

        inline const Vect3& normal() const { return n; }
        inline Vect3& normal() { return n; }

        inline int contains(int l) const {
            if(m_s1==l)
                return 1;
            if(m_s2==l)
                return 2;
            if(m_s3==l)
                return 3;
            return 0;
        }

        inline double getArea() const { return m_area; };
        inline void setArea( double a ) { m_area = a; };
        inline double& area() { return m_area; }

        inline int operator[] (const int i) const {
            switch(i)
            {
                case 0 : return m_s1;
                case 1 : return m_s2;
                case 2 : return m_s3;
                default : {std::cerr<<"Error in Triangle class: too large index\n"; exit(-1);}
            }
        }

        inline int& operator[] (const int i) {
            switch(i)
            {
                case 0 : return m_s1;
                case 1 : return m_s2;
                case 2 : return m_s3;
                default : {std::cerr<<"Error in Triangle class: too large index\n"; exit(-1);}
            }
        }

        friend std::istream& operator>>(std::istream &is,Triangle &t);
    };

    inline std::istream& operator>>(std::istream &is,Triangle &t)
    {
        return is >> t.m_s1 >> t.m_s2 >> t.m_s3;
    }

    inline std::ostream& operator<<(std::ostream &os,const Triangle &t)
    {
        return os << t[0] << " " << t[1] << " " << t[2];
    }
}

#endif  //! OPENMEEG_TRIANGLE_H

#ifndef H_TRIANGLE
#define H_TRIANGLE

#include <cstdlib>
#include "vect3.h"

/** \brief  Triangle
    
    Triangle class
    
© Copyright 2007 Odyssee INRIA . All Rights Reserved.

    \author $LastChangedBy$
    \date $Date$
    \version $Rev$  \sa
**/

class Triangle {

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

    inline int& som(int i) {
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

    inline int& next(int i) {
        return som(1+(i%3));
    }

    inline int& prev(int i) {
        return som(1+((1+i)%3));
    }

    inline const int& s1() { return m_s1; }
    inline const int& s2() { return m_s2; }
    inline const int& s3() { return m_s3; }

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

    inline void setNormal(const Vect3& v, bool normalize = false) {
        n = v;
        if (normalize) {
            double norm = n.norme();
            assert(norm > 0);
            n = n/norm;
        }
    }
};

inline std::istream& operator>>(std::istream &is,Triangle &t)
{
    return is >> t.m_s1 >> t.m_s2 >> t.m_s3;
}

inline std::ostream& operator<<(std::ostream &os,const Triangle &t)
{
    return os << t[0] << " " << t[1] << " " << t[2] ;
}

#endif



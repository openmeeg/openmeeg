#ifndef H_TRIANGLE
#define H_TRIANGLE

#include <cstdlib>
#include "vect3.h"

class triangle{

private:
    int s1,s2,s3;
    double aire;
    vect3 n;

public:
    inline triangle(int a, int b, int c, vect3 m) {
        s1=a; s2=b; s3=c; n=m;
    }

    inline triangle() {}

    inline ~triangle() {}

    inline int&  som(int i) {
        switch (i){
            case 1:
            return s1;
            case 2:
            return s2;
            case 3:
            return s3;
            default:
            static int foo;
            std::cerr << "bad idx in som\n";
            return foo;
        }
    }

    inline int&  next(int i) {
        return som(1+(i%3));
    }

    inline int&  prev(int i) {
        return som(1+((1+i)%3));
    }

    inline int&  som1() {
        return s1;
    }

    inline int&  som2(){
        return s2;
    }

    inline int&  som3(){
        return s3;
    }

    inline const vect3& normale() const{
        return n;
    }

    inline int appar(int l) const
    {
        if(s1==l)
            return 1;
        if(s2==l)
            return 2;
        if(s3==l)
            return 3;
        return 0;
    }

    inline void retour()
    {
        int tp;
        tp=s3;
        s3=s2;
        s2=tp;
    }

    inline int id1()const { return s1; }
    inline int id2()const { return s2; }
    inline int id3()const { return s3; }

    inline double getAire() const {return aire;};
    inline void setAire( double a ) {aire=a;};
    inline const int operator[] (const int i) const {
        switch(i)
        {
            case 0 : return s1;
            case 1 : return s2;
            case 2 : return s3;
            default : {std::cerr<<"Error in triangle class: too large index\n"; exit(-1);}
        }
    }

    inline int& operator[] (const int i) {
        switch(i)
        {
            case 0 : return s1;
            case 1 : return s2;
            case 2 : return s3;
            default : {std::cerr<<"Error in triangle class: too large index\n"; exit(-1);}
        }
    }

    friend std::istream& operator>>(std::istream &is,triangle &t);

    inline void set_normale(const vect3& v, bool normalize = false) {
        n=v;
        if (normalize) {
            double norm = n.norme();
            assert(norm > 0);
            n = n/norm;
        }
    }
};

inline std::istream& operator>>(std::istream &is,triangle &t)
{
    return is >> t.s1 >> t.s2 >> t.s3;
}

inline std::ostream& operator<<(std::ostream &os,const triangle &t)
{
    return os << t.id1() << " " << t.id2() << " " << t.id3() ;
}

#endif



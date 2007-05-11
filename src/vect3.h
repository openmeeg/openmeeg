#ifndef H_CLASS3VECT
#define H_CLASS3VECT

#if WIN32
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_DEPRECATE 1
#endif
#include <math.h>
#include <iostream>
#include <assert.h>

class vect3 {

private:
    double x,y,z;

public:
    inline vect3 (const double &a, const double &b, const double &c) : x(a),y(b),z(c) {}
    inline vect3 (const double &a) : x(a),y(a),z(a) {}
    inline vect3() {}
    inline ~vect3() {}
    inline double & _x(){return x;}
    inline double & _y(){return y;}
    inline double & _z(){return z;}
    inline const double& X() const {return x;}
    inline const double& Y() const {return y;}
    inline const double& Z() const {return z;}

    inline double operator* (const vect3 &v) const {return x*(v.x) + y*(v.y)+ z*(v.z);}
    inline double norme() const {return sqrt(x*x+y*y+z*z);}
    inline double norme2() const {return x*x+y*y+z*z;}
    inline bool operator== (const vect3 &v ) const {return (x==(v.x) &&  y==(v.y) && z==(v.z));}
    inline void operator+= (const vect3 &v) {x+=v.x; y+=v.y; z+=v.z;}
    inline void operator*= (const double &d) {x*=d; y*=d; z*=d;}
    inline void multadd(const double &d, const vect3&v) {x=x+d*v.x; y=y+d*v.y; z=z+d*v.z;}
    inline vect3 operator+ (const vect3&v) const {return vect3(x+(v.x), y+(v.y), z+(v.z));}
    inline vect3 operator- (const vect3 &v) const {return vect3(x-(v.x), y-(v.y), z-(v.z));}
    inline vect3 operator^ (const vect3 &v) const {return vect3( y*(v.z)-z*(v.y), z*(v.x)-x*(v.z), x*(v.y)-y*(v.x));}
    inline vect3 operator* (const double &d) const {return vect3(d*x, d*y, d*z);}
    inline vect3 operator/ (const double &d) const {double d2=1.0/d; return vect3(d2*x, d2*y, d2*z);}
    
    inline const double operator[] (const int i) const
    {
        assert(i>=0 && i<3);
        switch(i)
        {
            case 0 : return x;
            case 1 : return y;
            case 2 : return z;
        }
    }

    inline double& operator[] (const int i)
    {
        assert(i>=0 && i<3);
        switch(i)
        {
            case 0 : return x;
            case 1 : return y;
            case 2 : return z;
            default: return z;
        }
    }

    inline vect3 operator- () {return vect3(-x,-y,-z);}

    inline double det(const vect3 &y2, const vect3 &y3) const
    {
        return (*this)*(y2^y3); // y1.det(y2,y3):= y1/(y2^y3)
    }

    inline double solangl(const vect3 &v1,const vect3 &v2,const vect3 &v3,double *omega_i=NULL) const
    {
        double omega;

        // De Munck
        // Bon signe directement
        vect3 Y1=v1-*this;
        vect3 Y2=v2-*this;
        vect3 Y3=v3-*this;
        double y1=Y1.norme();
        double y2=Y2.norme();
        double y3=Y3.norme();
        double d=Y1*(Y2^Y3);
        omega=2*atan2(d,(y1*y2*y3+y1*(Y2*Y3)+y2*(Y3*Y1)+y3*(Y1*Y2)));

        if (omega_i==NULL)
            return omega;

        vect3 Z1=Y2^Y3;
        vect3 Z2=Y3^Y1;
        vect3 Z3=Y1^Y2;
        vect3 D1=Y2-Y1;
        vect3 D2=Y3-Y2;
        vect3 D3=Y1-Y3;
        double d1=D1.norme();
        double d2=D2.norme();
        double d3=D3.norme();
        double g1=-1/d1*log((y1*d1+Y1*D1)/(y2*d1+Y2*D1));
        double g2=-1/d2*log((y2*d2+Y2*D2)/(y3*d2+Y3*D2));
        double g3=-1/d3*log((y3*d3+Y3*D3)/(y1*d3+Y1*D3));
        vect3 N=Z1+Z2+Z3;
        double A=N.norme();
        vect3 S=D1*g1+D2*g2+D3*g3;
        omega_i[0]=1/A/A*(Z1*N*omega+d*(D2*S));
        omega_i[1]=1/A/A*(Z2*N*omega+d*(D3*S));
        omega_i[2]=1/A/A*(Z3*N*omega+d*(D1*S));

        return omega;
    }

    inline vect3 normale(const vect3 &v2, const vect3 &v3) const
    {
        vect3 v=*this;
        return ( (v-v2)^(v-v3) ) ;
    }

    inline void normalize() {*this=*this*(1/(*this).norme());}

    friend std::ostream& operator<<(std::ostream &os,const vect3 &v);
    friend std::istream& operator>>(std::istream &is,vect3 &v);

};

inline vect3 operator * (const double &d, const vect3 &v) {return v*d;}

inline std::istream& operator>>(std::istream &is,vect3 &v)
{
    return is >> v.x >> v.y >> v.z;
}

inline std::ostream& operator<<(std::ostream &os,const vect3 &v)
{
    return os << v.x << " " << v.y << " " << v.z ;
}

#endif


#ifndef POINTS_H
#define POINTS_H

#include <iostream>

namespace Geometry {

    template <unsigned N>
    struct Point {

        double& operator()(const unsigned int i)      { return coords[i-1]; }
        double  operator()(const unsigned int i) const{ return coords[i-1]; }

        static const unsigned Dim;

    private:

        double coords[N];
    };

    template <unsigned N> const unsigned Point<N>::Dim = N;

    template <unsigned N>
    std::istream&
    operator>>(std::istream& is,Point<N>& point) {
        for (unsigned i=1;i<=N;++i)
            is >> point(i);
        return is;
    }

    template <unsigned N>
    std::ostream&
    operator<<(std::ostream& os,const Point<N>& point) {
        for (unsigned i=1;i<=N;++i)
            os << point(i) << ' ';
        return os;
    }
}

#endif  //  ! POINTS_H

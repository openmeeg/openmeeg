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

#include <cmath>
#include <iostream>

#include <vertex.h>
#include <triangle.h>
#include <mesh.h>

namespace OpenMEEG {

    // light class containing d Vect3
    template <int d>
    class OPENMEEG_EXPORT Vect3array 
    {
        Vect3 t[d];

    public:
        Vect3array() {};
        inline Vect3array(double x) {
            for ( unsigned i = 0; i < d; ++i)
                t[i] = Vect3(x);
        }
        inline Vect3array<d> operator*(double x) const {
            Vect3array<d> r;
            for ( unsigned i = 0; i < d; ++i)
                r.t[i] = t[i]*x;
            return r;
        }
        inline Vect3 operator()(int i) const { return t[i]; }
        inline Vect3& operator()(int i) { return t[i]; }
    };

    template <int d>
    inline void multadd (Vect3array<d>& target, const double scale,  const Vect3array<d>& incr) 
    {
        for ( unsigned i = 0; i < d; ++i) {
            target(i) = target(i) + scale*incr(i);
        }
    }

    inline void multadd (double& target, const double scale, const double incr) 
    {
        target += scale*incr;
    }

    inline void multadd (Vect3& target, const double scale,  const Vect3& incr) 
    {
        target = target + scale*incr;
    }

    // Quadrature rules are from Marc Bonnet's book: Equations integrales..., Appendix B.3

    static const double cordBars[4][16][4] =
    {
        //parameters for N=3
        {
            {0.166666666666667, 0.166666666666667, 0.666666666666667, 0.166666666666667},
            {0.166666666666667, 0.666666666666667, 0.166666666666667, 0.166666666666667},
            {0.666666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        }
        ,
        // parameters for N=6
        {
            {0.445948490915965, 0.445948490915965, 0.108103018168070, 0.111690794839005},
            {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.111690794839005},
            {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.111690794839005},
            {0.091576213509771, 0.091576213509771, 0.816847572980458, 0.054975871827661},
            {0.091576213509771, 0.816847572980458, 0.091576213509771, 0.054975871827661},
            {0.816847572980458, 0.091576213509771, 0.091576213509771, 0.054975871827661},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        }
        ,
            // parameters for N=7
        {
            {0.333333333333333, 0.333333333333333, 0.333333333333333, 0.1125},
            {0.470142064105115, 0.470142064105115, 0.059715871789770, 0.066197076394253},
            {0.470142064105115, 0.059715871789770, 0.470142064105115, 0.066197076394253},
            {0.059715871789770, 0.470142064105115, 0.470142064105115, 0.066197076394253},
            {0.101286507323456, 0.101286507323456, 0.797426985353088, 0.062969590272414},
            {0.101286507323456, 0.797426985353088, 0.101286507323456, 0.062969590272414},
            {0.797426985353088, 0.101286507323456, 0.101286507323456, 0.062969590272414},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        }
        ,

            // parameters for N=16
        {
            {0.333333333333333, 0.333333333333333, 0.333333333333333, 0.072157803838893},
            {0.081414823414554, 0.459292588292722, 0.459292588292722, 0.047545817133642},
            {0.459292588292722, 0.081414823414554, 0.459292588292722, 0.047545817133642},
            {0.459292588292722, 0.459292588292722, 0.081414823414554, 0.047545817133642},
            {0.898905543365937, 0.050547228317031, 0.050547228317031, 0.016229248811599},
            {0.050547228317031, 0.898905543365937, 0.050547228317031, 0.016229248811599},
            {0.050547228317031, 0.050547228317031, 0.898905543365937, 0.016229248811599},
            {0.658861384496479, 0.170569307751760, 0.170569307751761, 0.051608685267359},
            {0.170569307751760, 0.658861384496479, 0.170569307751761, 0.051608685267359},
            {0.170569307751760, 0.170569307751761, 0.658861384496479, 0.051608685267359},
            {0.008394777409957, 0.728492392955404, 0.263112829634639, 0.013615157087217},
            {0.728492392955404, 0.008394777409957, 0.263112829634639, 0.013615157087217},
            {0.728492392955404, 0.263112829634639, 0.008394777409957, 0.013615157087217},
            {0.008394777409957, 0.263112829634639, 0.728492392955404, 0.013615157087217},
            {0.263112829634639, 0.008394777409957, 0.728492392955404, 0.013615157087217},
            {0.263112829634639, 0.728492392955404, 0.008394777409957, 0.013615157087217}
        }

    }; // end of gaussTriangleParams

    static const unsigned nbPts[4] = {3, 6, 7, 16};

    template <typename T, typename I>
    class OPENMEEG_EXPORT Integrator 
    {
        unsigned order;

    public:

        inline Integrator()             { setOrder(3);   }
        inline Integrator(unsigned ord) { setOrder(ord); }
        inline virtual ~Integrator() {}

        inline void setOrder(const unsigned n) 
        {
            if ( n < 4 ) {
                order = n;
            } else {
                std::cout << "Unavailable Gauss order: min is 1, max is 3" << n << std::endl;
                order = (n < 1) ? 1 : 3;
            }
        }

        virtual inline T integrate(const I& fc, const Triangle& Trg)
        {
            const Vect3 points[3] = { Trg.s1(), Trg.s2(), Trg.s3() };
            return triangle_integration(fc, points);
        }

    protected:

        inline T triangle_integration(const I& fc, const Vect3 points[3])
        {
            // compute double area of triangle defined by points
            Vect3 crossprod = (points[1] - points[0])^(points[2] - points[0]);
            double S = crossprod.norm();
            T result = 0;
            for ( unsigned i = 0; i < nbPts[order];++i) {
                Vect3 v(0.0, 0.0, 0.0);
                for ( unsigned j = 0; j < 3; ++j) {
                    v.multadd(cordBars[order][i][j], points[j]);
                }
                multadd(result, cordBars[order][i][3], fc.f(v));
            }
            return result*S;
        }
    };

    template <typename T, typename I>
    class OPENMEEG_EXPORT AdaptiveIntegrator: public Integrator<T, I>
    {
        typedef Integrator<T, I> base;

    public:

        inline AdaptiveIntegrator() : tolerance(0.0001) {}
        inline AdaptiveIntegrator(double tol) : tolerance(tol) {}
        inline ~AdaptiveIntegrator() {}

        inline double norm(const double a) { return fabs(a);  }
        inline double norm(const Vect3& a) { return a.norm(); }

        virtual inline T integrate(const I& fc, const Triangle& Trg) {
            const Vect3 points[3] = { Trg.s1(), Trg.s2(), Trg.s3() };
            T I0 = base::triangle_integration(fc, points);
            return adaptive_integration(fc, points, I0, 0);
        }

    private:

        double tolerance;

        inline T adaptive_integration(const I& fc, const Vect3 * points, T I0, unsigned n)
        {
            Vect3 newpoint0(0.0, 0.0, 0.0);
            multadd(newpoint0, 0.5, points[0]);
            multadd(newpoint0, 0.5, points[1]);
            Vect3 newpoint1(0.0, 0.0, 0.0);
            multadd(newpoint1, 0.5, points[1]);
            multadd(newpoint1, 0.5, points[2]);
            Vect3 newpoint2(0.0, 0.0, 0.0);
            multadd(newpoint2, 0.5, points[2]);
            multadd(newpoint2, 0.5, points[0]);
            Vect3 points1[3] = {points[0], newpoint0, newpoint2};
            Vect3 points2[3] = {points[1], newpoint1, newpoint0};
            Vect3 points3[3] = {points[2], newpoint2, newpoint1};
            Vect3 points4[3] = {newpoint0, newpoint1, newpoint2};
            T I1 = base::triangle_integration(fc, points1);
            T I2 = base::triangle_integration(fc, points2);
            T I3 = base::triangle_integration(fc, points3);
            T I4 = base::triangle_integration(fc, points4);
            T sum = I1+I2+I3+I4;
            if ( norm(I0-sum) > tolerance*norm(I0) ) {
                n = n+1;
                if ( n < 10 ) {
                    I1 = adaptive_integration(fc, points1, I1, n);
                    I2 = adaptive_integration(fc, points2, I2, n);
                    I3 = adaptive_integration(fc, points3, I3, n);
                    I4 = adaptive_integration(fc, points4, I4, n);
                    I0 = I1+I2+I3+I4;
                }
            }
            return I0;
        }
    };
}

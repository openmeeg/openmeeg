// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <cmath>
#include <iostream>

#include <vertex.h>
#include <triangle.h>
#include <mesh.h>

namespace OpenMEEG {

    class OPENMEEG_EXPORT Integrator {

        typedef Vect3 Point;
        typedef Point TrianglePoints[3];

        static unsigned safe_order(const unsigned n) {
            if (n>0 && n<4)
                return n;
            std::cout << "Unavailable Gauss order " << n << ": min is 1, max is 3" << std::endl;
            return (n<1) ? 1 : 3;
        }

    public:

        Integrator(const unsigned ord): Integrator(ord,0,0.0) { }
        Integrator(const unsigned ord,const double tol): Integrator(ord,10,tol) { }
        Integrator(const unsigned ord,const unsigned levels,const double tol=0.0001):
            order(safe_order(ord)),tolerance(tol),max_depth(levels)
        { }

        double norm(const double a) const { return fabs(a);  }
        double norm(const Vect3& a) const { return a.norm(); }

        // TODO: T can be deduced from Function.

        #ifndef SWIGPYTHON  // SWIG sees the integrate def as a syntax error
        template <typename Function>
        decltype(auto) integrate(const Function& function,const Triangle& triangle) const {
            const TrianglePoints tripts = { triangle.vertex(0), triangle.vertex(1), triangle.vertex(2) };
            const auto& coarse = triangle_integration(function,tripts);
            return (max_depth==0) ? coarse : adaptive_integration(function,tripts,coarse,max_depth);
        }
        #endif /* not SWIGPYTHON */

    private:

        #ifndef SWIGPYTHON  // SWIG sees the integrate def as a syntax error
        template <typename Function>
        decltype(auto) triangle_integration(const Function& function,const TrianglePoints& triangle) const {
            using T = decltype(function(Vect3()));
            T result = 0.0;
            for (unsigned i=0;i<nbPts[order];++i) {
                Vect3 v(0.0,0.0,0.0);
                for (unsigned j=0; j<3; ++j)
                    v.multadd(rules[order][i].barycentric_coordinates[j],triangle[j]);
                result += rules[order][i].weight*function(v);
            }

            // compute double area of triangle defined by points

            const double area2 = crossprod(triangle[1]-triangle[0],triangle[2]-triangle[0]).norm();
            return result*area2;
        }
        #endif /* not SWIGPYTHON */

        template <typename T,typename Function>
        T adaptive_integration(const Function& function,const TrianglePoints& triangle,const T& coarse,const unsigned level) const {
            const Point midpoints[] = { 0.5*(triangle[1]+triangle[2]), 0.5*(triangle[2]+triangle[0]), 0.5*(triangle[0]+triangle[1]) };
            const TrianglePoints new_triangles[] = {
                { triangle[0],  midpoints[1], midpoints[2] }, { midpoints[0], triangle[1],  midpoints[2] },
                { midpoints[0], midpoints[1], triangle[2]  }, { midpoints[0], midpoints[1], midpoints[2] }
            };

            T refined = 0.0;
            T integrals[4];
            for (unsigned i=0; i<4; ++i) {
                integrals[i] = triangle_integration(function,new_triangles[i]);
                refined += integrals[i];
            }

            if (norm(coarse-refined)<=tolerance*norm(coarse) || level==0)
                return refined;

            T sum = 0.0;
            for (unsigned i=0; i<4; ++i)
                sum += adaptive_integration(function,new_triangles[i],integrals[i],level-1);
            return sum;
        }

        static constexpr unsigned nbPts[4] = { 3, 6, 7, 16 };

        const unsigned order;
        const double   tolerance;
        const unsigned max_depth;

        // Quadrature rules are from Marc Bonnet's book: Equations integrales..., Appendix B.3

        struct QuadratureRule {
            double barycentric_coordinates[3];
            double weight;
        };

        static constexpr QuadratureRule rules[4][16] = {
            {   // Parameters for N=3
                {{ 0.166666666666667, 0.166666666666667, 0.666666666666667 }, 0.166666666666667 },
                {{ 0.166666666666667, 0.666666666666667, 0.166666666666667 }, 0.166666666666667 },
                {{ 0.666666666666667, 0.166666666666667, 0.166666666666667 }, 0.166666666666667 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }
            },
            {   // Parameters for N=6
                {{ 0.445948490915965, 0.445948490915965, 0.108103018168070 }, 0.111690794839005 },
                {{ 0.445948490915965, 0.108103018168070, 0.445948490915965 }, 0.111690794839005 },
                {{ 0.108103018168070, 0.445948490915965, 0.445948490915965 }, 0.111690794839005 },
                {{ 0.091576213509771, 0.091576213509771, 0.816847572980458 }, 0.054975871827661 },
                {{ 0.091576213509771, 0.816847572980458, 0.091576213509771 }, 0.054975871827661 },
                {{ 0.816847572980458, 0.091576213509771, 0.091576213509771 }, 0.054975871827661 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }
            },
            {   // Parameters for N=7
                {{ 0.333333333333333, 0.333333333333333, 0.333333333333333 }, 0.1125 },
                {{ 0.470142064105115, 0.470142064105115, 0.059715871789770 }, 0.066197076394253 },
                {{ 0.470142064105115, 0.059715871789770, 0.470142064105115 }, 0.066197076394253 },
                {{ 0.059715871789770, 0.470142064105115, 0.470142064105115 }, 0.066197076394253 },
                {{ 0.101286507323456, 0.101286507323456, 0.797426985353088 }, 0.062969590272414 },
                {{ 0.101286507323456, 0.797426985353088, 0.101286507323456 }, 0.062969590272414 },
                {{ 0.797426985353088, 0.101286507323456, 0.101286507323456 }, 0.062969590272414 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }
            },
            {   // Parameters for N=16
                {{ 0.333333333333333, 0.333333333333333, 0.333333333333333 }, 0.072157803838893 },
                {{ 0.081414823414554, 0.459292588292722, 0.459292588292722 }, 0.047545817133642 },
                {{ 0.459292588292722, 0.081414823414554, 0.459292588292722 }, 0.047545817133642 },
                {{ 0.459292588292722, 0.459292588292722, 0.081414823414554 }, 0.047545817133642 },
                {{ 0.898905543365937, 0.050547228317031, 0.050547228317031 }, 0.016229248811599 },
                {{ 0.050547228317031, 0.898905543365937, 0.050547228317031 }, 0.016229248811599 },
                {{ 0.050547228317031, 0.050547228317031, 0.898905543365937 }, 0.016229248811599 },
                {{ 0.658861384496479, 0.170569307751760, 0.170569307751761 }, 0.051608685267359 },
                {{ 0.170569307751760, 0.658861384496479, 0.170569307751761 }, 0.051608685267359 },
                {{ 0.170569307751760, 0.170569307751761, 0.658861384496479 }, 0.051608685267359 },
                {{ 0.008394777409957, 0.728492392955404, 0.263112829634639 }, 0.013615157087217 },
                {{ 0.728492392955404, 0.008394777409957, 0.263112829634639 }, 0.013615157087217 },
                {{ 0.728492392955404, 0.263112829634639, 0.008394777409957 }, 0.013615157087217 },
                {{ 0.008394777409957, 0.263112829634639, 0.728492392955404 }, 0.013615157087217 },
                {{ 0.263112829634639, 0.008394777409957, 0.728492392955404 }, 0.013615157087217 },
                {{ 0.263112829634639, 0.728492392955404, 0.008394777409957 }, 0.013615157087217 }
            }
        };
    };
}

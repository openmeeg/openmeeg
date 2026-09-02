// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <triangle.h>

namespace OpenMEEG {

    class OPENMEEG_EXPORT Integrator {

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

        #ifndef SWIG  // SWIG sees the integrate def as a syntax error
        template <typename Function>
        auto integrate(const Function& function,const Triangle& triangle) const {
            using ResType = decltype(function(Point3D()));
            const TrianglePoints tripts = { triangle.vertex(0), triangle.vertex(1), triangle.vertex(2) };
            const ResType coarse = triangle_integration<ResType>(function,tripts,triangle.area());
            return (max_depth==0) ? coarse : adaptive_integration<ResType>(function,tripts,triangle.area(),coarse,max_depth);
        }
        #endif

    private:

        template <typename ResType,typename Function>
        ResType triangle_integration(const Function& function,const TrianglePoints& triangle,const double area) const {
            ResType result = 0.0;
            for (unsigned i=0;i<nbPts[order];++i)
                result += rules[order][i].weight*function(cartesian_coordinates(TriangleCoords(rules[order][i].barycentric_coordinates),triangle));
            return area*result;
        }

        template <typename ResType,typename Function>
        ResType adaptive_integration(const Function& function,const TrianglePoints& triangle,const double area,const ResType& coarse,const unsigned level) const {
            const Point3D midpoints[] = { 0.5*(triangle[1]+triangle[2]), 0.5*(triangle[2]+triangle[0]), 0.5*(triangle[0]+triangle[1]) };
            const TrianglePoints new_triangles[] = {
                { triangle[0],  midpoints[1], midpoints[2] }, { midpoints[0], triangle[1],  midpoints[2] },
                { midpoints[0], midpoints[1], triangle[2]  }, { midpoints[0], midpoints[1], midpoints[2] }
            };

            const double subdivided_area = 0.25*area;

            ResType refined = 0.0;
            ResType integrals[4];
            for (unsigned i=0; i<4; ++i) {
                integrals[i] = triangle_integration<ResType>(function,new_triangles[i],subdivided_area);
                refined += integrals[i];
            }

            if (norm(coarse-refined)<=tolerance*norm(coarse) || level==0)
                return refined;

            ResType sum = 0.0;
            for (unsigned i=0; i<4; ++i)
                sum += adaptive_integration<ResType>(function,new_triangles[i],subdivided_area,integrals[i],level-1);
            return sum;
        }

        static constexpr unsigned nbPts[4] = { 3, 6, 7, 16 };

        const unsigned order;
        const double   tolerance;
        const unsigned max_depth;

        // Quadrature rules are from Marc Bonnet's book: Equations integrales..., Appendix B.3

        struct QuadraturePoint {
            double barycentric_coordinates[3];
            double weight;
        };

        using QuadratureRule = QuadraturePoint[16];

        static constexpr QuadratureRule rules[4] = {
            {   // Parameters for N=3
                {{ 0.166666666666667, 0.166666666666667, 0.666666666666667 }, 0.333333333333334 },
                {{ 0.166666666666667, 0.666666666666667, 0.166666666666667 }, 0.333333333333334 },
                {{ 0.666666666666667, 0.166666666666667, 0.166666666666667 }, 0.333333333333334 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }
            },
            {   // Parameters for N=6
                {{ 0.445948490915965, 0.445948490915965, 0.108103018168070 }, 0.223381589678010 },
                {{ 0.445948490915965, 0.108103018168070, 0.445948490915965 }, 0.223381589678010 },
                {{ 0.108103018168070, 0.445948490915965, 0.445948490915965 }, 0.223381589678010 },
                {{ 0.091576213509771, 0.091576213509771, 0.816847572980458 }, 0.109951743655322 },
                {{ 0.091576213509771, 0.816847572980458, 0.091576213509771 }, 0.109951743655322 },
                {{ 0.816847572980458, 0.091576213509771, 0.091576213509771 }, 0.109951743655322 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }
            },
            {   // Parameters for N=7
                {{ 0.333333333333333, 0.333333333333333, 0.333333333333333 }, 0.2250 },
                {{ 0.470142064105115, 0.470142064105115, 0.059715871789770 }, 0.132394152788506 },
                {{ 0.470142064105115, 0.059715871789770, 0.470142064105115 }, 0.132394152788506 },
                {{ 0.059715871789770, 0.470142064105115, 0.470142064105115 }, 0.132394152788506 },
                {{ 0.101286507323456, 0.101286507323456, 0.797426985353088 }, 0.125939180544828 },
                {{ 0.101286507323456, 0.797426985353088, 0.101286507323456 }, 0.125939180544828 },
                {{ 0.797426985353088, 0.101286507323456, 0.101286507323456 }, 0.125939180544828 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 }, {{ 0.0, 0.0, 0.0 }, 0.0 },
                {{ 0.0, 0.0, 0.0 }, 0.0 }
            },
            {   // Parameters for N=16
                {{ 0.333333333333333, 0.333333333333333, 0.333333333333333 }, 0.144315607677786 },
                {{ 0.081414823414554, 0.459292588292722, 0.459292588292722 }, 0.095091634267284 },
                {{ 0.459292588292722, 0.081414823414554, 0.459292588292722 }, 0.095091634267284 },
                {{ 0.459292588292722, 0.459292588292722, 0.081414823414554 }, 0.095091634267284 },
                {{ 0.898905543365937, 0.050547228317031, 0.050547228317031 }, 0.032458497623198 },
                {{ 0.050547228317031, 0.898905543365937, 0.050547228317031 }, 0.032458497623198 },
                {{ 0.050547228317031, 0.050547228317031, 0.898905543365937 }, 0.032458497623198 },
                {{ 0.658861384496479, 0.170569307751760, 0.170569307751761 }, 0.103217370534718 },
                {{ 0.170569307751760, 0.658861384496479, 0.170569307751761 }, 0.103217370534718 },
                {{ 0.170569307751760, 0.170569307751761, 0.658861384496479 }, 0.103217370534718 },
                {{ 0.008394777409957, 0.728492392955404, 0.263112829634639 }, 0.027230314174434 },
                {{ 0.728492392955404, 0.008394777409957, 0.263112829634639 }, 0.027230314174434 },
                {{ 0.728492392955404, 0.263112829634639, 0.008394777409957 }, 0.027230314174434 },
                {{ 0.008394777409957, 0.263112829634639, 0.728492392955404 }, 0.027230314174434 },
                {{ 0.263112829634639, 0.008394777409957, 0.728492392955404 }, 0.027230314174434 },
                {{ 0.263112829634639, 0.728492392955404, 0.008394777409957 }, 0.027230314174434 }
            }
        };
    };
}

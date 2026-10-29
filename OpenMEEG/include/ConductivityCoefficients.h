// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <mesh.h>
#include <domain.h>
#include <constants.h>

namespace OpenMEEG {

    class Conductivities {

        Conductivities(const DomainsRefs& domains,const int orientation) {
            const double factor = orientation*K;
            Ncond =  factor*eval_on_common_domains<IDENTITY>(domains); // Sum of the conductivity(ies) of the shared domain(s).
            Scond =  factor*eval_on_common_domains<INVERSE>(domains);  // Sum of the inverse of conductivity(ies) of the shared domain(s).
            Dcond = -factor*domains.size();                            // Opposite of the sum of the indicator function of the shared domain(s).
        }

    public:

        Conductivities(const DomainsRefs& doms): Conductivities(doms,1)                                        { } // Single mesh constructor
        Conductivities(const MeshPair& mp):      Conductivities(mp.common_domains(),mp.relative_orientation()) { } // Adjacent meshes constructor.

        double N() const { return Ncond; }
        double S() const { return Scond; }
        double D() const { return Dcond; }

        static double D(const int orientation) { return orientation*K; }

    private:

        //  Accumulate a function over the domain common to two meshes.

        static double IDENTITY(const Domain& domain) { return domain.conductivity();        }
        static double INVERSE(const Domain& domain)  { return 1.0/domain.conductivity();    }

        template <double Function(const Domain&)>
        double eval_on_common_domains(const DomainsRefs& domains) const {
            double result = 0.0;
            for (const Domain* domainptr : domains)
                result += Function(*domainptr);
            return result;
        }

        static constexpr double K = 1.0/(4*Pi);

        double Ncond;
        double Scond;
        double Dcond;
    };
}

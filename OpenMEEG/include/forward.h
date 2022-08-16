// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include "matrix.h"

namespace OpenMEEG {

    class Forward : public virtual Matrix {
    public:

        Forward(const Matrix& GainMatrix,const Matrix& RealSourcesData,const double NoiseLevel) {

            Matrix& SimulatedData = *this;
            SimulatedData = GainMatrix*RealSourcesData;

            if (NoiseLevel>0) {

                std::random_device rd{};
                std::mt19937 gen{rd()};
                std::normal_distribution<> noise{0,NoiseLevel};

                for (unsigned i=0; i<SimulatedData.nlin(); ++i)
                    for (unsigned j=0; j<SimulatedData.ncol(); ++j)
                        SimulatedData(i,j) += noise(gen);
            }
        }

        virtual ~Forward() { }
    };
}

/*
Project Name: OpenMEEG

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

#include <vector>

#include <vector.h>
#include <matrix.h>
#include <symmatrix.h>
#include <geometry.h>
#include <sensors.h>

namespace OpenMEEG {

    class OPENMEEG_EXPORT HeadMat: public virtual SymMatrix {
    public:
        HeadMat (const Geometry& geo, const unsigned gauss_order=3);
        virtual ~HeadMat () {};
    };

    class OPENMEEG_EXPORT SurfSourceMat: public virtual Matrix {
    public:
        SurfSourceMat (const Geometry& geo, Mesh& sources, const unsigned gauss_order=3);
        virtual ~SurfSourceMat () {};
    };

    class OPENMEEG_EXPORT DipSourceMat: public virtual Matrix {
    public:
        DipSourceMat (const Geometry& geo, const Matrix& dipoles, const unsigned gauss_order=3,
                      const bool adapt_rhs = true, const std::string& domain_name = "");
        virtual ~DipSourceMat () {};
    };

    class OPENMEEG_EXPORT EITSourceMat: public virtual Matrix {
    public:
        EITSourceMat(const Geometry& geo, const Sensors& electrodes, const unsigned gauss_order=3);
        virtual ~EITSourceMat () {};
    };

    class OPENMEEG_EXPORT Surf2VolMat: public virtual Matrix {
    public:
        using Matrix::operator=;
        Surf2VolMat(const Geometry& geo, const Matrix& points);
        virtual ~Surf2VolMat () {};
    };

    class OPENMEEG_EXPORT Head2EEGMat: public virtual SparseMatrix {
    public:
        Head2EEGMat (const Geometry& geo, const Sensors& electrodes);
        virtual ~Head2EEGMat () {};
    };

    class OPENMEEG_EXPORT Head2ECoGMat: public virtual SparseMatrix {
    public:
        Head2ECoGMat (const Geometry& geo, const Sensors& electrodes, const Interface& i);
        Head2ECoGMat (const Geometry& geo, const Sensors& electrodes, const std::string& id); // mainly for SWIG
        virtual ~Head2ECoGMat () {};
    };

    class OPENMEEG_EXPORT Head2MEGMat: public virtual Matrix {
    public:
        Head2MEGMat (const Geometry& geo, const Sensors& sensors);
        virtual ~Head2MEGMat () {};
    };

    class OPENMEEG_EXPORT SurfSource2MEGMat: public virtual Matrix {
    public:
        SurfSource2MEGMat (const Mesh& sources, const Sensors& sensors);
        virtual ~SurfSource2MEGMat () {};
    };

    class OPENMEEG_EXPORT DipSource2MEGMat: public virtual Matrix {
    public:
        DipSource2MEGMat(const Matrix& dipoles, const Sensors& sensors);
        virtual ~DipSource2MEGMat () {};
    };

    class OPENMEEG_EXPORT DipSource2InternalPotMat: public virtual Matrix {
    public:
        DipSource2InternalPotMat(const Geometry& geo, const Matrix& dipoles,
                                 const Matrix& points, const std::string& domain_name = "");
        virtual ~DipSource2InternalPotMat () {};
    };

    class OPENMEEG_EXPORT CorticalMat: public virtual Matrix {
    public:
        CorticalMat (const Geometry& geo, const Head2EEGMat& M, const std::string& domain_name = "CORTEX",
                const unsigned gauss_order=3, double alpha=-1., double beta=-1., const std::string &filename="");
        virtual ~CorticalMat () {};
    };

    class OPENMEEG_EXPORT CorticalMat2: public virtual Matrix {
    public:
        CorticalMat2(const Geometry& geo, const Head2EEGMat& M, const std::string& domain_name = "CORTEX",
                const unsigned gauss_order=3, double gamma=1., const std::string &filename="");
        virtual ~CorticalMat2() {};
    };
}

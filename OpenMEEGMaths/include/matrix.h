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

// #include <OpenMEEGMathsConfig.h>

namespace OpenMEEG {

    // namespace maths {
    //     struct OPENMEEGMATHS_EXPORT MathsIO;
    // }

    // to properly convert a int int to an int
    // OPENMEEGMATHS_EXPORT inline BLAS_INT sizet_to_int(const int& num)
    // {
    //     BLAS_INT num_out = static_cast<BLAS_INT>(num);
    //     // om_assert(num_out >= 0);
    //     return num_out;
    // }

    class LinOpInfo {
    public:

        // typedef maths::MathsIO* IO;


        LinOpInfo() { }
        LinOpInfo(const int m,const int n):
            num_lines(m),num_cols(n)  { }

        virtual ~LinOpInfo() {};

        // LinOpInfo& operator=(const LinOpInfo& l) {
        //     num_lines = l.num_lines;
        //     num_cols  = l.num_cols;
        //     storage   = l.storage;
        //     dim       = l.dim;
        //     return *this;
        // }

        // int  nlin() const { return num_lines; }
        // int& nlin()       { return num_lines; }

        // virtual int  ncol() const { return num_cols; }
        //         int& ncol()       { return num_cols; }

        // StorageType  storageType() const { return storage; }
        // StorageType& storageType()       { return storage; }

        // Dimension   dimension()   const { return dim;     }
        // Dimension&  dimension()         { return dim;     }

        // IO& default_io() { return DefaultIO; }

    protected:

        int            num_lines;
        int            num_cols;
        // IO                DefaultIO;
    };

    class LinOp: public LinOpInfo {

        typedef LinOpInfo base;

    public:

        LinOp() { }
        LinOp(const int m,const int n): base(m,n) { }

    };

    // typedef enum { DEEP_COPY } DeepCopy;

}

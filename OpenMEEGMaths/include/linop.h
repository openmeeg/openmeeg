// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <cstdlib>
#include <cmath>
#include <memory>

#include <OpenMEEGConfigure.h>
#include "OpenMEEGMathsConfig.h"
#include <OMassert.H>

namespace OpenMEEG {

    namespace maths {
        struct OPENMEEGMATHS_EXPORT MathsIO;
    }

    // Properly convert an unsigned int to a BLAS_INT

    inline BLAS_INT sizet_to_int(const unsigned& num) {
        const BLAS_INT num_out = static_cast<BLAS_INT>(num);
        om_assert(num_out>=0);
        return num_out;
    }

    typedef unsigned Dimension;
    typedef unsigned Index;

    class OPENMEEGMATHS_EXPORT LinOpInfo {
    public:

        typedef maths::MathsIO* IO;

        typedef enum { FULL, SYMMETRIC, BLOCK, BLOCK_SYMMETRIC, SPARSE } StorageType;

        LinOpInfo() { }
        LinOpInfo(const Dimension m,const Dimension n,const StorageType st,const unsigned d):
            num_lines(m),num_cols(n),storage(st),dim(d)  { }

        virtual ~LinOpInfo() {};

        Dimension  nlin() const { return num_lines; }
        Dimension& nlin()       { return num_lines; }

        virtual Dimension  ncol() const { return num_cols; }
                Dimension& ncol()       { return num_cols; }

        StorageType  storageType() const { return storage; }
        StorageType& storageType()       { return storage; }

        unsigned   dimension()   const { return dim;     }
        unsigned&  dimension()         { return dim;     }

        IO& default_io() { return DefaultIO; }

    protected:

        Dimension   num_lines;
        Dimension   num_cols;
        StorageType storage;
        unsigned    dim;
        IO          DefaultIO = nullptr;
    };

    class OPENMEEGMATHS_EXPORT LinOp: public LinOpInfo {

        typedef LinOpInfo base;

    public:

        LinOp() { }
        LinOp(const Dimension m,const Dimension n,const StorageType st,const unsigned d): base(m,n,st,d) { }

        virtual size_t size() const = 0;
        virtual void   info() const = 0;
    };

    typedef enum { DEEP_COPY } DeepCopy;

    struct OPENMEEGMATHS_EXPORT LinOpValue: public std::shared_ptr<double[]> {
        typedef std::shared_ptr<double[]> base;

        LinOpValue(): base(0) { }
        LinOpValue(const size_t n): base(new double[n]) { }
        LinOpValue(const size_t n,const double* initval): LinOpValue(n) { std::copy(initval,initval+n,&(*this)[0]); }
        LinOpValue(const size_t n,const LinOpValue& v): LinOpValue(n,&(v[0])) { }

        ~LinOpValue() { }

        bool empty() const { return static_cast<bool>(*this); }
    };
}

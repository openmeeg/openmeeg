/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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

#ifndef OPENMEEG_TRIANGULAR_MATRIX_H
#define OPENMEEG_TRIANGULAR_MATRIX_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>

#include "om_utils.h"
#include "linop.h"
#include "vector.h"
#include "matrix.h"
#include "diagmatrix.h"

namespace OpenMEEG {

    class LowerTriangularMatrix;

    class OPENMEEGMATHS_EXPORT TriangularMatrix : public LinOp {

        friend class Vector;

        utils::RCPtr<LinOpValue> value;

    public:

        TriangularMatrix(): LinOp(0,0,SYMMETRIC,2),value() {}

        // TriangularMatrix(const char* fname): LinOp(0,0,SYMMETRIC,2),value() { this->load(fname); }
        TriangularMatrix(size_t N): LinOp(N,N,SYMMETRIC,2),value(new LinOpValue(size())) { }
        TriangularMatrix(size_t M,size_t N): LinOp(N,N,SYMMETRIC,2),value(new LinOpValue(size())) { assert(N==M); }
        TriangularMatrix(const TriangularMatrix& S,const DeepCopy): LinOp(S.nlin(),S.nlin(),SYMMETRIC,2),value(new LinOpValue(S.size(),S.data())) { }

        ~TriangularMatrix() {};

        size_t size() const { return nlin()*(nlin()+1)/2; };

        size_t  ncol() const { return nlin(); } 
        size_t& ncol()       { return nlin(); }

        bool empty() const { return value->empty(); }

        double* data() const { return value->data; }
    };

    class OPENMEEGMATHS_EXPORT UpperTriangularMatrix : public TriangularMatrix {

    public:
        UpperTriangularMatrix(): TriangularMatrix() {}

        // UpperTriangularMatrix(const char* fname): TriangularMatrix(fname) {}
        UpperTriangularMatrix(size_t N)         : TriangularMatrix(N) {}
        UpperTriangularMatrix(size_t M,size_t N): TriangularMatrix(M,N) {}
        UpperTriangularMatrix(const UpperTriangularMatrix& S,const DeepCopy): TriangularMatrix(S,DEEP_COPY) { }

        ~UpperTriangularMatrix() {};

        inline double operator()( size_t i, size_t j ) const {
            assert(i < nlin() && j < ncol());
            if(i<=j)
                return data()[i+j*(j+1)/2];
            else
                return 0.;
        }

        inline double& operator()( size_t i, size_t j ) {
            assert(i < nlin() && j < ncol() && i <= j);
            return data()[i+j*(j+1)/2];
        }

        inline UpperTriangularMatrix operator*(const DiagMatrix& M) const;

        inline SymMatrix operator*(const LowerTriangularMatrix& L) const;

        inline Vector operator*(const Vector& v) const;

        inline UpperTriangularMatrix inverse() const;

        inline LowerTriangularMatrix transpose() const;
        
        inline void info() const;
    };

    class OPENMEEGMATHS_EXPORT LowerTriangularMatrix : public TriangularMatrix {

    public:
        LowerTriangularMatrix(): TriangularMatrix() {}

        // LowerTriangularMatrix(const char* fname): TriangularMatrix(fname) {}
        LowerTriangularMatrix(size_t N): TriangularMatrix(N) {}
        LowerTriangularMatrix(size_t M,size_t N): TriangularMatrix(M,N) {}
        LowerTriangularMatrix(const LowerTriangularMatrix& S,const DeepCopy): TriangularMatrix(S,DEEP_COPY) { }

        ~LowerTriangularMatrix() {};

        inline double operator()( size_t i, size_t j ) const {
            assert(i < nlin() && j < ncol());
            if(i>=j)
                return data()[j+i*(i+1)/2];
            else
                return 0.;
        } 

        inline double& operator()( size_t i, size_t j ) {
            assert(i < nlin() && j < ncol() && i >= j);
            return data()[j+i*(i+1)/2];
        }

        inline Vector operator*(const Vector& v) const;

        UpperTriangularMatrix transpose() const {
            UpperTriangularMatrix trans(nlin());
            for(size_t i=0;i<size();i++) 
                trans.data()[i]=this->data()[i];
            return trans;
        }

        inline LowerTriangularMatrix inverse() const ;

        inline void info() const;
    };

    Vector UpperTriangularMatrix::operator*(const Vector& v) const {
        assert(v.nlin()==nlin());
    #ifdef HAVE_LAPACK
        Vector y(v);
        DTPMV(CblasUpper,CblasNoTrans,CblasNonUnit,(int)nlin(),data(),y.data(),1);
    #else
        Vector y(nlin());
        for (size_t i=0;i<nlin();i++) {
            y(i)=0;
            for (size_t j=i;j<nlin();j++)
                y(i)+=(*this)(i,j)*v(j);
        }
    #endif
        return y;
    }

    Vector LowerTriangularMatrix::operator*(const Vector& v) const {
        assert(v.nlin()==nlin());
    #ifdef HAVE_LAPACK
        Vector y(v);
        DTPMV(CblasUpper,CblasTrans,CblasNonUnit,(int)nlin(),data(),y.data(),1);
    #else
        Vector y(nlin());
        for (size_t i=0;i<nlin();i++) {
            y(i)=0;
            for (size_t j=0;j<=i;j++)
                y(i)+=(*this)(i,j)*v(j);
        }
    #endif
        return y;
    }

    UpperTriangularMatrix UpperTriangularMatrix::operator*(const DiagMatrix& M) const {
        assert(M.nlin()==nlin());
    // #ifdef HAVE_LAPACK
        // UpperTriangularMatrix Y(M);
        // 
    // #else
        UpperTriangularMatrix Y(*this,DEEP_COPY);
        for (size_t i=0;i<nlin();i++) {
            for (size_t j=i;j<nlin();j++)
                Y(i,j)*=M(j);
        }
    // #endif
        return Y;
    }

    SymMatrix UpperTriangularMatrix::operator*(const LowerTriangularMatrix& L) const {
        assert(L.nlin()==nlin());
    #ifdef HAVE_LAPACK
        Matrix A(*this);
        Matrix B(L);
        DTRMM(CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,(int)nlin(),(int)nlin(),1.,A.data(),(int)nlin(),B.data(),(int)nlin());
        SymMatrix S(B);
    #else
        SymMatrix S(nlin());
        S.set(0.0);
        for (size_t i=0;i<nlin();i++)
            for (size_t j=i;j<nlin();j++) {
                for (size_t k=j;k<nlin();k++)
                    S(i,j)+=(*this)(i,k)*L(k,j);
            }
    #endif
        return S;
    }

    LowerTriangularMatrix UpperTriangularMatrix::transpose() const{
        LowerTriangularMatrix trans(nlin());
        for(size_t i=0;i<size();i++) 
            trans.data()[i]=this->data()[i];
        return trans;
    }

    // inverse
    inline LowerTriangularMatrix LowerTriangularMatrix::inverse() const {
        UpperTriangularMatrix invAt=this->transpose();
        return invAt.inverse().transpose();
    #if 0  // with LAPACK we need to change the storage type, that's why it is computed as above
    #endif
    }

    // inverse
    inline UpperTriangularMatrix UpperTriangularMatrix::inverse() const {
    #ifdef HAVE_LAPACK
        UpperTriangularMatrix invA(*this,DEEP_COPY);
        int Info;
        int sz=(int)invA.nlin()*64;
        DTPTRI('U','N',invA.nlin(),invA.data(),sz,Info);
        return invA;
    #else
        std::cerr << "!!!!! Inverse not implemented !!!!!" << std::endl;
        exit(1);
    #endif
    }

    void UpperTriangularMatrix::info() const {
        if (size() == 0) {
            std::cout << "TriangularMatrix Empty" << std::endl;
            return;
        }

        std::cout << "Size : " << nlin() << std::endl;

        double maxv = this->operator()(0,0);
        size_t maxi = 0;
        size_t maxj = 0;

        for(size_t i = 0; i < nlin(); ++i)
        {
            for(size_t j = i; j < ncol(); ++j)
            {
               if (maxv < this->operator()(i,j)) {
                    maxv = this->operator()(i,j);
                    maxi = i;
                    maxj = j;
                }
            }
        }
        std::cout << "Min Value : " << "0." << " (" << "1" << "," << "0" << ")" << std::endl;
        std::cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << std::endl;
        std::cout << "First Values" << std::endl;
        for(size_t i = 0; i < std::min(nlin(),(size_t) 5); ++i)
        {
            for(size_t j = 0; j < std::min(ncol(),(size_t) 5); ++j)
            {
                std::cout << this->operator()(i,j) << " " ;
            }
            std::cout << std::endl ;
        }
    }


    void LowerTriangularMatrix::info() const {
        if (size() == 0) {
            std::cout << "TriangularMatrix Empty" << std::endl;
            return;
        }

        std::cout << "Size : " << nlin() << std::endl;

        double maxv = this->operator()(0,0);
        size_t maxi = 0;
        size_t maxj = 0;

        for(size_t i = 0; i < nlin(); ++i)
        {
            for(size_t j = 0; j <= i; ++j)
            {
               if (maxv < this->operator()(i,j)) {
                    maxv = this->operator()(i,j);
                    maxi = i;
                    maxj = j;
                }
            }
        }
        std::cout << "Min Value : " << "0." << " (" << "0" << "," << "1" << ")" << std::endl;
        std::cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << std::endl;
        std::cout << "First Values" << std::endl;
        for(size_t i = 0; i < std::min(nlin(),(size_t) 5); ++i)
        {
            for(size_t j = 0; j < std::min(ncol(),(size_t) 5); ++j)
            {
                std::cout << this->operator()(i,j) << " " ;
            }
            std::cout << std::endl ;
        }
    }
}
#endif  //! OPENMEEG_TRIANGULAR_MATRIX_H


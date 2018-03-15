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

#include <sstream>
#include <cmath>
#include <cstdlib>

#include "OpenMEEGMathsConfig.h"
#include "matrix.h"
#include "symmatrix.h"
#include "om_utils.h"

namespace OpenMEEG {

    const SymMatrix& SymMatrix::operator=(const double d) {
        for(size_t i=0;i<size();i++) data()[i]=d;
        return *this;
    }

    SymMatrix::SymMatrix(const Vector& v) {
        size_t N = v.size();
        nlin() = (size_t)((sqrt((double)(1+8*N))-1)/2+0.1);
        om_assert(nlin()*(nlin()+1)/2==N);
        value = v.value;
    }

    SymMatrix::SymMatrix(const Matrix& M): LinOp(M.nlin(),M.nlin(),SYMMETRIC,2),value(new LinOpValue(size())){
        om_assert(nlin() == M.nlin());
        for (size_t i=0; i<nlin();++i)
            for (size_t j=i; j<nlin();++j)
                (*this)(i,j) = M(i,j);
    }

    void SymMatrix::set(double x) {
        for (size_t i=0;i<(nlin()*(nlin()+1))/2;i++)
            data()[i]=x;
    }

    SymMatrix SymMatrix::operator *(double x) const {
        SymMatrix C(nlin());
        for (size_t k=0; k<nlin()*(nlin()+1)/2; k++) C.data()[k] = data()[k]*x;
        return C;
    }

    void SymMatrix::operator *=(double x) {
        for (size_t k=0; k<nlin()*(nlin()+1)/2; k++) data()[k] *= x;
    }

    Matrix SymMatrix::operator()(size_t i_start, size_t i_end, size_t j_start, size_t j_end) const {
        Matrix retMat(i_end-i_start+1,j_end-j_start+1);
        for(size_t i=0;i<=i_end-i_start;i++)
            for(size_t j=0;j<=j_end-j_start;j++)
                retMat(i,j)=this->operator()(i_start+i,j_start+j);

        return retMat;
    }

    Matrix SymMatrix::submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const {
        om_assert ( istart+isize<=nlin() && jstart+jsize<=nlin() );
        return (*this)(istart,istart+isize-1,jstart,jstart+jsize-1);
    }

    SymMatrix SymMatrix::submat(size_t istart, size_t iend) const {
        om_assert( iend > istart);
        size_t isize = iend - istart + 1;
        om_assert ( istart+isize<=nlin() );

        SymMatrix mat(isize);
        for(size_t i=istart;i<=iend;i++)
            for(size_t j=i;j<=iend;j++)
                mat(i,j)=this->operator()(i,j);

        return mat;
    }

    SymMatrix SymMatrix::operator*(const SymMatrix &m) const
    {
        om_assert(nlin()==m.nlin());
    #ifdef HAVE_BLAS
        Matrix D(*this);
        Matrix B(m);
        Matrix C(nlin(),nlin());
        DSYMM(CblasLeft,CblasUpper,sizet_to_int(nlin()),sizet_to_int(B.ncol()),1.,D.data(),sizet_to_int(D.ncol()),B.data(),sizet_to_int(B.nlin()),0,C.data(),sizet_to_int(C.nlin()));
        return SymMatrix(C);
    #else
        SymMatrix C(nlin());
        for ( size_t j = 0; j < m.ncol(); ++j) {
            for ( size_t i = 0; i < ncol(); ++i) {
                C(i, j) = 0;
                for ( size_t k = 0; k < ncol(); ++k) {
                    C(i, j) += (*this)(i, k) * m(k, j);
                }
            }
        }
        return C;
    #endif
    }

    Matrix SymMatrix::operator*(const Matrix &B) const
    {
        om_assert(ncol()==B.nlin());
        Matrix C(nlin(),B.ncol());
    #ifdef HAVE_BLAS
        Matrix D(*this);
        DSYMM(CblasLeft,CblasUpper, sizet_to_int(nlin()), sizet_to_int(B.ncol()), 1. , D.data(), sizet_to_int(D.ncol()), B.data(), sizet_to_int(B.nlin()), 0, C.data(),sizet_to_int(C.nlin()));
    #else
        for ( size_t j = 0; j < B.ncol(); ++j) {
            for ( size_t i = 0; i < ncol(); ++i) {
                C(i, j) = 0;
                for ( size_t k = 0; k < ncol(); ++k) {
                    C(i, j) += (*this)(i, k) * B(k, j);
                }
            }
        }
    #endif
        return C;
    }

    Matrix SymMatrix::solveLin(Matrix &RHS) const
    {
    #ifdef HAVE_LAPACK
        SymMatrix A(*this,DEEP_COPY);
        // LU
        BLAS_INT *pivots = new BLAS_INT[nlin()];
        int Info = 0;
        DSPTRF('U',sizet_to_int(A.nlin()),A.data(),pivots,Info);
        // Solve the linear system AX=B
        DSPTRS('U',sizet_to_int(A.nlin()),sizet_to_int(RHS.ncol()),A.data(),pivots,RHS.data(),sizet_to_int(A.nlin()),Info);
        om_assert(Info == 0);
        return RHS;
    #else
        std::cerr << "!!!!! solveLin not defined : Try a GMres !!!!!" << std::endl;
        exit(1);
    #endif
    }

    void SymMatrix::info() const {
        if (nlin() == 0) {
            std::cout << "Matrix Empty" << std::endl;
            return;
        }

        std::cout << "Dimensions : " << nlin() << " x " << ncol() << std::endl;

        double minv = this->operator()(0,0);
        double maxv = this->operator()(0,0);
        size_t mini = 0;
        size_t maxi = 0;
        size_t minj = 0;
        size_t maxj = 0;

        for (size_t i=0;i<nlin();++i)
            for (size_t j=i;j<ncol();++j)
                if (minv>this->operator()(i,j)) {
                    minv = this->operator()(i,j);
                    mini = i;
                    minj = j;
                } else if (maxv<this->operator()(i,j)) {
                    maxv = this->operator()(i,j);
                    maxi = i;
                    maxj = j;
                }

        std::cout << "Min Value : " << minv << " (" << mini << "," << minj << ")" << std::endl;
        std::cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << std::endl;
        std::cout << "First Values" << std::endl;

        for (size_t i=0;i<std::min(nlin(),(size_t) 5);++i) {
            for (size_t j=i;j<std::min(ncol(),(size_t) 5);++j)
                std::cout << this->operator()(i,j) << " " ;
            std::cout << std::endl ;
        }
    }

    // =======
    // = IOs =
    // =======

    void SymMatrix::load(const char *filename) {
        maths::ifstream ifs(filename);
        try {
            ifs >> maths::format(filename,maths::format::FromSuffix) >> *this;
        }
        catch (maths::Exception& e) {
            ifs >> *this;
        }
    }

    void SymMatrix::save(const char *filename) const {
        maths::ofstream ofs(filename);
        try {
            ofs << maths::format(filename,maths::format::FromSuffix) << *this;
        }
        catch (maths::Exception& e) {
            ofs << *this;
        }
    }
}

/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
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

#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <cfloat>

#include "om_utils.h"
#include "matrix.h"
#include "symmatrix.h"
#include "sparse_matrix.h"
#include "vector.h"

namespace OpenMEEG {

    const Matrix& Matrix::set(const double d) {
        for(size_t i=0;i<size();i++) data()[i]=d;
        return *this;
    }

    Matrix::Matrix(const SymMatrix& A): LinOp(A.nlin(),A.ncol(),FULL,TWO),value(new LinOpValue(size())) {
        for (size_t j=0; j<ncol();++j)
            for (size_t i=0; i<nlin();++i)
                (*this)(i,j) = A(i,j);
    }

    Matrix::Matrix(const Vector& v,const size_t M,const size_t N): LinOp(M,N,FULL,TWO) {
        assert(M*N==v.size());
        value = v.value;
    }

    // Matrix Matrix::pinverse(double tolrel) const {
    // #if defined(HAVE_BLAS) && defined(HAVE_LAPACK)
    //     if(ncol() > nlin()) return transpose().pinverse().transpose();
    //     else {
    //         Matrix result(ncol(),nlin());
    //         Matrix U,S,V;
    //         svd(U,S,V);
    //         double maxs=0;
    //         int mimi=(int)std::min(S.nlin(),S.ncol());
    //         for(int i=0;i<mimi;i++) maxs=std::max(S(i,i),maxs);
    //         if (tolrel==0) tolrel=DBL_EPSILON;
    //         double tol = std::max(nlin(),ncol()) * maxs * tolrel;
    //         int r=0; for(int i=0;i<mimi;i++) if(S(i,i)>tol) r++;
    //         if (r == 0) {
    //             result.set(0.);
    //             return result;
    //         } else {
    //             Matrix s(r,r); s.set(0);
    //             for(int i=0;i<r;i++) s(i,i)=1.0/S(i,i);
    //             const Matrix Vbis(V,r);
    //             const Matrix Ubis(U,r);
    //             return Vbis*s*Ubis.transpose();
    //         }
    //     }
    // #else
    //     std::cerr << "pinv not implemented without blas/lapack" << std::endl;
    //     exit(1);
    // #endif
    // }

    Matrix Matrix::transpose() const {
        Matrix result(ncol(),nlin());
        for(size_t i=0;i<nlin();i++) for(size_t j=0;j<ncol();j++) result(j,i)=(*this)(i,j);
        return result;
    }

    // void Matrix::svd(Matrix &U,Matrix &S, Matrix &V) const {
    // #ifdef HAVE_LAPACK
    //     Matrix cpy(*this,DEEP_COPY);
    //     int mimi = (int)std::min(nlin(),ncol());
    //     U = Matrix(nlin(),ncol()); U.set(0);
    //     V = Matrix(ncol(),ncol()); V.set(0);
    //     S = Matrix(ncol(),ncol()); S.set(0);
    //     double *s=new double[mimi];
    //     int lwork=4 *mimi*mimi + (int)std::max(nlin(),ncol()) + 9*mimi;
    //     double *work=new double[lwork];
    //     int *iwork=new int[8*mimi];
    //     int info;
    //     DGESDD('S',nlin(),ncol(),cpy.data(),nlin(),s,U.data(),U.nlin(),V.data(),V.nlin(),work,lwork,iwork,info);
    //     for(int i=0;i<mimi;i++) S(i,i)=s[i];
    //     V=V.transpose();
    //     delete[] s;
    //     delete[] work;
    //     delete[] iwork;
    // #else
    //     std::cerr<<"svd not implemented without blas/lapack"<<std::endl;
    // #endif
    // }

    Matrix Matrix::operator *(const SparseMatrix &mat) const
    {
        assert(ncol()==mat.nlin());
        Matrix out(nlin(),mat.ncol());
        out.set(0.0);

        SparseMatrix::const_iterator it;
        for(it = mat.begin(); it != mat.end(); ++it) {
            size_t i = it->first.first;
            size_t j = it->first.second;
            double val = it->second;
            for(size_t k = 0; k < nlin(); ++k) {
                out(k,j) += this->operator()(k,i) * val;
            }
        }
        return out;
    }

    Matrix Matrix::operator*(double x) const {
        Matrix C(nlin(),ncol());
        for (size_t k=0; k<nlin()*ncol(); k++) C.data()[k] = data()[k]*x;
        return C;
    }

    Matrix Matrix::operator/(double x) const {
        Matrix C(nlin(),ncol());
        for (size_t k=0; k<nlin()*ncol(); k++) C.data()[k] = data()[k]/x;
        return C;
    }

    void Matrix::operator*=(double x) {
        for (size_t k=0; k<nlin()*ncol(); k++) data()[k] *= x;
    }

    void Matrix::operator/=(double x) {
        for (size_t k=0; k<nlin()*ncol(); k++) data()[k] /= x;
    }

    Vector Matrix::mean() const {
        Vector v(ncol()); v.set(0);
        for(size_t j = 0; j < ncol(); ++j) {
            for(size_t i = 0; i < nlin(); ++i) {
                v(j) += this->operator()(i,j);
            }
        }
        for(size_t j = 0; j < ncol(); ++j) {
            v(j) = v(j) / nlin();
        }
        return v;
    }

    Vector Matrix::tmean() const {
        Vector v(nlin()); v.set(0);
        for(size_t j = 0; j < ncol(); ++j) {
            for(size_t i = 0; i < nlin(); ++i) {
                v(i) += this->operator()(i,j);
            }
        }
        for(size_t i = 0; i < nlin(); ++i) {
            v(i) = v(i) / ncol();
        }
        return v;
    }

    void Matrix::info() const {
        if ((nlin() == 0) && (ncol() == 0)) {
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

        for(size_t i = 0; i < nlin(); ++i)
        {
            for(size_t j = 0; j < ncol(); ++j)
            {
                if (minv > this->operator()(i,j)) {
                    minv = this->operator()(i,j);
                    mini = i;
                    minj = j;
                } else if (maxv < this->operator()(i,j)) {
                    maxv = this->operator()(i,j);
                    maxi = i;
                    maxj = j;
                }
            }
        }
        std::cout << "Min Value : " << minv << " (" << mini << "," << minj << ")" << std::endl;
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

    // =======
    // = IOs =
    // =======

    void Matrix::loadBin( const char *filename )
    {
        maths::ifstream ifs(filename);
        ifs >> maths::format("binary") >> *this;
    }

    void Matrix::saveBin( const char *filename ) const
    {
        maths::ofstream ofs(filename);
        ofs << maths::format("binary") << *this;
    }

    void Matrix::loadTxt( const char *filename )
    {
        maths::ifstream ifs(filename);
        ifs >> maths::format("ascii") >> *this;
    }

    void Matrix::saveTxt( const char *filename ) const
    {
        maths::ofstream ofs(filename);
        ofs << maths::format("ascii") << *this;
    }

    void Matrix::loadMat(const char *filename)
    {
        maths::ifstream ifs(filename);
        ifs >> maths::format("matlab") >> *this;
    }

    void Matrix::saveMat( const char *filename ) const
    {
        maths::ofstream ofs(filename);
        ofs << maths::format("matlab") << *this;
    }

    void Matrix::loadBrainvisa(const char *filename)
    {
        maths::ifstream ifs(filename);
        ifs >> maths::format("tex") >> *this;
    }

    void Matrix::saveBrainvisa( const char *filename ) const
    {
        maths::ofstream ofs(filename);
        ofs << maths::format("tex") << *this;
    }

    void Matrix::load( const char *filename ) {
        try {
            maths::ifstream ifs(filename);
            ifs >> *this;
        }
        catch (std::string s) {
            std::cout << s << std::endl;
        }
    }

    void Matrix::save( const char *filename ) const {
        try {
            maths::ofstream ofs(filename);
            ofs << *this;
        }
        catch (std::string s) {
            std::cout << s << std::endl;
        }
    }
}

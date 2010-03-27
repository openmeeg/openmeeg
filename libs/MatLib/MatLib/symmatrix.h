/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

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

#ifndef SYMMATRIX_H
#define SYMMATRIX_H

#include <fstream>
#include <cassert>

#include "vector.h"
#include "linop.h"

class Matrix;

class OPENMEEGMATHS_EXPORT SymMatrix : public LinOp {

    friend class Vector;

    utils::RCPtr<LinOpValue> value;

    std::string identity() const;

public:

    SymMatrix(): LinOp(0,0,SYMMETRIC,TWO),value() {}

    SymMatrix(const char* fname): LinOp(0,0,SYMMETRIC,TWO),value() { this->load(fname); }
    SymMatrix(size_t N): LinOp(N,N,SYMMETRIC,TWO),value(new LinOpValue(size())) { }
    SymMatrix(const SymMatrix& S,const DeepCopy): LinOp(S.nlin(),S.nlin(),SYMMETRIC,TWO),value(new LinOpValue(S.size(),S.data())) { }

    explicit SymMatrix(const Vector& v);

    size_t size() const { return nlin()*(nlin()+1)/2; };
    void info() const ;

    size_t  ncol() const { return nlin(); } // SymMatrix only need num_lines
    size_t& ncol()       { return nlin(); }

    void alloc_data() { value = new LinOpValue(size()); }

    bool empty() const { return value->empty(); }
    void set(double x) ;
    double* data() const { return value->data; }

    inline double operator()(size_t i,size_t j) const;
    inline double& operator()(size_t i,size_t j) ;

    Matrix operator()(size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;
    Matrix submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const;
    SymMatrix submat(size_t istart, size_t iend) const;
    Vector solveLin(const Vector &B) const;
    void solveLin(Vector * B, int nbvect);

    const SymMatrix& operator=(const double d);

    SymMatrix operator+(const SymMatrix& B) const;
    SymMatrix operator-(const SymMatrix& B) const;
    SymMatrix operator*(double x) const;
    SymMatrix operator/(double x) const {return (*this)*(1/x);}
    void operator +=(const SymMatrix& B);
    void operator -=(const SymMatrix& B);
    void operator *=(double x);
    void operator /=(double x) { (*this)*=(1/x); }
    Matrix operator*(const Matrix& B) const; // faux !!
    Vector operator*(const Vector& v) const; // faux ?

    SymMatrix inverse() const;
    SymMatrix posdefinverse() const;
    double det();
    void eigen(Matrix & Z, Vector & D );

    void save( const char *filename ) const;
    void load( const char *filename );
    void saveTxt( const char *filename ) const;
    void saveBin( const char *filename ) const;
    void loadTxt( const char *filename );
    void loadBin( const char *filename );

    friend class Matrix;
};

inline double SymMatrix::operator()(size_t i,size_t j) const {
    assert(i<nlin() && j<nlin());
    if(i<=j)
        return data()[i+j*(j+1)/2];
    else
        return data()[j+i*(i+1)/2];
}

inline double& SymMatrix::operator()(size_t i,size_t j) {
    assert(i<nlin() && j<nlin());
    if(i<=j)
        return data()[i+j*(j+1)/2];
    else
        return data()[j+i*(i+1)/2];
}

#endif

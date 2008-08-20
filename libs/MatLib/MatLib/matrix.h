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

#ifndef MATRIX_H
#define MATRIX_H

#include "MatLibConfig.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>

#ifdef USE_MATIO
#include <matio.h>
#endif

#include "linop.h"
#include "MathsIO.H"

class SparseMatrix;
class SymMatrix;
class Vector;

class OPENMEEGMATHS_EXPORT Matrix: public LinOp {
protected:
    double *t;
    int* count;

    std::string identity() const;
    void copy(const Matrix& A);
    void destroy();
    void copyout(double * p) const;
    void copyin(const double * p);
public:
    Matrix();
    Matrix(const char*);
    Matrix(size_t M,size_t N);
    Matrix(const Matrix& A);
    explicit Matrix(const SymMatrix& A);
    Matrix(const Vector& v, size_t M, size_t N); // violent access to v.t and v.count
    Matrix(double* T, int* COUNT, size_t M, size_t N);

    void alloc_data();

    ~Matrix() { destroy(); }
    bool empty() const;
    size_t size() const { return nlin()*ncol(); };
    void DangerousBuild( double *, size_t i, size_t j);
    void DangerousKill ();
    double* data() const;
    int* DangerousGetCount () const;

    inline double operator()(size_t i,size_t j) const ;
    inline double& operator()(size_t i,size_t j) ;

    Matrix duplicate() const;
    void copyin(const Matrix& A) ;

    Matrix submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const;
    Vector getcol(size_t j) const;
    void setcol(size_t j, const Vector& v);
    Vector getlin(size_t i) const;
    void setlin(size_t i, const Vector& v);

    // inline double operator[](size_t i) const ;
    // inline double& operator[](size_t i) ;

    const Matrix& operator=(const Matrix& A);
    const Matrix& set(const double d);

    Matrix operator*(const Matrix& B) const;
    Matrix operator*(const SymMatrix& B) const;
    Matrix operator*(const SparseMatrix& B) const;
    Matrix operator+(const Matrix& B) const;
    Matrix operator-(const Matrix& B) const;
    Matrix operator*(double x) const;
    Matrix operator/(double x) const;
    void operator+=(const Matrix& B);
    void operator-=(const Matrix& B);
    void operator*=(double x);
    void operator/=(double x);

    Vector operator*(const Vector& v) const;
    Vector tmult(const Vector &v) const;
    Matrix tmult(const Matrix &m) const;
    Matrix multt(const Matrix &m) const;
    Matrix tmultt(const Matrix &m) const;

    Vector mean() const;
    Vector tmean() const;

    Matrix transpose () const;
    Matrix inverse() const;
    Matrix pinverse(double reltol=0) const;
    void svd(Matrix &U,Matrix &S, Matrix &V) const;

    double frobenius_norm() const;
    double dot(const Matrix& B) const;

    static void readDimsBin( const char* filename, size_t& mm, size_t& nn)
    {
        FILE *infile=fopen(filename,"rb");
        if(infile == NULL) {
            std::cerr<<"Error Opening Matrix File "<<filename<<std::endl;
            exit(1);
        }
        unsigned int ui;
        fread(&ui,sizeof(unsigned int),1,infile);
        mm=ui;
        fread(&ui,sizeof(unsigned int),1,infile);
        nn=ui;
        fclose(infile);
    }

    void save( const char *filename ) const;
    void saveTxt( const char *filename ) const;
    void saveBin( const char *filename ) const;
    void saveMat( const char *filename ) const;

    void load( const char *filename );
    void loadTxt( const char *filename );
    void loadBin( const char *filename );
    void loadMat( const char *filename );

    void info() const;

    friend class SparseMatrix;
    friend class SymMatrix;
};

inline double Matrix::operator()(size_t i,size_t j) const
{
    assert(i<nlin() && j<ncol());
    return t[i+nlin()*j];
}
inline double& Matrix::operator()(size_t i,size_t j)
{
    assert(i<nlin() && j<ncol());
    return t[i+nlin()*j];
}

#endif


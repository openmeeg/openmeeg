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

/** \brief  Matrix class

    Matrix class
**/
class OPENMEEGMATHS_EXPORT Matrix: public LinOp {
protected:

    friend class Vector;

    utils::RCPtr<LinOpValue> value;

    std::string identity() const;

    explicit Matrix(const Matrix& A,const size_t M): LinOp(A.nlin(),M,FULL,TWO),value(A.value) { }

public:

    Matrix(): LinOp(0,0,FULL,TWO),value() { }
    Matrix(const char* fname): LinOp(0,0,FULL,TWO),value() { this->load(fname); }
    Matrix(size_t M,size_t N): LinOp(M,N,FULL,TWO),value(new LinOpValue(N*M)) { }
    Matrix(const Matrix& A,const DeepCopy): LinOp(A.nlin(),A.ncol(),FULL,TWO),value(new LinOpValue(A.size(),A.data())) { }

    explicit Matrix(const SymMatrix& A);
    Matrix(const Vector& v,const size_t M,const size_t N);

    void alloc_data() { value = new LinOpValue(size()); }

    /** \brief Test if Matrix is empty
        \return true if Matrix is empty
        \sa
    **/
    bool empty() const { return value->empty(); }

    /** \brief Get Matrix size
        \return number of values (nb lines x nb columns)
        \sa
    **/
    size_t size() const { return nlin()*ncol(); };

    /** \brief Get Matrix data
        \return pointer to Matrix values
        \sa
    **/
    double* data() const { return value->data; }

    /** \brief Get Matrix value
        \return value in Matrix
        \sa
    **/
    inline double operator()(size_t i,size_t j) const ;

    /** \brief Get Matrix value
        \return reference to value in Matrix
        \sa
    **/
    inline double& operator()(size_t i,size_t j) ;

    Matrix submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const;
    Vector getcol(size_t j) const;
    void setcol(size_t j, const Vector& v);
    Vector getlin(size_t i) const;
    void setlin(size_t i, const Vector& v);

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

    /** \brief Get Matrix Frobenius norm
        \return norm value
        \sa
    **/
    double frobenius_norm() const;
    double dot(const Matrix& B) const;

    /** \brief Read Matrix dimensions for raw binary file without loading the full data
        \sa
    **/
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

    /** \brief Save Matrix to file (Format set using file name extension)
        \sa
    **/
    void save( const char *filename ) const;

    /** \brief Save Matrix in ascii file
        \sa
    **/
    void saveTxt( const char *filename ) const;

    /** \brief Save Matrix in raw binary file
        \sa
    **/
    void saveBin( const char *filename ) const;

    /** \brief Save Matrix in Matlab file
        \sa
    **/
    void saveMat( const char *filename ) const;

    /** \brief Save Matrix in Brainvisa ascii texture file
        \sa
    **/
    void saveBrainvisa( const char *filename ) const;

    /** \brief Load Matrix from file (Format set using file name extension)
        \sa
    **/
    void load( const char *filename );

    /** \brief Load Matrix from ascii file
        \sa
    **/
    void loadTxt( const char *filename );

    /** \brief Load Matrix from raw binary file
        \sa
    **/
    void loadBin( const char *filename );

    /** \brief Load Matrix from Matlab file
        \sa
    **/
    void loadMat( const char *filename );

    /** \brief Load Matrix from Brainvisa ascii texture file
        \sa
    **/
    void loadBrainvisa( const char *filename );

    /** \brief Print info on Matrix
        \sa
    **/
    void info() const;

    friend class SparseMatrix;
    friend class SymMatrix;
};

inline double Matrix::operator()(size_t i,size_t j) const
{
    assert(i<nlin() && j<ncol());
    return value->data[i+nlin()*j];
}
inline double& Matrix::operator()(size_t i,size_t j)
{
    assert(i<nlin() && j<ncol());
    return value->data[i+nlin()*j];
}

#endif

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

#ifndef VECTOR_H
#define VECTOR_H

#include "MatLibConfig.h"
#include "linop.h"
#include "MathsIO.H"

#include <iostream>
#include <cassert>

class Matrix;
class SymMatrix;

class Vector: public LinOp {
    double *t;
    int *count;

    void destroy();
    void copy(const Vector& A) ;
    void copyout(double * p) const;
    void copyin(const double * p);
public:
    Vector();
    Vector(size_t N);
    Vector(const Vector& A);
    explicit Vector(Matrix& A);
    explicit Vector(SymMatrix& A);
    Vector(double* T, int* COUNT, size_t N);
     ~Vector() { destroy(); }

    void alloc_data();

    size_t size() const ;
    bool empty() const ;
    void DangerousBuild(double *t, size_t size);
    void DangerousKill();
    double* data() const ;
    int* DangerousGetCount() const ;

    Vector duplicate() const ;
    void copyin(const Vector& A) ;

    inline double operator()(size_t i) const ;
    inline double& operator()(size_t i) ;

    const Vector& operator=(const Vector& A);

    Vector operator+(const Vector& v) const;
    Vector operator-(const Vector& v) const;
    void operator+=(const Vector& v);
    void operator-=(const Vector& v);
    void operator*=(double x);
    void operator/=(double x);
    Vector operator+(double i) const;
    Vector operator-(double i) const;
    Vector operator*(double x) const;
    Vector operator/(double x) const ;
    double operator*(const Vector& v) const;

    Vector kmult(const Vector& x) const;
    Vector conv(const Vector& v) const;
    Vector conv_trunc(const Vector& v) const;
    Matrix outer_product(const Vector& v) const;

    double norm() const;
    double sum() const;
    double mean() const;

    void set(double x);
    void saveTxt( const char* filename) const;
    void saveBin( const char* filename) const;
    void loadTxt( const char* filename);
    void loadBin( const char* filename);
    void save( const char *filename ) const ;
    void load( const char *filename ) ;

    void info() const;

    friend class SymMatrix;
    friend class Matrix;
};

Vector operator * (const double &d, const Vector &v) ;
std::ostream& operator<<(std::ostream& f,const Vector &M);
std::istream& operator>>(std::istream& f,Vector &M);

inline double Vector::operator()(size_t i) const {
    assert(i<nlin());
    return t[i];
}
inline double& Vector::operator()(size_t i) {
    assert(i<nlin());
    return t[i];
}

#endif


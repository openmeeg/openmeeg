/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
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

#ifndef H_VECTEUR_DCL
#define H_VECTEUR_DCL

#include "MatLibConfig.h"
#include "base_matrix.H"
#include "MatrixIO.H"

#include <iostream>
#include <cassert>

class matrice;
class symmatrice;

class vecteur: public MatrixBase {
    double *t;
    int *count;

    void destroy();
    void copy(const vecteur& A) ;
    void copyout(double * p) const;
    void copyin(const double * p);
public:
    vecteur();
    vecteur(size_t N);
    vecteur(const vecteur& A);
    explicit vecteur(matrice& A);
    explicit vecteur(symmatrice& A);
    vecteur(double* T, int* COUNT, size_t N);
     ~vecteur() { destroy(); }

    void alloc_data();

    size_t size() const ;
    bool empty() const ;
    void DangerousBuild(double *t, size_t size);
    void DangerousKill();
    double* data() const ;
    int* DangerousGetCount() const ;

    vecteur duplicate() const ;
    void copyin(const vecteur& A) ;

    inline double operator()(size_t i) const ;
    inline double& operator()(size_t i) ;

    const vecteur& operator=(const vecteur& A);

    vecteur operator+(const vecteur& v) const;
    vecteur operator-(const vecteur& v) const;
    void operator+=(const vecteur& v);
    void operator-=(const vecteur& v);
    void operator*=(double x);
    void operator/=(double x);
    vecteur operator+(double i) const;
    vecteur operator-(double i) const;
    vecteur operator*(double x) const;
    vecteur operator/(double x) const ;
    double operator*(const vecteur& v) const;

    vecteur kmult(const vecteur& x) const;
    vecteur conv(const vecteur& v) const;
    vecteur conv_trunc(const vecteur& v) const;
    matrice outer_product(const vecteur& v) const;

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

    friend class symmatrice;
    friend class matrice;
};

vecteur operator * (const double &d, const vecteur &v) ;
std::ostream& operator<<(std::ostream& f,const vecteur &M);
std::istream& operator>>(std::istream& f,vecteur &M);

inline double vecteur::operator()(size_t i) const {
    assert(i<nlin());
    return t[i];
}
inline double& vecteur::operator()(size_t i) {
    assert(i<nlin());
    return t[i];
}

#endif


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

#ifndef H_SYMMATRICE_DCL
#define H_SYMMATRICE_DCL

#include <fstream>
#include <cassert>

#include "vecteur.h"
#include "base_matrix.H"

class matrice;

class symmatrice : public MatrixBase {
    double *t;
    int *count;

    std::string identity() const;
    void copy(const symmatrice& A);
    void destroy();
public:
    symmatrice();
    symmatrice(const char* fname);
    symmatrice(size_t N) ;
    symmatrice(const symmatrice& A);
    explicit symmatrice(const vecteur& v);
    explicit symmatrice(const matrice& A); // upper triangle based (lower triangle need not to be set)
    symmatrice(double* T, int* COUNT, size_t N);
     ~symmatrice() { destroy(); }
    symmatrice duplicate() const;
    size_t size() const { return nlin()*(nlin()+1)/2; };
    void info() const ;

    void alloc_data();

    bool empty() const ;
    void set(double x) ;
    double* data() const ;
    int* DangerousGetCount() ;

    inline double operator()(size_t i,size_t j) const;
    inline double& operator()(size_t i,size_t j) ;

    matrice operator()(size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;
    matrice submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const;
    symmatrice submat(size_t istart, size_t iend) const;
    vecteur solveLin(const vecteur &B) const;
    void solveLin(vecteur * B, int nbvect);

    const symmatrice& operator=(const double d);
    const symmatrice& operator=(const symmatrice& A);

    symmatrice operator+(const symmatrice& B) const;
    symmatrice operator-(const symmatrice& B) const;
    symmatrice operator*(double x) const;
    symmatrice operator/(double x) const {return (*this)*(1/x);}
    void operator +=(const symmatrice& B);
    void operator -=(const symmatrice& B);
    void operator *=(double x);
    void operator /=(double x) ;
    matrice operator*(const matrice& B) const; // faux !!
    vecteur operator*(const vecteur& v) const; // faux ?

    symmatrice inverse() const;
    symmatrice posdefinverse() const;
    double det();
    void eigen(matrice & Z, vecteur & D );

    void save( const char *filename ) const;
    void load( const char *filename );
    void saveTxt( const char *filename ) const;
    void saveBin( const char *filename ) const;
    void loadTxt( const char *filename );
    void loadBin( const char *filename );

    friend class matrice;
};

inline double symmatrice::operator()(size_t i,size_t j) const {
    assert(i<nlin() && j<nlin());
    if(i<=j)
        return t[i+j*(j+1)/2];
    else
        return t[j+i*(i+1)/2];
}

inline double& symmatrice::operator()(size_t i,size_t j) {
    assert(i<nlin() && j<nlin());
    if(i<=j)
        return t[i+j*(j+1)/2];
    else
        return t[j+i*(i+1)/2];
}

std::ostream& operator<<(std::ostream& f,const symmatrice &M);

#endif


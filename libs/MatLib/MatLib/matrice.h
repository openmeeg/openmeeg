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

#ifndef H_MATRICE_DCL
#define H_MATRICE_DCL

#include "MatLibConfig.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>

#ifdef USE_MATIO
#include <matio.h>
#endif

#include "base_matrix.H"
#include "MatrixIO.H"

class sparse_matrice;
class symmatrice;
class vecteur;

class matrice: public MatrixBase {
protected:
    double *t;
    int* count;

    std::string identity() const;
    void copy(const matrice& A);
    void destroy();
    void copyout(double * p) const;
    void copyin(const double * p);
public:
    matrice();
    matrice(const char*);
    matrice(size_t M,size_t N);
    matrice(const matrice& A);
    explicit matrice(const symmatrice& A);
    matrice(const vecteur& v, size_t M, size_t N); // violent access to v.t and v.count
    matrice(double* T, int* COUNT, size_t M, size_t N);

    void alloc_data();

    ~matrice() { destroy(); }
    bool empty() const;
    size_t size() const { return nlin()*ncol(); };
    void DangerousBuild( double *, size_t i, size_t j);
    void DangerousKill ();
    double* data() const;
    int* DangerousGetCount () const;

    inline double operator()(size_t i,size_t j) const ;
    inline double& operator()(size_t i,size_t j) ;

    matrice duplicate() const;
    void copyin(const matrice& A) ;

    matrice submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const;
    vecteur getcol(size_t j) const;
    void setcol(size_t j, const vecteur& v);
    vecteur getlin(size_t i) const;
    void setlin(size_t i, const vecteur& v);

    // inline double operator[](size_t i) const ;
    // inline double& operator[](size_t i) ;

    const matrice& operator=(const matrice& A);
    const matrice& set(const double d);

    matrice operator*(const matrice& B) const;
    matrice operator*(const symmatrice& B) const;
    matrice operator*(const sparse_matrice& B) const;
    matrice operator+(const matrice& B) const;
    matrice operator-(const matrice& B) const;
    matrice operator*(double x) const;
    matrice operator/(double x) const;
    void operator+=(const matrice& B);
    void operator-=(const matrice& B);
    void operator*=(double x);
    void operator/=(double x);

    vecteur operator*(const vecteur& v) const;
    vecteur tmult(const vecteur &v) const;
    matrice tmult(const matrice &m) const;
    matrice multt(const matrice &m) const;
    matrice tmultt(const matrice &m) const;

    vecteur mean() const;
    vecteur tmean() const;

    matrice transpose () const;
    matrice inverse() const;
    matrice pinverse(double reltol=0) const;
    void svd(matrice &U,matrice &S, matrice &V) const;

    double frobenius_norm() const;
    double dot(const matrice& B) const;

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

    friend class sparse_matrice;
    friend class symmatrice;
};

inline double matrice::operator()(size_t i,size_t j) const
{
    assert(i<nlin() && j<ncol());
    return t[i+nlin()*j];
}
inline double& matrice::operator()(size_t i,size_t j)
{
    assert(i<nlin() && j<ncol());
    return t[i+nlin()*j];
}

std::ostream& operator<<(std::ostream& f,const matrice &M);

#endif


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

#ifdef USE_MATIO
#include <matio.h>
#endif

class symmatrice;
class vecteur;

class matrice {
protected:
    size_t m,n;
    double *t;
    int* count;

    inline void alloc(size_t M,size_t N);
    inline void copy(const matrice& A);
    inline void destroy();
    inline void copyout(double * p) const;
    inline void copyin(const double * p);
public:
    inline matrice();
    inline matrice(size_t M,size_t N);
    inline matrice(const matrice& A);
    inline explicit matrice(const symmatrice& A);
    inline matrice(const vecteur& v, size_t M, size_t N); // accede violemment a v.t et v.count
    inline matrice(double* T, int* COUNT, size_t M, size_t N);
    inline matrice(const char *filename, char c='t');

    inline ~matrice() { destroy(); }
    inline size_t nlin() const;
    inline size_t ncol() const;
    inline bool empty() const;
    inline void DangerousBuild( double *, size_t i, size_t j);
    inline void DangerousKill ();
    inline double* DangerousGetData () const;
    inline int* DangerousGetCount () const;
    inline void DangerousReshape(size_t M, size_t N);

    // inline void print() const { std::cout << *this; };

    inline matrice duplicate() const;
    inline void copyin(const matrice& A) ;

    inline matrice colsref(size_t jstart, size_t jsize) const ;
    inline vecteur colref(size_t j) const;
    inline vecteur subcolref(size_t j, size_t istart, size_t isize) const ;
    inline matrice getsubmat(size_t istart, size_t isize, size_t jstart, size_t jsize) const;
    inline vecteur getcol(size_t j) const;
    inline void setcol(size_t j, const vecteur& v);
    inline vecteur getlin(size_t i) const;
    inline void setlin(size_t i, const vecteur& v);

    inline double operator[](size_t i) const ;
    inline double& operator[](size_t i) ;

    inline double operator()(size_t i,size_t j) const ;
    inline double& operator()(size_t i,size_t j) ;

    inline matrice operator()(size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;

    inline const matrice& operator=(const matrice& A);
    inline const matrice& set(const double d);

    inline matrice operator*(const matrice& B) const;
    inline matrice operator*(const symmatrice& B) const;
    inline matrice operator+(const matrice& B) const;
    inline matrice operator-(const matrice& B) const;
    inline matrice operator*(double x) const;
    inline matrice operator/(double x) const;
    inline void operator+=(const matrice& B);
    inline void operator-=(const matrice& B);
    inline void operator*=(double x);
    inline void operator/=(double x);

    inline vecteur operator*(const vecteur& v) const;
    inline vecteur tmult(const vecteur &v) const;
    inline matrice tmult(const matrice &m) const;
    inline matrice multt(const matrice &m) const;
    inline matrice tmultt(const matrice &m) const;

    inline vecteur mean() const;
    inline vecteur tmean() const;
    // inline matrice operator.*(const matrice&) const;

    inline matrice transpose () const;
    inline matrice inverse() const;
    inline matrice pinverse(double reltol=0) const;
    inline void svd(matrice &U,matrice &S, matrice &V) const;
    inline double det() const;

    inline double frobenius_norm() const;
    inline double dot(const matrice& B) const;
    inline matrice lowercholesky() const;

    inline void write(std::ostream& f) const {
        f.write(reinterpret_cast<const char*>(&m),(std::streamsize)sizeof(int));
        f.write(reinterpret_cast<const char*>(&n),(std::streamsize)sizeof(int));
        f.write(reinterpret_cast<const char*>(t),(std::streamsize)(n*m*sizeof(double)));
    }

    inline void read(std::istream& f) {
        destroy();
        f.read((char*)&m,(std::streamsize)sizeof(int));
        f.read((char*)&n,(std::streamsize)sizeof(int));
        alloc(m,n);
        f.read((char*)t,(std::streamsize)(n*m*sizeof(double)));
    }

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

    virtual std::ostream& operator>>(std::ostream& f) const {
        for (size_t i=0;i<this->nlin();i++) {
            for (size_t j=0;j<this->ncol();j++) {
                f << (*this)(i,j) << " ";
            }
            f << std::endl;
        }
        return f;
    }

    inline void save( const char *filename ) const;
    inline void saveTxt( const char *filename ) const;
    inline void saveSubTxt( const char *filename , size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;
    inline void saveBin( const char *filename ) const;
    inline void saveSubBin( const char *filename , size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;
    inline void saveMat( const char *filename ) const;

    inline void load( const char *filename );
    inline void loadTxt( const char *filename );
    inline void loadBin( const char *filename );
    inline void loadMat( const char *filename );

    inline void info() const;

    friend class symmatrice;
};

#endif


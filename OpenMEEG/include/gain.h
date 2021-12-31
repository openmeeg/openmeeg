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

#pragma once

#include "matrix.h"
#include "sparse_matrix.h"
#include "symmatrix.h"
#include "geometry.h"
#include "progressbar.h"
#include "assemble.h"

#define USE_GMRES 0
#if USE_GMRES
#include "gmres.h"
#endif

namespace OpenMEEG {

#if USE_GMRES

    // Consider the GMRes solver for problems with dimension>15,000 (3,000 vertices per interface)

    template <typename MATRIX>
    Matrix linsolve(const SymMatrix& H,const MATRIX& S) {
        Matrix res(S.nlin(),H.nlin());
        Jacobi<SymMatrix> M(H);    // Jacobi preconditionner
        #pragma omp parallel for
        #ifdef OPENMP_UNSIGNED
        for (unsigned i=0; i<S.nlin(); ++i) {
        #else
        for (int i=0; i<static_cast<int>(S.nlin()); ++i) {
        #endif
            Vector vtemp(H.nlin());
            GMRes(H,M,vtemp,S.getlin(i),1000,1e-7,H.nlin()); // max number of iteration=1000, and precision=1e-7 (1e-5 for faster resolution)
            res.setlin(i,vtemp);
            #pragma omp critical
            PROGRESSBAR(i,S.nlin());
        }
        return res
    }
#else
    template <typename SelectionMatrix>
    Matrix linsolve(const SymMatrix& H,const SelectionMatrix& S) {
        Matrix res(S.transpose());
        H.solveLin(res); // solving the system AX=B with LAPACK
        return res.transpose();
    }
#endif

    class GainMEG: public Matrix {
    public:
        using Matrix::operator=;
        GainMEG (const Matrix& GainMat): Matrix(GainMat) {}
        GainMEG(const SymMatrix& HeadMatInv,const Matrix& SourceMat,const Matrix& Head2MEGMat,const Matrix& Source2MEGMat):
            Matrix(Source2MEGMat+(Head2MEGMat*HeadMatInv)*SourceMat)
        { }
    };

    class GainEEG: public Matrix {
    public:
        using Matrix::operator=;
        GainEEG (const Matrix& GainMat): Matrix(GainMat) {}
        GainEEG (const SymMatrix& HeadMatInv,const Matrix& SourceMat,const SparseMatrix& Head2EEGMat):
            Matrix((Head2EEGMat*HeadMatInv)*SourceMat)
        { }
    };

    class GainEEGadjoint: public Matrix {
    public:

        using Matrix::operator=;

        GainEEGadjoint(const Geometry& geo,const Matrix& dipoles,const SymMatrix& HeadMat,const SparseMatrix& Head2EEGMat): Matrix(Head2EEGMat.nlin(),dipoles.nlin()) {
            const Matrix& Hinv = linsolve(HeadMat,Head2EEGMat);
            ProgressBar pb(ncol());
            for (unsigned i=0; i<ncol(); ++i,++pb)
                setcol(i,Hinv*DipSourceMat(geo,dipoles.submat(i,1,0,dipoles.ncol())).getcol(0)); // TODO ugly
        }
    };

    class GainMEGadjoint: public Matrix {
    public:

        using Matrix::operator=;

        GainMEGadjoint(const Geometry& geo,const Matrix& dipoles,const SymMatrix& HeadMat,const Matrix& Head2MEGMat,const Matrix& Source2MEGMat):
            Matrix(Head2MEGMat.nlin(),dipoles.nlin()) 
        {
            const Matrix& Hinv = linsolve(HeadMat,Head2MEGMat);
            ProgressBar pb(ncol());
            for (unsigned i=0; i<ncol(); ++i,++pb)
                setcol(i,Hinv*DipSourceMat(geo,dipoles.submat(i,1,0,dipoles.ncol())).getcol(0)+Source2MEGMat.getcol(i)); // TODO ugly
        }
    };

    class GainEEGMEGadjoint {
    public:
        GainEEGMEGadjoint(const Geometry& geo,const Matrix& dipoles,const SymMatrix& HeadMat,const SparseMatrix& Head2EEGMat,const Matrix& Head2MEGMat,const Matrix& Source2MEGMat):
            EEGleadfield(Head2EEGMat.nlin(),dipoles.nlin()),MEGleadfield(Head2MEGMat.nlin(),dipoles.nlin())
        {
            Matrix RHS(Head2EEGMat.nlin()+Head2MEGMat.nlin(),HeadMat.nlin());
            for (unsigned i=0; i<Head2EEGMat.nlin(); ++i)
                RHS.setlin(i,Head2EEGMat.getlin(i));
            for (unsigned i=0; i<Head2MEGMat.nlin(); ++i)
                RHS.setlin(i+Head2EEGMat.nlin(),Head2MEGMat.getlin(i));

            const Matrix& Hinv = linsolve(HeadMat,RHS);

            ProgressBar pb(dipoles.nlin());
            for (unsigned i=0; i<dipoles.nlin(); ++i,++pb) {
                const Vector& dsm = DipSourceMat(geo,dipoles.submat(i,1,0,dipoles.ncol())).getcol(0); // TODO ugly
                EEGleadfield.setcol(i,Hinv.submat(0,Head2EEGMat.nlin(),0,HeadMat.nlin())*dsm);
                MEGleadfield.setcol(i,Hinv.submat(Head2EEGMat.nlin(),Head2MEGMat.nlin(),0,HeadMat.nlin())*dsm+Source2MEGMat.getcol(i));
            }
        }
        
        void saveEEG( const std::string filename ) const { EEGleadfield.save(filename); }
        void saveMEG( const std::string filename ) const { MEGleadfield.save(filename); }
        
    private:

        Matrix EEGleadfield;
        Matrix MEGleadfield;
    };

    class GainInternalPot : public Matrix {
    public:
        using Matrix::operator=;
        GainInternalPot (const SymMatrix& HeadMatInv,const Matrix& SourceMat,const Matrix& Head2IPMat,const Matrix& Source2IPMat):
            Matrix(Source2IPMat+(Head2IPMat*HeadMatInv)*SourceMat)
        { }
    };

    class GainEITInternalPot : public Matrix {
    public:
        using Matrix::operator=;
        GainEITInternalPot (const SymMatrix& HeadMatInv,const Matrix& SourceMat,const Matrix& Head2IPMat):
            Matrix((Head2IPMat*HeadMatInv)*SourceMat)
        { }
    };
}

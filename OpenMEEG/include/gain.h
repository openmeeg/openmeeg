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

#define USE_GMRES 0

#include "matrix.h"
#include "sparse_matrix.h"
#include "symmatrix.h"
#include "geometry.h"
#include "assemble.h"
#include "gmres.h"

namespace OpenMEEG {

    class GainMEG : public Matrix {
    public:
        using Matrix::operator=;
        GainMEG (const Matrix& GainMat): Matrix(GainMat) {}
        GainMEG (const SymMatrix& HeadMatInv,const Matrix& SourceMat, const Matrix& Head2MEGMat, const Matrix& Source2MEGMat) {
            *this = Source2MEGMat+(Head2MEGMat*HeadMatInv)*SourceMat;
        }
        ~GainMEG () {};
    };

    class GainEEG : public Matrix {
    public:
        using Matrix::operator=;
        GainEEG (const Matrix& GainMat): Matrix(GainMat) {}
        GainEEG (const SymMatrix& HeadMatInv,const Matrix& SourceMat, const SparseMatrix& Head2EEGMat) {
            *this = (Head2EEGMat*HeadMatInv)*SourceMat;
        }
        ~GainEEG () {};
    };

    class GainEEGadjoint : public Matrix {
        public:
            using Matrix::operator=;
            GainEEGadjoint (const Geometry& geo,const Matrix& dipoles,const SymMatrix& HeadMat, const SparseMatrix& Head2EEGMat) {
                Matrix LeadField(Head2EEGMat.nlin(),dipoles.nlin());
                int gauss_order = 3;
                // Consider the GMRes solver for problem with dimension > 15,000 (3,000 vertices per interface) else use LAPACK solver
                #if USE_GMRES
                Matrix mtemp(Head2EEGMat.nlin(),HeadMat.nlin());
                Jacobi<SymMatrix> M(HeadMat);    // Jacobi preconditionner
                #pragma omp parallel for
                #ifndef OPENMP_3_0
                for (int i=0;i<static_cast<int>(LeadField.nlin());++i) {
                #else
                for (unsigned i=0;i<LeadField.nlin();++i) {
                #endif
                    Vector vtemp(HeadMat.nlin());
                    GMRes(HeadMat, M, vtemp, Head2EEGMat.getlin(i), 1e3, 1e-7, HeadMat.nlin()); // max number of iteration = 1000, and precision=1e-7 (1e-5 for faster resolution)
                    mtemp.setlin(i, vtemp);
                    #pragma omp critical
                    PROGRESSBAR(i, LeadField.nlin());
                }
                #else
                Matrix mtemp(Head2EEGMat.transpose());
                HeadMat.solveLin(mtemp); // solving the system AX=B with LAPACK
                mtemp=mtemp.transpose();
                #endif
                for ( unsigned i = 0; i < LeadField.ncol(); ++i) {
                    LeadField.setcol(i,mtemp * DipSourceMat(geo, dipoles.submat(i, 1, 0, dipoles.ncol()), gauss_order, true, "").getcol(0)); // TODO ugly
                    PROGRESSBAR(i,LeadField.ncol());
                }
                *this = LeadField;
            }
            ~GainEEGadjoint () {};
    };

    class GainMEGadjoint : public Matrix {
        public:
            using Matrix::operator=;
            GainMEGadjoint (const Geometry& geo, const Matrix& dipoles,
                            const SymMatrix& HeadMat,
                            const Matrix& Head2MEGMat,
                            const Matrix& Source2MEGMat) {
                Matrix LeadField(Head2MEGMat.nlin(),dipoles.nlin());
                int gauss_order = 3;
                // Consider the GMRes solver for problem with dimension > 15,000 (3,000 vertices per interface) else use LAPACK solver
                #if USE_GMRES
                Matrix mtemp(Head2MEGMat.nlin(),HeadMat.nlin()); 
                Jacobi<SymMatrix> M(HeadMat);    // Jacobi preconditionner
                #pragma omp parallel for
                #ifndef OPENMP_3_0
                for (int i=0;i<static_cast<int>(LeadField.nlin());++i) {
                #else
                for (unsigned i=0;i<LeadField.nlin();++i) {
                #endif
                    Vector vtemp(HeadMat.nlin());
                    GMRes(HeadMat, M, vtemp, Head2MEGMat.getlin(i), 1e3, 1e-7, HeadMat.nlin()); // max number of iteration = 1000, and precision=1e-7 (1e-5 for faster resolution)
                    mtemp.setlin(i,vtemp);
                    #pragma omp critical
                    PROGRESSBAR(i,LeadField.nlin());
                }
                #else
                Matrix mtemp(Head2MEGMat.transpose());
                HeadMat.solveLin(mtemp); // solving the system AX=B with LAPACK
                mtemp=mtemp.transpose();
                #endif
                for (unsigned i=0;i<LeadField.ncol();i++) {
                    LeadField.setcol(i, mtemp * DipSourceMat(geo, dipoles.submat(i, 1, 0, dipoles.ncol()), gauss_order, true, "").getcol(0)+Source2MEGMat.getcol(i)); // TODO ugly
                    PROGRESSBAR(i,LeadField.ncol());
                }
                *this = LeadField;
            }
            ~GainMEGadjoint () {};
    };

    class GainEEGMEGadjoint {
        public:
            GainEEGMEGadjoint (const Geometry& geo,const Matrix& dipoles,const SymMatrix& HeadMat, const SparseMatrix& Head2EEGMat, const Matrix& Head2MEGMat, const Matrix& Source2MEGMat) {
                unsigned gauss_order = 3;
                this->EEGleadfield = Matrix(Head2EEGMat.nlin(), dipoles.nlin());
                this->MEGleadfield = Matrix(Head2MEGMat.nlin(), dipoles.nlin());
                Matrix RHS(Head2EEGMat.nlin()+Head2MEGMat.nlin(), HeadMat.nlin());

                for ( unsigned i = 0; i < Head2EEGMat.nlin(); ++i) {
                    RHS.setlin(i, Head2EEGMat.getlin(i));
                }
                for ( unsigned i = 0; i < Head2MEGMat.nlin(); ++i) {
                    RHS.setlin(i + Head2EEGMat.nlin(), Head2MEGMat.getlin(i));
                }
                // Consider the GMRes solver for problem with dimension > 15,000 (3,000 vertices per interface) else use LAPACK solver
                #if USE_GMRES
                Matrix mtemp(RHS.nlin(), HeadMat.nlin()); 
                Jacobi<SymMatrix> M(HeadMat); // Jacobi preconditionner
                #pragma omp parallel for
                #ifndef OPENMP_3_0
                for (int i=0;i<static_cast<int>(RHS.nlin());++i) {
                #else
                for (unsigned i=0;i<RHS.nlin();++i) {
                #endif
                    Vector vtemp(HeadMat.nlin());
                    GMRes(HeadMat, M, vtemp, RHS.getlin(i), 1e3, 1e-7, HeadMat.nlin()); // max number of iteration = 1000, and precision=1e-7 (1e-5 for faster resolution)
                    mtemp.setlin(i, vtemp);
                    #pragma omp critical
                    PROGRESSBAR(i, RHS.nlin());
                }
                #else
                Matrix mtemp(RHS.transpose());
                HeadMat.solveLin(mtemp); // solving the system AX=B with LAPACK
                mtemp = mtemp.transpose();
                #endif
                for ( unsigned i = 0; i < dipoles.nlin(); ++i) {
                    Vector dsm = DipSourceMat(geo,dipoles.submat(i, 1, 0, dipoles.ncol()), gauss_order, true, "").getcol(0); // TODO ugly
                    EEGleadfield.setcol(i, mtemp.submat(0, Head2EEGMat.nlin(), 0, HeadMat.nlin()) * dsm);
                    MEGleadfield.setcol(i, mtemp.submat(Head2EEGMat.nlin(), Head2MEGMat.nlin(), 0, HeadMat.nlin()) * dsm + Source2MEGMat.getcol(i));
                    PROGRESSBAR(i, dipoles.nlin());
                }
            }
            
            void saveEEG( const std::string filename ) const { EEGleadfield.save(filename); }

            void saveMEG( const std::string filename ) const { MEGleadfield.save(filename); }
            
            ~GainEEGMEGadjoint () {};
        private:
            Matrix EEGleadfield;
            Matrix MEGleadfield;
    };

    class GainInternalPot : public Matrix {
    public:
        using Matrix::operator=;
        GainInternalPot (const SymMatrix& HeadMatInv, const Matrix& SourceMat, const Matrix& Head2IPMat, const Matrix& Source2IPMat) {
            *this = Source2IPMat + (Head2IPMat * HeadMatInv) * SourceMat;
        }
        ~GainInternalPot () {};
    };

   class GainEITInternalPot : public Matrix {
    public:
        using Matrix::operator=;
        GainEITInternalPot (const SymMatrix& HeadMatInv, const Matrix& SourceMat, const Matrix& Head2IPMat) {
            *this = (Head2IPMat * HeadMatInv) * SourceMat;
        }
        ~GainEITInternalPot () {};
    };
}

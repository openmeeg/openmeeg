// FILE: $Id$

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

#ifndef TRIVIALBINIO_H
#define TRIVIALBINIO_H

#include "MatrixIO.H"
#include "sparse_matrice_dcl.h"
#include "matrice_dcl.h"
#include "symmatrice_dcl.h"

namespace Maths {

    struct TrivialBinIO: public MathsIOBase {

        void read_header(std::istream& is,MatrixBase& matrix) const {
            unsigned int ui;
            is.read(reinterpret_cast<char*>(&ui),sizeof(unsigned int));
            matrix.nlin() = ui;
            if (matrix.storageType()!=MatrixBase::SYMMETRIC && matrix.dimension()!=MatrixBase::ONE)
                is.read(reinterpret_cast<char*>(&ui),sizeof(unsigned int));
            if (matrix.dimension()==MatrixBase::TWO) {
                matrix.ncol() = ui;
            }
        }

        const std::string& identity() const { return Identity; }
        const Suffixes&    suffixes() const { return suffs;    }

        bool identify(const std::string&)    const { return true; }
        bool known(const MatrixBase& matrix) const {
            return matrix.dimension()==MatrixBase::TWO
                   || (matrix.dimension()==MatrixBase::ONE && matrix.storageType()==MatrixBase::FULL);
        }

        void read_sparse(std::ifstream& is,MatrixBase& matrix) const {
            sparse_matrice* M = dynamic_cast<sparse_matrice*>(&matrix);

            while (! is.eof() ) {
                unsigned int ui;
                size_t i,j;
                double val;
                if (is.read(reinterpret_cast<char*>(&ui),sizeof(unsigned int)) == 0) break;
                i = ui;
                if (is.read(reinterpret_cast<char*>(&ui),sizeof(unsigned int)) == 0) break;
                j = ui;
                if (is.read(reinterpret_cast<char*>(&val),sizeof(double)) == 0) break;

                M->operator()(i,j) = val;
            }
        }

        void read_full(std::ifstream& is,MatrixBase& matrix) const {
            if(matrix.dimension()==MatrixBase::ONE) {
                vecteur* M = dynamic_cast<vecteur*>(&matrix);
                M->alloc_data();
                is.read(reinterpret_cast<char*>(M->data()),M->size()*sizeof(double));
            } else {
                matrice* M = dynamic_cast<matrice*>(&matrix);
                M->alloc_data();
                is.read(reinterpret_cast<char*>(M->data()),M->size()*sizeof(double));
            }
        }

        void read_symmetric(std::ifstream& is,MatrixBase& matrix) const {
            symmatrice* M = dynamic_cast<symmatrice*>(&matrix);
            M->alloc_data();
            is.read(reinterpret_cast<char*>(M->data()),M->size()*sizeof(double));
        }

        void read(std::ifstream& is,MatrixBase& matrix) const {
            read_header(is,matrix);
            switch (matrix.storageType()) {
                case MatrixBase::SPARSE :
                    read_sparse(is,matrix);
                    return;
                case MatrixBase::FULL :
                    read_full(is,matrix);
                    return;
                case MatrixBase::SYMMETRIC :
                    read_symmetric(is,matrix);
                    return;
                default:
                    return;
            }
        }

        void write_sparse(std::ofstream& os, const MatrixBase& matrix) const {
            sparse_matrice* M = dynamic_cast<sparse_matrice*>(const_cast<MatrixBase*>(&matrix));
            sparse_matrice::const_iterator it;
            unsigned int ui;
            for(it = M->begin(); it != M->end(); ++it) {
                size_t i = it->first.first;
                size_t j = it->first.second;
                double val = it->second;
                ui = i; os.write(reinterpret_cast<const char*>(&ui),sizeof(unsigned int));
                ui = j; os.write(reinterpret_cast<const char*>(&ui),sizeof(unsigned int));
                os.write(reinterpret_cast<const char*>(&val),sizeof(double));
            }
        }

        void write_symmetric(std::ofstream& os, const MatrixBase& matrix) const {
            symmatrice* M = dynamic_cast<symmatrice*>(const_cast<MatrixBase*>(&matrix));
            os.write(reinterpret_cast<const char*>(M->data()),M->size()*sizeof(double));
        }

        void write_full(std::ofstream& os, const MatrixBase& matrix) const {
            if(matrix.dimension()==MatrixBase::ONE) {
                vecteur* M = dynamic_cast<vecteur*>(const_cast<MatrixBase*>(&matrix));
                os.write(reinterpret_cast<const char*>(M->data()),M->size()*sizeof(double));
            } else {
                matrice* M = dynamic_cast<matrice*>(const_cast<MatrixBase*>(&matrix));
                os.write(reinterpret_cast<const char*>(M->data()),M->size()*sizeof(double));
            }
        }

        void write_header(std::ofstream& os, const MatrixBase& matrix) const {
            unsigned int ui;
            ui = matrix.nlin();
            os.write(reinterpret_cast<const char*>(&ui),sizeof(unsigned int));
            if (matrix.storageType()!=MatrixBase::SYMMETRIC && matrix.dimension()!=MatrixBase::ONE) {
                ui = matrix.ncol();
                os.write(reinterpret_cast<const char*>(&ui),sizeof(unsigned int));
            }
        }

        void write(std::ofstream& os, const MatrixBase& matrix) const {
            write_header(os,matrix);
            switch (matrix.storageType()) {
                case MatrixBase::SPARSE :
                    write_sparse(os,matrix);
                    return;
                case MatrixBase::FULL :
                    write_full(os,matrix);
                    return;
                case MatrixBase::SYMMETRIC :
                    write_symmetric(os,matrix);
                    return;
                default:
                    return;
            }
        }

    private:

        TrivialBinIO(): MathsIOBase(MathsIO::ios()) { }
        ~TrivialBinIO() {};

        static Suffixes init() {
            Suffixes suffs;
            suffs.push_back("bin");
            return suffs;
        }

        static const TrivialBinIO prototype;
        static const Suffixes     suffs;
        static const std::string  Identity;
    };
}

#endif  //! TRIVIALBINIO_H

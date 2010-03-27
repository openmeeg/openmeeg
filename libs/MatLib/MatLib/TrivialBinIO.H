// FILE: $Id$

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
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

#include "MathsIO.H"
#include "sparse_matrix.h"
#include "matrix.h"
#include "symmatrix.h"

namespace maths {

    struct OPENMEEGMATHS_EXPORT TrivialBinIO: public MathsIOBase {

        TrivialBinIO() { }
        ~TrivialBinIO() {};

        void read_header(std::istream& is,LinOp& linop) const {
            unsigned int ui;
            is.read(reinterpret_cast<char*>(&ui),sizeof(unsigned int));
            linop.nlin() = ui;
            if (linop.storageType()!=LinOp::SYMMETRIC && linop.dimension()!=LinOp::ONE)
                is.read(reinterpret_cast<char*>(&ui),sizeof(unsigned int));
            if (linop.dimension()==LinOp::TWO) {
                linop.ncol() = ui;
            }
        }

        const std::string& identity() const { return Identity; }
        const Suffixes&    suffixes() const { return suffs;    }

        bool identify(const std::string&)    const { return true; }
        bool known(const LinOp& linop) const {
            return linop.dimension()==LinOp::TWO
                   || (linop.dimension()==LinOp::ONE && linop.storageType()==LinOp::FULL);
        }

        void read(std::ifstream& is,LinOp& linop) const {
            read_header(is,linop);
            switch (linop.storageType()) {
                case LinOp::SPARSE :
                    read_sparse(is,linop);
                    return;
                case LinOp::FULL :
                    if(linop.dimension()==LinOp::ONE) {
                        read_internal<Vector>(is,linop);
                    } else {
                        read_internal<Matrix>(is,linop);
                    }
                    return;
                case LinOp::SYMMETRIC :
                    read_internal<SymMatrix>(is,linop);
                    return;
                default:
                    return;
            }
        }

        void write(std::ofstream& os, const LinOp& linop) const {
            write_header(os,linop);
            switch (linop.storageType()) {
                case LinOp::SPARSE :
                    write_sparse(os,linop);
                    return;
                case LinOp::FULL :
                    if(linop.dimension()==LinOp::ONE) {
                        write_internal<Vector>(os,linop);
                    } else {
                        write_internal<Matrix>(os,linop);
                    }
                    return;
                case LinOp::SYMMETRIC :
                    write_internal<SymMatrix>(os,linop);
                    return;
                default:
                    return;
            }
        }

    private:

        void read_sparse(std::ifstream& is,LinOp& linop) const {
            SparseMatrix& m = dynamic_cast<SparseMatrix&>(linop);

            while (!is.eof()) {
                unsigned int ui;
                size_t i,j;
                double val;
                if (is.read(reinterpret_cast<char*>(&ui),sizeof(unsigned int)) == 0) break;
                i = ui;
                if (is.read(reinterpret_cast<char*>(&ui),sizeof(unsigned int)) == 0) break;
                j = ui;
                if (is.read(reinterpret_cast<char*>(&val),sizeof(double)) == 0) break;

                m(i,j) = val;
            }
        }

        void write_sparse(std::ofstream& os, const LinOp& linop) const {
            const SparseMatrix& m = dynamic_cast<const SparseMatrix&>(linop);
            for (SparseMatrix::const_iterator it=m.begin();it!=m.end();++it) {
                const unsigned i = it->first.first;
                const unsigned j = it->first.second;
                const double val = it->second;
                os.write(reinterpret_cast<const char*>(&i),sizeof(unsigned int));
                os.write(reinterpret_cast<const char*>(&j),sizeof(unsigned int));
                os.write(reinterpret_cast<const char*>(&val),sizeof(double));
            }
        }

        template <typename LINOP>
        static void read_internal(std::ifstream& is,LinOp& linop) {
            LINOP& l = dynamic_cast<LINOP&>(linop);
            l.alloc_data();
            is.read(reinterpret_cast<char*>(l.data()),l.size()*sizeof(double));
        }

        template <typename LINOP>
        static void write_internal(std::ofstream& os, const LinOp& linop) {
            const LINOP& l = dynamic_cast<const LINOP&>(linop);
            os.write(reinterpret_cast<const char*>(l.data()),l.size()*sizeof(double));
        }

        void write_header(std::ofstream& os,const LinOp& linop) const {
            const unsigned nlin = linop.nlin();
            os.write(reinterpret_cast<const char*>(&nlin),sizeof(unsigned));
            if (linop.storageType()!=LinOp::SYMMETRIC && linop.dimension()!=LinOp::ONE) {
                const unsigned ncol = linop.ncol();
                os.write(reinterpret_cast<const char*>(&ncol),sizeof(unsigned));
            }
        }

        static Suffixes init() {
            Suffixes suffs;
            suffs.push_back("bin");
            return suffs;
        }

        static const Suffixes     suffs;
        static const std::string  Identity;
    };
}

#endif  //! TRIVIALBINIO_H

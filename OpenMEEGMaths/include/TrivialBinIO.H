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

#include "MathsIO.H"
#include "sparse_matrix.h"
#include "matrix.h"
#include "symmatrix.h"

namespace OpenMEEG {

    namespace maths {

        struct OPENMEEGMATHS_EXPORT TrivialBinIO: public MathsIOBase {

            const std::string& identity() const { return Identity; }
            const Suffixes&    suffixes() const { return suffs;    }

            bool identify(const std::string&)    const { return true; }
            bool known(const LinOp& linop) const {
                return linop.dimension()==2
                       || (linop.dimension()==1 && linop.storageType()==LinOp::FULL);
            }

            /** \brief Read Matrix dimensions for raw binary file without loading the full data
                \sa
            **/

            LinOpInfo info(std::ifstream& is) const {

                //  Get the size of the file.

                std::ifstream::pos_type pos = is.tellg();
                is.seekg(0, std::ios::end);
                const std::streamoff size = static_cast<std::streamoff>(is.tellg())-sizeof(unsigned int);
                is.seekg(pos);

                //  Attempt to find the type of the file from its size.

                LinOpInfo linop;
                unsigned ui;
                is.read(reinterpret_cast<char*>(&ui),sizeof(unsigned));
                linop.nlin() = ui;

                if (size==static_cast<std::streamoff>(ui*sizeof(double))) {
                    linop.storageType() = LinOp::FULL;
                    linop.dimension()   = 1;
                    linop.ncol()        = 1;
                } else {
                    linop.dimension()   = 2;
                    if (size==static_cast<std::streamoff>((long)ui*(ui+1)*sizeof(double)/2)) {
                        linop.storageType() = LinOp::SYMMETRIC;
                        linop.ncol()        = ui;
                    } else {
                        is.read(reinterpret_cast<char*>(&ui),sizeof(unsigned int));
                        linop.ncol()        = ui;
                        const std::streamoff n_elmts = linop.nlin()*linop.ncol();
                        linop.storageType() = (size-sizeof(unsigned int)==n_elmts*sizeof(double)) ?
                                                LinOp::FULL : LinOp::SPARSE;
                    }
                }

                return linop;
            }

            void read(std::ifstream& is,LinOp& linop) const {
                const LinOpInfo& inforead = info(is);

                if (linop.storageType()!=inforead.storageType())
                    throw BadStorageType(name());

                if (linop.dimension()!=inforead.dimension())
                    throw BadStorageType(name());

                linop.nlin() = inforead.nlin();
                linop.ncol() = inforead.ncol();

                switch (linop.storageType()) {
                    case LinOp::SPARSE :
                        read_sparse(is,linop);
                        return;
                    case LinOp::FULL :
                        if (linop.dimension()==1) {
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

                //  Write the header.

                const unsigned nlin = linop.nlin();
                os.write(reinterpret_cast<const char*>(&nlin),sizeof(unsigned));
                if (linop.storageType()!=LinOp::SYMMETRIC && linop.dimension()!=1) {
                    const unsigned ncol = linop.ncol();
                    os.write(reinterpret_cast<const char*>(&ncol),sizeof(unsigned));
                }

                //  Write the data.

                switch (linop.storageType()) {
                    case LinOp::SPARSE :
                        write_sparse(os,linop);
                        return;
                    case LinOp::FULL :
                        if(linop.dimension()==1) {
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

                    if (!is.read(reinterpret_cast<char*>(&ui), sizeof(unsigned int))) break;
                    const size_t i = ui;

                    if (!is.read(reinterpret_cast<char*>(&ui), sizeof(unsigned int))) break;
                    const size_t j = ui;

                    double val;
                    if (!is.read(reinterpret_cast<char*>(&val), sizeof(double))) break;

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

#ifdef NOBUG
                is.read(reinterpret_cast<char*>(l.data()),l.size()*sizeof(double));
#else
                // There seems to be a bug in the C++ library on MAC
                // XXX : Read line by line to avoid overflow on MAC Grrr !!!
                size_t buffer_size = l.ncol()*sizeof(double);
                for(size_t i = 0; i < linop.nlin(); ++i)
                    is.read(reinterpret_cast<char*>(l.data()) + i*buffer_size, buffer_size);
#endif
            }

            template <typename LINOP>
            static void write_internal(std::ofstream& os, const LinOp& linop) {
                const LINOP& l = dynamic_cast<const LINOP&>(linop);
                os.write(reinterpret_cast<const char*>(l.data()),l.size()*sizeof(double));
            }

            TrivialBinIO(): MathsIOBase(20) { }
            ~TrivialBinIO() {};

            static Suffixes init() {
                Suffixes suffixes;
                suffixes.push_back("bin");
                return suffixes;
            }

            static const TrivialBinIO prototype;
            static const Suffixes     suffs;
            static const std::string  Identity;
        };
    }
}

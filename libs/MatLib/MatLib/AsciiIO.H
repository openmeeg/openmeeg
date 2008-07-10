// FILE: $Id: AsciiIO.H 252 2008-06-27 14:43:17Z papadop $

/*
Project Name : OpenMEEG

version           : $Revision: 252 $
last revision     : $Date: 2008-06-27 16:43:17 +0200 (Ven, 27 jui 2008) $
modified by       : $LastChangedBy: papadop $
last modified     : $LastChangedDate: 2008-06-27 16:43:17 +0200 (Ven, 27 jui 2008) $

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

#ifndef ASCIIIO_H
#define ASCIIIO_H

#include <sstream>

#ifdef USE_MATIO
#include <matio.h>
#endif

#include <Exceptions.H>
#include "MathsIO.H"
#include "sparse_matrix.h"
#include "matrix.h"
#include "symmatrix.h"

namespace maths {

    struct AsciiIO: public MathsIOBase {

        const std::string& identity() const { return Identity; }
        const Suffixes&    suffixes() const { return suffs;    }

        bool identify(const std::string& buffer) const {
            double tmp;
            std::stringstream ss(buffer);
            return (ss >> tmp) ? true : false; // test if file starts with a proper float value
        }

        bool known(const LinOp& linop) const { return true; }

        void read_sparse(std::ifstream& is,LinOp& linop) const {

            SparseMatrix& m = dynamic_cast<SparseMatrix&>(linop);

            is.clear();
            is.seekg(0,std::ios::beg);

            std::string line;

            // Start by reading dimensions
            std::getline(is, line, '\n');
            std::stringstream buffer(line);
            buffer >> m.nlin();
            buffer >> m.ncol();

            while (std::getline(is, line, '\n')) {
                std::stringstream buffer(line);
                size_t i,j;
                buffer >> i >> j;
                buffer >> m(i,j);
            }
            // throw ?
        }

        void read_full(std::ifstream& is,LinOp& linop) const {
            is.clear();
            is.seekg(0,std::ios::beg);

            std::string line;
            if (linop.dimension()==LinOp::ONE) {
                Vector& v = dynamic_cast<Vector&>(linop);
                v.alloc_data();

                for (size_t i=0;i<v.nlin();++i) {
                    std::getline(is,line,'\n');
                    std::stringstream buffer(line);
                    buffer >> v(i);
                    if (buffer.fail())
                        throw("Error parsing"); //  FIXME : throw proper exception
                }
            } else {
                Matrix& m = dynamic_cast<Matrix&>(linop);
                m.alloc_data();

                for (size_t i=0;i<m.nlin();++i) {
                    std::getline(is,line,'\n');
                    std::stringstream buffer(line);
                    for (size_t j=0;j<m.ncol();++j)
                        buffer >> m(i,j);
                    if (buffer.fail())
                        throw("Error parsing"); //  FIXME : throw proper exception
                }
            }
        }

        void read_symmetric(std::ifstream& is,LinOp& linop) const {
            SymMatrix& m = dynamic_cast<SymMatrix&>(linop);
            m.alloc_data();

            // Read data

            is.clear();
            is.seekg(0,std::ios::beg);

            std::string line;
            for (size_t i=0;i<m.nlin();++i) {
                std::getline(is,line,'\n');
                std::stringstream buffer(line);
                for (size_t j=i;j<m.ncol();++j)
                    buffer >> m(i,j);
                if (buffer.fail())
                    throw("Error parsing"); //  FIXME : throw proper exception
            }
        }

        void read(std::ifstream& is,LinOp& linop) const {

            if (linop.storageType()==LinOp::SPARSE) {

                read_sparse(is,linop);

            } else {

                // Determine the number of columns on the first line

                is.clear();
                is.seekg(0,std::ios::beg);

                unsigned ncol = 0;
                std::string line;
                std::getline(is, line, '\n');
                std::stringstream buffer(line);
                while(!buffer.fail()) {
                    double d;
                    buffer >> d;
                    ++ncol;
                }
                --ncol;

                if (linop.dimension()==LinOp::ONE) {
                    assert(ncol == 1); // TODO : throw proper exception
                }

                // Determine the number of lines in the file

                is.clear();
                is.seekg(0,std::ios::beg);

                unsigned nlin = 0;
                while (std::getline(is, line, '\n'))
                    if (!line.empty())
                        ++nlin;

                //  Read the data according to the type of the matrix.

                if (linop.storageType()==LinOp::SYMMETRIC) {
                    assert(nlin==ncol); // TODO : throw proper exception
                    linop.ncol() = ncol;
                    read_symmetric(is,linop);
                } else {
                    linop.nlin() = nlin;
                    linop.ncol() = ncol;
                    read_full(is,linop);
                }
            }
        }

        void write_sparse(std::ofstream& os, const LinOp& linop) const {
            const SparseMatrix& spm = dynamic_cast<const SparseMatrix&>(linop);
            SparseMatrix::const_iterator it;
            os << spm.nlin() << " " << spm.ncol() << std::endl;
            for(it = spm.begin(); it != spm.end(); ++it) {
                size_t i = it->first.first;
                size_t j = it->first.second;
                double val = it->second;
                os << i << " " << j << " " << val << std::endl;
            }
        }

        void write_symmetric(std::ofstream& os, const LinOp& linop) const {
            const SymMatrix& m = dynamic_cast<const SymMatrix&>(linop);
            for (size_t i=0;i<m.nlin();++i)
                for (size_t j=i;j<m.ncol();++j)
                    os << m(i,j) << ((j!=m.ncol()-1) ? "\t" : "\n");
        }

        void write_full(std::ofstream& os, const LinOp& linop) const {
            if (linop.dimension()==LinOp::ONE) {
                const Vector& v = dynamic_cast<const Vector&>(linop);
                for (size_t i=0;i<v.nlin();++i)
                    os << v(i) << std::endl;
            } else {
                const Matrix& m = dynamic_cast<const Matrix&>(linop);
                for (size_t i=0;i<m.nlin();++i)
                    for (size_t j=0;j<m.ncol();++j)
                        os << m(i,j) << ((j!=m.ncol()-1) ? "\t" : "\n");
            }
        }

        void write(std::ofstream& os, const LinOp& linop) const {
            switch (linop.storageType()) {
                case LinOp::SPARSE :
                    write_sparse(os,linop);
                    return;
                case LinOp::FULL :
                    write_full(os,linop);
                    return;
                case LinOp::SYMMETRIC :
                    write_symmetric(os,linop);
                    return;
                default:
                    return;
            }
        }

    private:

        AsciiIO(): MathsIOBase(MathsIO::ios()) { }
        ~AsciiIO() {};

        static Suffixes init() {
            Suffixes suffs;
            suffs.push_back("txt");
            return suffs;
        }

        static const AsciiIO     prototype;
        static const Suffixes    suffs;
        static const std::string Identity;
    };
}

#endif  //! ASCIIIO_H

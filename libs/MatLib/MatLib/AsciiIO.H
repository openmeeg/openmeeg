// FILE: $Id: AsciiIO.H 252 2008-06-27 14:43:17Z papadop $

/*
Project Name : OpenMEEG

author            : $Author: papadop $
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

#ifdef USE_MATIO
#include <matio.h>
#endif

#include <sstream>

#include "MatrixIO.H"
#include "sparse_matrice.h"
#include "matrice.h"
#include "symmatrice.h"

namespace Maths {

    struct AsciiIO: public MathsIOBase {

        const std::string& identity() const { return Identity; }
        const Suffixes&    suffixes() const { return suffs;    }

        bool identify(const std::string& buffer) const {
            double tmp;
            std::stringstream ss(buffer);
            if(ss >> tmp) { // test if file starts with a proper float value
                return true;
            } else {
                return false;
            }
        }

        bool known(const MatrixBase& matrix) const {
            return true;
        }

        void read_sparse(std::ifstream& is,MatrixBase& matrix) const {
            size_t i,j;
            double val;
            std::stringstream buffer;
            std::string line;

            sparse_matrice* M = dynamic_cast<sparse_matrice*>(const_cast<MatrixBase*>(&matrix));

            is.clear();
            is.seekg(0,std::ios::beg);

            while (std::getline(is, line, '\n')) {
                buffer.clear();
                buffer.str(line);
                buffer >> i >> j >> val;
                M->operator()(i,j) = val;
            }
        }

        void read_full(std::ifstream& is,MatrixBase& matrix) const {
            std::stringstream buffer;
            std::string line;
            is.clear();
            is.seekg(0,std::ios::beg);

            if(matrix.dimension()==MatrixBase::ONE) {
                vecteur* M = dynamic_cast<vecteur*>(&matrix);
                M->alloc_data();

                for(size_t i=0;i<M->nlin();i++)
                {
                    std::getline(is,line,'\n');
                    buffer.clear();
                    buffer.str(line);
                    buffer>>M->data()[i];
                    if(buffer.fail())
                    {
                        std::cerr << "Error parsing" << std::endl;
                        return;
                    }
                }
            } else {
                matrice* M = dynamic_cast<matrice*>(&matrix);
                M->alloc_data();

                for(size_t i=0;i<M->nlin();i++)
                {
                    std::getline(is,line,'\n');
                    buffer.clear();
                    buffer.str(line);
                    for(size_t j=0;j<M->ncol();j++) {
                        buffer>>M->data()[i+j*M->nlin()];
                    }
                    if(buffer.fail())
                    {
                        std::cerr << "Error parsing" << std::endl;
                        return;
                    }
                }
            }
        }

        void read_symmetric(std::ifstream& is,MatrixBase& matrix) const {
            symmatrice* M = dynamic_cast<symmatrice*>(&matrix);
            M->alloc_data();

            std::stringstream buffer;
            std::string line;

            // Read data
            is.clear();
            is.seekg(0,std::ios::beg);
            for(size_t i=0;i<M->nlin();i++)
            {
                std::getline(is,line,'\n');
                buffer.clear();
                buffer.str(line);
                for(size_t j=i;j<M->ncol();j++) {
                    buffer>>M->data()[i+j*(j+1)/2];
                }
                if(buffer.fail())
                {
                    std::cerr << "Error parsing" << std::endl; //  FIXME : throw proper exception
                    return;
                }
            }
        }

        void read(std::ifstream& is,MatrixBase& matrix) const {

            if(matrix.storageType()==MatrixBase::SPARSE) {
                read_sparse(is,matrix);
                return;
            }

            std::stringstream buffer;
            std::string line;
            std::string s;

            matrix.ncol() = 0;
            matrix.nlin() = 0;
            // Determine the number of columns on the first line
            is.clear();
            is.seekg(0,std::ios::beg);
            std::getline(is, line, '\n');
            buffer.clear();
            buffer.str(line);
            while(!buffer.fail()) {
                buffer >> s;
                matrix.ncol()++;
            }
            matrix.ncol()--;

            if(matrix.dimension()==MatrixBase::ONE) {
                assert(matrix.ncol() == 1); // TODO : throw proper exception
            }

            is.clear();
            is.seekg(0,std::ios::beg);
            // Determine the number of lines in the file
            while (std::getline(is, line, '\n')) {
                if (!line.empty()) {
                    matrix.nlin()++;
                }
            }

            if(matrix.storageType()==MatrixBase::SYMMETRIC) {
                read_symmetric(is,matrix);
                return;
            }

            read_full(is,matrix);
            return;
        }

        void write_sparse(std::ofstream& os, const MatrixBase& matrix) const {
            sparse_matrice* spm = dynamic_cast<sparse_matrice*>(const_cast<MatrixBase*>(&matrix));
            sparse_matrice::const_iterator it;
            for(it = spm->begin(); it != spm->end(); ++it) {
                size_t i = it->first.first;
                size_t j = it->first.second;
                double val = it->second;
                os << i << " " << j << " " << val << std::endl;
            }
        }

        void write_symmetric(std::ofstream& os, const MatrixBase& matrix) const {
            symmatrice* M = dynamic_cast<symmatrice*>(const_cast<MatrixBase*>(&matrix));
            for(size_t i=0;i<M->nlin();i++) {
                for(size_t j=i;j<M->ncol();j++) {
                    os << M->data()[i+j*(j+1)/2];
                    if(j!=M->ncol()-1) os<<"\t"; else os<<"\n";
                }
            }
        }

        void write_full(std::ofstream& os, const MatrixBase& matrix) const {
            if (matrix.dimension()==MatrixBase::ONE) {
                vecteur* M = dynamic_cast<vecteur*>(const_cast<MatrixBase*>(&matrix));
                for(size_t i=0;i<M->nlin();i++) {
                    os << M->data()[i];
                    os<<"\n";
                }
            } else {
                matrice* M = dynamic_cast<matrice*>(const_cast<MatrixBase*>(&matrix));
                for(size_t i=0;i<M->nlin();i++) {
                    for(size_t j=0;j<M->ncol();j++) {
                        os << M->data()[i+j*M->nlin()];
                        if(j!=M->ncol()-1) os<<"\t"; else os<<"\n";
                    }
                }
            }
        }

        void write(std::ofstream& os, const MatrixBase& matrix) const {
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

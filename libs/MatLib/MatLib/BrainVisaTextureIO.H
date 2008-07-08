// FILE: $Id: BrainVisaTextureIO.H 252 2008-06-27 14:43:17Z papadop $

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

#ifndef TRIVIALBINIO_H
#define TRIVIALBINIO_H

#include <sstream>

#include "MathsIO.H"
#include "matrix.h"

namespace maths {

    struct BrainVisaTextureIO: public MathsIOBase {

        const std::string& identity() const { return Identity; }
        const Suffixes&    suffixes() const { return suffs;    }

        bool identify(const std::string& buffer)    const {
            if (buffer.size()<MagicTag.size())
                return false;
            return strncmp(buffer.c_str(),MagicTag.c_str(),MagicTag.size()) == 0;
        }

        bool known(const LinOp& linop) const {
            return linop.storageType()==LinOp::FULL;
        }

        void read(std::ifstream& is,LinOp& linop) const {
            Matrix* m = dynamic_cast<Matrix*>(&linop);

            std::string line;
            std::getline(is, line, '\n');
            std::istringstream buffer(line);
            std::string s;
            buffer >> s;
            assert(s == "ascii");

            std::getline(is, line, '\n'); // skip keyword FLOAT

            std::getline(is, line, '\n');
            buffer.clear();
            buffer.str(line);
            buffer >> m->ncol();

            if(linop.dimension()==LinOp::ONE) {
                assert(m->ncol() == 1); // TODO : throw proper exception
            }

            size_t j = 0; // column index
            do {
                std::getline(is, line, '\n'); // skip column index
                std::getline(is, line, '\n');
                buffer.clear();
                buffer.str(line);
                buffer >> m->nlin();
                if(j==0) {
                    m->alloc_data();
                }
                for(size_t i = 0; i < m->nlin(); ++i) {
                    buffer >> m->data()[i+j*m->nlin()];
                }
                j++;
            } while(j<m->ncol());
        }

        void write(std::ofstream& os, const LinOp& linop) const {
            Matrix* m = dynamic_cast<Matrix*>(const_cast<LinOp*>(&linop));

            os << "ascii" << std::endl << "FLOAT" << std::endl << (unsigned int)m->ncol() << std::endl;
            for(size_t j=0; j<m->ncol(); j++) {
                os << (unsigned int) j << std::endl;
                os << (unsigned int) m->nlin();
                for(size_t i=0; i<m->nlin(); i++) {
                    os << " " << m->data()[i+j*m->nlin()];
                }
                os << std::endl;
            }
        }

    private:

        BrainVisaTextureIO(): MathsIOBase(MathsIO::ios()) { }
        ~BrainVisaTextureIO() {};

        static Suffixes init() {
            Suffixes suffs;
            suffs.push_back("tex");
            return suffs;
        }

        static const BrainVisaTextureIO prototype;
        static const std::string        MagicTag;
        static const Suffixes           suffs;
        static const std::string        Identity;
    };
}

#endif  //! TRIVIALBINIO_H

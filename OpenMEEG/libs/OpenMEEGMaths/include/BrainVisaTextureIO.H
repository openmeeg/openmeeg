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

#include <sstream>

#include <Exceptions.H>
#include <IOUtils.H>
#include "MathsIO.H"
#include "matrix.h"

namespace OpenMEEG {
    namespace maths {

        struct OPENMEEGMATHS_EXPORT BrainVisaTextureIO: public MathsIOBase {

            const std::string& identity() const { return Identity; }
            const Suffixes&    suffixes() const { return suffs;    }

            bool identify(const std::string& buffer)    const {
                if (buffer.size()<MagicTag.size())
                    return false;
                return strncmp(buffer.c_str(),MagicTag.c_str(),MagicTag.size()) == 0;
            }

            bool known(const LinOp& linop) const {
                return linop.storageType()==LinOp::FULL && linop.dimension()==2;
            }

            LinOpInfo info(std::ifstream& is) const {
                size_t nlin,ncol;
                read_header(is,nlin,ncol);
                is.seekg(0,std::ios_base::beg);
                return LinOpInfo(nlin,ncol,LinOpInfo::FULL,2);
            }

            void read(std::ifstream& is,LinOp& linop) const {
                Matrix& m = dynamic_cast<Matrix&>(linop);
                read_header(is,m.nlin(),m.ncol());
                if (linop.dimension()==1 && m.ncol()!=1)
                    throw maths::BadVector(m.ncol());
                m.alloc_data();

                for (size_t j=0;j<m.ncol();++j) {
                    unsigned ind;
                    is >> io_utils::skip_line;   // skip column index
                    is >> ind;
                    if(ind==0)
                        is >> ind; // Hack for certain tex files

                    for (size_t i=0;i<m.nlin();++i)
                        is >> m(i,j);

                    is >> io_utils::skip_line;
                }
            }

            void write(std::ofstream& os,const LinOp& linop) const {
                const Matrix& m = dynamic_cast<const Matrix&>(linop);

                os << "ascii" << std::endl << "FLOAT" << std::endl << static_cast<unsigned>(m.ncol()) << std::endl;
                for (size_t j=0; j<m.ncol(); ++j) {
                    os << static_cast<unsigned>(j) << std::endl << static_cast<unsigned>(m.nlin()) << std::endl;
                    for (size_t i=0; i<m.nlin(); ++i)
                        os << " " << m(i,j);
                    os << std::endl;
                }
            }

        private:

            BrainVisaTextureIO(): MathsIOBase(8) { }
            ~BrainVisaTextureIO() {};

            static void read_header(std::ifstream& is,size_t& nlin,size_t& ncol) {
                is >> io_utils::match(MagicTag) >> io_utils::skip_line;
                is >> io_utils::skip_line;         // skip keyword FLOAT
                is >> ncol >> io_utils::skip_line;
                const std::ifstream::pos_type pos = is.tellg();
                is >> io_utils::skip_line;         // skip column index
                is >> nlin; 
                if (nlin==0)
                    is >> nlin;
                is.seekg(pos);
            }

            static Suffixes init() {
                Suffixes suffixes;
                suffixes.push_back("tex");
                return suffixes;
            }

            static const BrainVisaTextureIO prototype;
            static const std::string        MagicTag;
            static const Suffixes           suffs;
            static const std::string        Identity;
        };
    }
}

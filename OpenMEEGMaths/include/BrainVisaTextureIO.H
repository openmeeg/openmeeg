// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

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
                Dimension nlin,ncol;
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

                for (Index j=0; j<m.ncol(); ++j) {
                    unsigned ind;
                    is >> io_utils::skip_line;   // skip column index
                    is >> ind;
                    if(ind==0)
                        is >> ind; // Hack for certain tex files

                    for (Index i=0; i<m.nlin(); ++i)
                        is >> m(i,j);

                    is >> io_utils::skip_line;
                }
            }

            void write(std::ofstream& os,const LinOp& linop) const {
                const Matrix& m = dynamic_cast<const Matrix&>(linop);

                os << "ascii" << std::endl << "FLOAT" << std::endl << static_cast<unsigned>(m.ncol()) << std::endl;
                for (Index j=0; j<m.ncol(); ++j) {
                    os << static_cast<unsigned>(j) << std::endl << static_cast<unsigned>(m.nlin()) << std::endl;
                    for (Index i=0; i<m.nlin(); ++i)
                        os << " " << m(i,j);
                    os << std::endl;
                }
            }

        private:

            BrainVisaTextureIO(): MathsIOBase(8) { }
            ~BrainVisaTextureIO() {};

            static void read_header(std::ifstream& is,Dimension& nlin,Dimension& ncol) {
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

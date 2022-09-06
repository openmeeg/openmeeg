// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <sstream>

#include <OMMathExceptions.H>
#include <MathsIO.H>
#include <sparse_matrix.h>
#include <matrix.h>
#include <symmatrix.h>

namespace OpenMEEG {
    namespace maths {

        struct OPENMEEGMATHS_EXPORT AsciiIO: public MathsIOBase {

            const std::string& identity() const { return Identity; }
            const Suffixes&    suffixes() const { return suffs;    }

            bool identify(const std::string& buffer) const {
                double tmp;
                std::stringstream ss(buffer);
                return (ss >> tmp) ? true : false; // test if file starts with a proper float value
            }

            bool known(const LinOp&) const { return true; }

            LinOpInfo info(std::ifstream& is) const {
                is.clear();
                is.seekg(0,std::ios::beg);

                std::string line;
                std::getline(is,line,'\n');

                LinOpInfo linop;
                set_type(is,length(line),linop);

                if (linop.storageType()==LinOp::SPARSE) {
                    is.clear();
                    is.seekg(0,std::ios::beg);

                    std::stringstream buffer(line);
                    buffer >> linop.nlin() >> linop.ncol();
                    return linop;
                }

                // Determine the number of lines in the file

                is.clear();
                is.seekg(0,std::ios::beg);

                unsigned nlin = 0;
                while (std::getline(is,line,'\n'))
                    if (!line.empty())
                        ++nlin;
                linop.nlin() = nlin;

                if (linop.storageType()==LinOp::SYMMETRIC && linop.nlin()!=linop.ncol())
                    throw BadSymmMatrix(linop.nlin(),linop.ncol());

                return linop;
            }

            void read(std::ifstream& is,LinOp& linop) const {

                const LinOpInfo& inforead = info(is);

                if (linop.storageType()!=inforead.storageType())
                    throw BadStorageType("AsciiIO reading of " + name());

                if (linop.dimension()!=inforead.dimension())
                    throw BadVector(inforead.ncol());

                linop.nlin() = inforead.nlin();
                linop.ncol() = inforead.ncol();

                //  Read the data according to the type of the matrix.

                if (linop.storageType()==LinOp::SPARSE) {
                    read_sparse(is,linop);
                } else if (linop.storageType()==LinOp::SYMMETRIC) {
                    read_symmetric(is,linop);
                } else {
                    read_full(is,linop);
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

            AsciiIO(): MathsIOBase(9) {}
            ~AsciiIO() {};

            static unsigned length(const std::string& line) {
                std::stringstream buffer(line);
                unsigned ncol = 0;
                while(!buffer.fail()) {
                    double d;
                    buffer >> d;
                    ++ncol;
                }
                return ncol-1;
            }

            static void set_type(std::ifstream& is,const unsigned len,LinOpInfo& linop) {
                linop.ncol() = len;

                std::string line;
                std::getline(is,line,'\n');

                if (is.eof()) {
                    linop.storageType() = LinOp::FULL;
                    linop.dimension() = 2;
                    linop.nlin()      = 1;
                    return;
                }

                const unsigned len1 = length(line);

                if (len1==len)
                    linop.storageType() = LinOp::FULL;
                else if (len1==len-1)
                    linop.storageType() = LinOp::SYMMETRIC;
                else if (len ==2 && len1==3)
                    linop.storageType() = LinOp::SPARSE;
                else {
                    throw UnexpectedException("Should never happen",__FILE__,__LINE__);
                }

                linop.dimension() = (linop.storageType()==LinOp::FULL && len==1) ? 1 : 2;
            }

            void read_sparse(std::ifstream& is,LinOp& linop) const {

                is.clear();
                is.seekg(0,std::ios::beg);

                std::string line;

                // Start by reading dimensions
                std::getline(is, line, '\n');
                std::stringstream buffer(line);

                SparseMatrix& m = dynamic_cast<SparseMatrix&>(linop);

                buffer >> m.nlin();
                buffer >> m.ncol();

                while (std::getline(is, line, '\n')) {
                    std::stringstream buff(line);
                    size_t i,j;
                    buff >> i >> j;
                    buff >> m(i,j);
                }

                //  Check that there is nothing more in the file (otherwise the file is corrupted somehow).

                is >> std::ws;
                if (!is.eof())
                    throw BadData(identity()+" sparse matrix");
            }

            void read_full(std::ifstream& is,LinOp& linop) const {
                is.clear();
                is.seekg(0,std::ios::beg);

                std::string line;
                if (linop.dimension()==1) {
                    Vector& v = dynamic_cast<Vector&>(linop);
                    v.alloc_data();

                    for (size_t i=0;i<v.nlin();++i) {
                        std::getline(is,line,'\n');
                        std::stringstream buffer(line);
                        buffer >> v(i);
                        if (buffer.fail())
                            throw BadData(identity()+" vector");
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
                            throw BadData(identity()+" matrix");
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
                        throw BadData(identity()+" symmetric matrix");
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
                if (linop.dimension()==1) {
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

            static Suffixes init() {
                Suffixes suffixes;
                suffixes.push_back("txt");
                return suffixes;
            }

            static const AsciiIO     prototype;
            static const Suffixes    suffs;
            static const std::string Identity;
        };
    }
}

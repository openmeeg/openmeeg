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

#include <OpenMEEGMathsConfig.h>
#include <iostream>
#include <map>
#include <algorithm>

#include <linop.h>
#include <range.h>
#include <matrix.h>
#include <Exceptions.H>

namespace OpenMEEG::maths {

    /// \brief  Block symmetric matrix class
    /// Block symmetric matrix class

    class OPENMEEGMATHS_EXPORT SymmetricBlockMatrix: public LinOp {

        typedef std::pair<unsigned,unsigned> Index;
        typedef std::map<Index,Matrix>       Blocks;

    public:

        SymmetricBlockMatrix(): LinOp(0,0,BLOCK,2) { }
        SymmetricBlockMatrix(const size_t N): LinOp(N,N,BLOCK_SYMMETRIC,2) { }

        size_t size() const override {
            unsigned sz = 0;
            for (const auto& block : blocks)
                sz += block.second.size();
            return sz;
        };

        void info() const override {
            if ((nlin()==0) && (ncol()==0)) {
                std::cout << "Empty matrix" << std::endl;
                return;
            }

            std::cout << "Symmetric block matrix" << std::endl;
            std::cout << "Dimensions: " << nlin() << " x " << ncol() << std::endl;
            std::cout << "Number of blocks: " << blocks.size() << std::endl;
            std::cout << "Number of coefficients: " << size() << std::endl;
        }

        Matrix& block(const unsigned i,const unsigned j) {
            bool transposed;
            const Index& ind = find_block_indices(i,j,transposed);
            if (transposed)
                throw 1;
            return blocks[ind];
        }

        const Matrix& block(const unsigned i,const unsigned j) const {
            bool transposed;
            const Index& ind = find_block_indices(i,j,transposed);
            if (transposed)
                throw 1;
            return blocks.at(ind);
        }

        void add_block(const Range& ir,const Range& jr) {
            bool transposed;
            const Index& ind  = create_block_indices(ir,jr,transposed);
            const Index  size = (transposed) ? Index({jr.length(),ir.length()}) : Index({ir.length(),jr.length()});
            blocks[ind] = Matrix(size.first,size.second);
        }

        void add_blocks(const Ranges& ir,const Ranges& jr) {
            for (unsigned i=0; i<ir.size(); ++i)
                for (unsigned j=i; j<jr.size(); ++j)
                    add_block(ir[i],jr[j]);
        }

        void set_blocks(const Ranges& r) {
            blocks.clear();
            ranges.clear();
            for (unsigned i=0; i<ranges.size(); ++i)
                for (unsigned j=i; j<ranges.size(); ++j)
                    add_block(ranges[i],ranges[j]);
        }

        double& operator()(const size_t i,const size_t j) {
            bool transposed;
            const Index& ind = find_block_indices(i,j,transposed);
            const size_t inblockindex_i = ((!transposed) ? i : j)-ranges[ind.first].start();
            const size_t inblockindex_j = ((!transposed) ? j : i)-ranges[ind.second].start();
            return blocks[ind](inblockindex_i,inblockindex_j);
        }

        double operator()(const size_t i,const size_t j) const {
            bool transposed;
            const Index& ind = find_block_indices(i,j,transposed);
            const size_t inblockindex_i = ((!transposed) ? i : j)-ranges[ind.first].start();
            const size_t inblockindex_j = ((!transposed) ? j : i)-ranges[ind.second].start();
            return blocks.at(ind)(inblockindex_i,inblockindex_j);
        }

    private:

        unsigned find_block_index(const Range& r) const {
            const unsigned ind = ranges.find_index(r);
            if (ind!=-1)
                throw 1;
            return ind;
        }

        unsigned find_block_index(const unsigned i) const {
            const unsigned ind = ranges.find_index(i);
            if (ind!=-1)
                throw 2;
            return ind;
        }

        unsigned create_block_index(const Range& r) {
            const unsigned ind = ranges.find_index(r);
            if (ind!=-1)
                return ind;
            ranges.push_back(r);
            return ranges.size()-1;
        }

        Index create_block_indices(const Range& ir,const Range& jr,bool& transposed) {
            transposed = ir.start()>jr.start();
            const unsigned iind = create_block_index(ir);
            const unsigned jind = create_block_index(jr);
            return (transposed) ? Index({jind,iind}) : Index({iind,jind});
        }

        Index find_block_indices(const Range& ir,const Range& jr,bool& transposed) const {
            transposed = ir.start()>jr.start();
            const unsigned iind = find_block_index(ir);
            const unsigned jind = find_block_index(jr);
            return (transposed) ? Index({jind,iind}) : Index({iind,jind});
        }

        Index find_block_indices(const unsigned i,const unsigned j,bool& transposed) const {
            transposed = i>j;
            const unsigned iind = find_block_index(i);
            const unsigned jind = find_block_index(j);
            return (transposed) ? Index({jind,iind}) : Index({iind,jind});
        }

        Ranges ranges;
        Blocks blocks;
    };
}

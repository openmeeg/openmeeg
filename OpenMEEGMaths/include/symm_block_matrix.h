// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <iostream>
#include <map>
#include <algorithm>

#include <linop.h>
#include <range.h>
#include <ranges.h>
#include <matrix.h>
#include <OMMathExceptions.H>

namespace OpenMEEG::maths {

    /// \brief Block symmetric matrix class
    /// Block symmetric matrix class

    class SymmetricBlockMatrix: public LinOp {

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

        void set_blocks(const Ranges& r) {
            blocks.clear();
            ranges.clear();
            for (unsigned i=0; i<r.size(); ++i)
                for (unsigned j=i; j<r.size(); ++j)
                    add_block(r[i],r[j]);
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

        unsigned find_block_index(const Range& r)   const { return ranges.find_index(r); }
        unsigned find_block_index(const unsigned i) const { return ranges.find_index(i); }

        unsigned create_block_index(const Range& r) try {
            const unsigned ind = ranges.find_index(r);
            return ind;
        } catch(...) {
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

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

namespace OpenMEEG::maths {

    /// \brief  Block matrix class
    /// Block matrix class

    class BlockMatrix: public LinOp {

        typedef std::pair<unsigned,unsigned> Index;
        typedef std::map<Index,Matrix>       Blocks;

    public:

        BlockMatrix(): LinOp(0,0,BLOCK,2) { }
        BlockMatrix(const size_t M,const size_t N): LinOp(M,N,BLOCK,2) { }

        size_t size() const override {
            unsigned sz = 0;
            for (const auto& block : all_blocks)
                sz += block.second.size();
            return sz;
        };

        void info() const override {
            if ((nlin()==0) && (ncol()==0)) {
                std::cout << "Empty matrix" << std::endl;
                return;
            }

            std::cout << "Block matrix" << std::endl;
            std::cout << "Dimensions: " << nlin() << " x " << ncol() << std::endl;
            std::cout << "Number of blocks: " << all_blocks.size() << std::endl;
            std::cout << "Number of coefficients: " << size() << std::endl;
        }

              Matrix& block(const unsigned i,const unsigned j)       { return all_blocks[{i,j}];    }
        const Matrix& block(const unsigned i,const unsigned j) const { return all_blocks.at({i,j}); }

        const Blocks& blocks() const { return all_blocks; }

        void add_block(const Range& ir,const Range& jr) {
            const unsigned iind = row_ranges.add(ir);
            const unsigned jind = col_ranges.add(jr);
            const Index inds = { iind, jind };
            all_blocks[inds] = Matrix(ir.length(),jr.length());
        }

        void set_blocks(const Ranges& rows,const Ranges& cols) {
            row_ranges = rows;
            col_ranges = cols;
            for (const auto& ir : row_ranges)
                for (const auto& jr : col_ranges)
                    add_block(ir,jr);
        }

        double& operator()(const size_t i,const size_t j) {
            const Index& ind = find_block_indices(i,j);
            const size_t inblockindex_i = i-row_ranges[ind.first].start();
            const size_t inblockindex_j = j-col_ranges[ind.second].start();
            return all_blocks[ind](inblockindex_i,inblockindex_j);
        }

        double  operator()(const size_t i,const size_t j) const {
            const Index& ind = find_block_indices(i,j);
            const size_t inblockindex_i = i-row_ranges[ind.first].start();
            const size_t inblockindex_j = j-col_ranges[ind.second].start();
            return all_blocks.at(ind)(inblockindex_i,inblockindex_j);
        }

    private:

        Index find_block_indices(const Range& ir,const Range& jr) const {
            const unsigned iind = row_ranges.find_index(ir);
            const unsigned jind = col_ranges.find_index(jr);
            return {iind,jind};
        }

        Index find_block_indices(const unsigned i,const unsigned j) const {
            const unsigned iind = row_ranges.find_index(i);
            const unsigned jind = col_ranges.find_index(j);
            return { iind, jind };
        }

        Ranges row_ranges;
        Ranges col_ranges;
        Blocks all_blocks;
    };

    inline std::ostream& operator<<(std::ostream& os,const BlockMatrix& bm) {
        for (const auto& block : bm.blocks())
            os << "Block " << block.first.first << ',' << block.first.second << std::endl;
            #if 0
            os << "Block " << block.first.first << ',' << block.first.second << std::endl
               << block.second << std::endl;
            #endif
        return os;
    }
}

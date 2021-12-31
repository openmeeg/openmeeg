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

    /// \brief  Block matrix class
    /// Block matrix class

    class OPENMEEGMATHS_EXPORT BlockMatrix: public LinOp {

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
            Index inds = find_block_indices(ir,jr);
            if (inds.first==-1) {
                row_ranges.push_back(ir);
                inds.first = row_ranges.size()-1;
            }
            if (inds.second==-1) {
                col_ranges.push_back(jr);
                inds.second = col_ranges.size()-1;
            }
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

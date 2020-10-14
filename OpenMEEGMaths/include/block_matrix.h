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
#include <vector>
#include <map>

#include <linop.h>
#include <range.h>
#include <Exceptions.H>

namespace OpenMEEG::maths {

    /// \brief  Block matrix class
    /// Block matrix class

    class OPENMEEGMATHS_EXPORT BlockMatrix: public LinOp {

        typedef std::vector<Range>           Ranges;
        typedef std::pair<unsigned,unsigned> Index;
        typedef std::map<Index,Matrix>       Blocks;

    public:

        BlockMatrix(): LinOp(0,0,BLOCK,2) { }
        BlockMatrix(const size_t M,const size_t N): LinOp(M,N,BLOCK,2) { }

              Matrix& block(const unsigned i,const unsigned j)       { return blocks[{i,j}];    }
        const Matrix& block(const unsigned i,const unsigned j) const { return blocks.at({i,j}); }

        void add_block(const Range& ir,const Range& jr) {
            Index inds = find_block_indices(ir,jr);
            if (inds.first==-1) {
                iranges.push_back(ir);
                inds.first = iranges.size()-1;
            }
            if (inds.second==-1) {
                jranges.push_back(jr);
                inds.second = jranges.size()-1;
            }
            blocks[inds] = Matrix(ir.length(),jr.length());
        }

        void set_blocks(const Ranges& rs) {
            iranges = rs;
            jranges = rs;
            for (const auto& ir : iranges)
                for (const auto& jr : jranges)
                    add_block(ir,jr);
        }

        double& operator()(const size_t i,const size_t j) {
            const Index& ind = find_block_indices(i,j);
            const size_t inblockindex_i = i-iranges[ind.first].start();
            const size_t inblockindex_j = j-jranges[ind.second].start();
            return blocks[ind](inblockindex_i,inblockindex_j);
        }

        double  operator()(const size_t i,const size_t j) const {
            const Index& ind = find_block_indices(i,j);
            const size_t inblockindex_i = i-iranges[ind.first].start();
            const size_t inblockindex_j = j-jranges[ind.second].start();
            return blocks.at(ind)(inblockindex_i,inblockindex_j);
        }

    private:

        Index find_block_indices(const Range& ir,const Range& jr) const {
            const unsigned iind = find_overlapping_range(ir,iranges);
            if (iind!=-1 && ir!=iranges[iind])
                throw OverlappingRanges(ir,iranges[iind]);
            const unsigned jind = find_overlapping_range(jr,jranges);
            if (jind!=-1 && jr!=iranges[jind])
                throw OverlappingRanges(jr,iranges[jind]);
            return {iind,jind};
        }

        static unsigned find_overlapping_range(const Range& r,const Ranges& ranges) {
            for (unsigned i=0;i<ranges.size();++i)
                if (ranges[i].intersect(r))
                    return i;
            return -1;
        }

        static unsigned find_block_index(const size_t ind,const Ranges& ranges) {
            for (unsigned i=0;i<ranges.size();++i)
                if (ranges[i].contains(ind))
                    return i;
            throw NonExistingBlock(ind);
        }
        
        Index find_block_indices(const size_t i,const unsigned j) const {
            const unsigned iind = find_block_index(i,iranges);
            const unsigned jind = find_block_index(j,jranges);
            return { iind, jind };
        }

        Ranges iranges;
        Ranges jranges;
        Blocks blocks;
    };
}

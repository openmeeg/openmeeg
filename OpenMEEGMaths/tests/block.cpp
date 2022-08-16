// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <iostream>
#include <range.h>
#include <ranges.h>
#include <block_matrix.h>
#include <symm_block_matrix.h>

int main() {

    using namespace OpenMEEG::maths;

    Ranges row_ranges;
    Ranges col_ranges;

    row_ranges.push_back(Range(1,4));
    row_ranges.push_back(Range(5,6));
    row_ranges.push_back(Range(7,9));

    col_ranges.push_back(Range(1,2));
    col_ranges.push_back(Range(5,9));

    BlockMatrix bm(10,10);
    bm.set_blocks(row_ranges,col_ranges);
    bm.info();

    std::cout << bm << std::endl;

    // section SymMatrix

    SymmetricBlockMatrix sbm(10);
    sbm.set_blocks(row_ranges);
    sbm.info();
}

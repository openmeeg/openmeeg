// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <map>
#include <functional>
#include <matrix.h>
#include <symmatrix.h>

namespace OpenMEEG {

    namespace Details {

        // We need the following function to be able to index a matrix in its transposed form.
        // This is needed for when we fill Dstar through the D operator.

        template <bool Transposed> struct Indexer;

        template <>
        struct Indexer<false> {
            using IndexFunc = std::function<std::size_t(std::size_t)>;

            template <typename Matrix>
            static double&
            value(Matrix& M,const IndexFunc& ind1,const IndexFunc& ind2,const unsigned i,const unsigned j) {
                return M(ind1(i),ind2(j));
            }
            template <typename Matrix>
            static double
            value(const Matrix& M,const IndexFunc& ind1,const IndexFunc& ind2,const unsigned i,const unsigned j) {
                return M(ind1(i),ind2(j));
            }
        };

        template <>
        struct Indexer<true> {
            using IndexFunc = std::function<std::size_t(std::size_t)>;

            template <typename Matrix>
            static double&
            value(Matrix& M,const IndexFunc& ind1,const IndexFunc& ind2,const unsigned i,const unsigned j) {
                return M(ind1(j),ind2(i));
            }

            template <typename Matrix>
            static double
            value(const Matrix& M,const IndexFunc& ind1,const IndexFunc& ind2,const unsigned i,const unsigned j) {
                return M(ind1(j),ind2(i));
            }
        };
    }

    // This class in an adapter that allows to use a matrix after a conversion of its indices.
    // It allows to represent a block in a bigger matrix, while keeping the indexing of that bigger matrix.

    template <typename Matrix,bool Transposed=false>
    class Block {

        using IndexFunc = std::function<std::size_t(std::size_t)>;

    public:

        Block(Matrix& M,const IndexFunc& f1,const IndexFunc& f2): matrix(M),ind1(f1),ind2(f2) { }

        double  operator()(const size_t i,const size_t j) const { return Details::Indexer<Transposed>::value(matrix,ind1,ind2,i,j); }
        double& operator()(const size_t i,const size_t j)       { return Details::Indexer<Transposed>::value(matrix,ind1,ind2,i,j); }

        Block<Matrix,!Transposed> transpose() const { return Block<Matrix,!Transposed>(matrix,ind1,ind2); }

        bool empty() const { return matrix.empty(); }

    protected:
        
        Matrix&   matrix;
        IndexFunc ind1;
        IndexFunc ind2;
    };
}


// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <OMassert.H>
#include <map>
#include <utility>

#include <linop.h>
#include <vector.h>
#include <matrix.h>

// #ifdef WIN32
// #pragma warning( disable : 4251)    //MSVC warning C4251 : DLL exports of STL templates
// #endif

#ifdef WIN32
    template class OPENMEEGMATHS_EXPORT std::map< std::pair< size_t, size_t >, double >;
#endif

namespace OpenMEEG {

    class SymMatrix;

    class OPENMEEGMATHS_EXPORT SparseMatrix : public LinOp {

    public:

        typedef std::map< std::pair< size_t, size_t >, double > Tank;
        typedef std::map< std::pair< size_t, size_t >, double >::const_iterator const_iterator;
        typedef std::map< std::pair< size_t, size_t >, double >::iterator iterator;

        SparseMatrix() : LinOp(0,0,SPARSE,2) {};
        SparseMatrix(const char* fname) : LinOp(0,0,SPARSE,2) { this->load(fname); }
        SparseMatrix(size_t N,size_t M) : LinOp(N,M,SPARSE,2) {};
        ~SparseMatrix() {};

        inline double operator()( size_t i, size_t j ) const {
            om_assert(i < nlin());
            om_assert(j < ncol());
            const_iterator it = m_tank.find(std::make_pair(i, j));
            if (it != m_tank.end()) return it->second;
            else return 0.0;
        }

        inline double& operator()( size_t i, size_t j ) {
            om_assert(i < nlin());
            om_assert(j < ncol());
            return m_tank[ std::make_pair( i, j ) ];
        }

        size_t size() const {
            return m_tank.size();
        }

        const_iterator begin() const { return m_tank.begin(); }
        const_iterator end()   const { return m_tank.end(); }

        SparseMatrix transpose() const;

        const Tank& tank() const {return m_tank;}

        void set( double t);
        Vector getlin(size_t i) const;
        void setlin(Vector v, size_t i);

        void save(const char *filename) const;
        void load(const char *filename);

        void save(const std::string& s) const { save(s.c_str()); }
        void load(const std::string& s)       { load(s.c_str()); }

        void info() const;
        double frobenius_norm() const;

        Vector       operator*( const Vector &x ) const;
        Matrix       operator*( const SymMatrix &m ) const;
        Matrix       operator*( const Matrix &m ) const;
        SparseMatrix operator*( const SparseMatrix &m ) const;
        SparseMatrix operator+( const SparseMatrix &m ) const;

    private:

        Tank m_tank;
    };

    inline Vector SparseMatrix::getlin(size_t i) const {
        om_assert(i<nlin());
        Vector v(ncol());
        for (size_t j=0;j<ncol();j++){
            const_iterator it = m_tank.find(std::make_pair(i, j));
            if (it != m_tank.end()) v(j)=it->second;
            else v(j)=0.0;
        }
        return v;
    }

    inline void SparseMatrix::setlin(Vector v, size_t i) {
        om_assert(i<nlin());
        for (size_t j=0;j<v.nlin();j++){
            (*this)(i,j) = v(j);
        }
    }
}

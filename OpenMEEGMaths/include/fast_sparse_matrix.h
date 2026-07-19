// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include "OpenMEEGMathsConfig.h"
#include "vector.h"
#include "sparse_matrix.h"

namespace OpenMEEG {

    class OPENMEEGMATHS_EXPORT FastSparseMatrix {
    public:

        friend std::ostream& operator<<(std::ostream& f,const FastSparseMatrix &M);

    public:

        FastSparseMatrix(const size_t n,const size_t p,const size_t sp) { alloc(n,p,sp); }

        FastSparseMatrix(): FastSparseMatrix(1,1,1) { }

        FastSparseMatrix(const SparseMatrix& M): FastSparseMatrix(M.nlin(),M.ncol(),M.size()) {

            unsigned cnt = 0;
            size_t current_line = static_cast<size_t>(-1);
            for (SparseMatrix::const_iterator it=M.begin(); it!=M.end(); ++it,++cnt) {
                const size_t i = it->first.first;
                const size_t j = it->first.second;
                const double val = it->second;
                tank[cnt] = val;
                js[cnt]   = j;
                if (i!=current_line) {
                    for (size_t k=current_line+1; k<=i; ++k)
                        rowindex[k] = cnt;
                    current_line = i;
                }
            }

            for (size_t k=current_line+1; k<=M.nlin(); ++k)
                rowindex[k] = M.size();
        }

        FastSparseMatrix(const FastSparseMatrix& M) {
            alloc(M.m_nlin,M.m_ncol,M.rowindex[M.m_nlin]);
            copy(M);
        }

        void operator=(const FastSparseMatrix& M) {
            destroy();
            alloc(M.m_nlin,M.m_ncol,M.rowindex[M.m_nlin]);
            copy(M);
        }

        ~FastSparseMatrix() { destroy(); }

        size_t nlin() const { return static_cast<size_t>(m_nlin); }
        size_t ncol() const { return static_cast<size_t>(m_ncol); }

        void write(std::ostream& f) const {
            const size_t nz = rowindex[m_nlin];
            f.write(reinterpret_cast<const char*>(&m_nlin), static_cast<std::streamsize>(sizeof(size_t)));
            f.write(reinterpret_cast<const char*>(&m_ncol), static_cast<std::streamsize>(sizeof(size_t)));
            f.write(reinterpret_cast<const char*>(&nz),     static_cast<std::streamsize>(sizeof(size_t)));
            f.write(reinterpret_cast<const char*>(tank),    static_cast<std::streamsize>(sizeof(double)*nz));
            f.write(reinterpret_cast<const char*>(js),      static_cast<std::streamsize>(sizeof(size_t)*nz));
            f.write(reinterpret_cast<const char*>(rowindex),static_cast<std::streamsize>(sizeof(size_t)*m_nlin));
        }

        void read(std::istream& f) {
            destroy();
            size_t nz;
            f.read(reinterpret_cast<char*>(&m_nlin), static_cast<std::streamsize>(sizeof(size_t)));
            f.read(reinterpret_cast<char*>(&m_ncol), static_cast<std::streamsize>(sizeof(size_t)));
            f.read(reinterpret_cast<char*>(&nz),     static_cast<std::streamsize>(sizeof(size_t)));
            alloc(m_nlin,m_ncol,nz);
            f.read(reinterpret_cast<char*>(tank),    static_cast<std::streamsize>(sizeof(double)*nz));
            f.read(reinterpret_cast<char*>(js),      static_cast<std::streamsize>(sizeof(size_t)*nz));
            f.read(reinterpret_cast<char*>(rowindex),static_cast<std::streamsize>(sizeof(size_t)*m_nlin));
        }

        double operator()(const size_t i,const size_t j) const {
            for (size_t k=rowindex[i]; k<rowindex[i+1]; ++k) {
                if (js[k]<j)
                    continue;
                if (js[k]==j)
                    return tank[k];
                break;
            }

            return 0.0;
        }

        double& operator()(const size_t i,const size_t j) {
            for (size_t k=rowindex[i]; k<rowindex[i+1]; ++k) {
                if (js[k]<j)
                    continue;
                if (js[k]==j)
                    return tank[k];
                break;
            }

            throw OpenMEEG::maths::BadSparseOperation("FastSparseMatrix : double& operator()(size_t i,size_t j) can't add element");
        }

        Vector operator*(const Vector& v) const {
            Vector result(m_nlin);
            result.set(0.0);
            double* pt_result = result.data();
            double* pt_vect   = v.data();

            for (size_t i=0; i<m_nlin; ++i) {
                double& total = pt_result[i];
                for (size_t j=rowindex[i]; j<rowindex[i+1]; ++j)
                    total += tank[j]*pt_vect[js[j]];
            }
            return result;
        }

        double& operator[](const size_t i) { return tank[i]; }

        void info() const {
            if (nlin()==0 && ncol()==0) {
                std::cout << "Empty Matrix" << std::endl;
                return;
            }

            std::cout << "Dimensions : " << nlin() << " x " << ncol() << std::endl;
            std::cout << *this;
        }

    private:

        void alloc(const size_t nl,const size_t nc,const size_t nz) {
            m_nlin       = nl;
            m_ncol       = nc;
            tank         = new double[nz];
            js           = new size_t[nz];
            rowindex     = new size_t[nl+1];
            rowindex[nl] = nz;
        }

        void destroy() {
            delete[] tank;
            delete[] js;
            delete[] rowindex;
        }

        void copy(const FastSparseMatrix& M) {
            alloc(M.m_nlin,M.m_ncol,M.rowindex[M.m_nlin]);
            memcpy(tank,M.tank,sizeof(double)*M.rowindex[M.m_nlin]);
            memcpy(js,M.js,sizeof(size_t)*M.rowindex[M.m_nlin]);
            memcpy(rowindex,M.rowindex,sizeof(size_t)*(M.m_nlin+1));
        }

        double* tank;
        size_t* js;
        size_t* rowindex;
        size_t  m_nlin;
        size_t  m_ncol;
    };

    inline std::ostream& operator<<(std::ostream& f,const FastSparseMatrix &M) {
        const size_t nz = M.rowindex[M.nlin()];
        f << M.nlin() << " " << M.ncol() << std::endl;
        f << nz << std::endl;
        for (size_t i=0; i<M.nlin(); ++i)
            for (size_t j=M.rowindex[i]; j<M.rowindex[i+1]; ++j)
                f << static_cast<unsigned long>(i) << "\t" << static_cast<unsigned long>(M.js[j]) << "\t" << M.tank[j] << std::endl;
        return f;
    }
}

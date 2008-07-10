// FILE: $Id$

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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

#ifndef MATLABIO_H
#define MATLABIO_H

#ifdef USE_MATIO
#include <matio.h>
#endif

#include <Exceptions.H>
#include "MathsIO.H"
#include "sparse_matrix.h"
#include "matrix.h"

namespace maths {

    struct MatlabIO: public MathsIOBase {

        const std::string& identity() const { return Identity; }
        const Suffixes&    suffixes() const { return suffs;    }

        bool identify(const std::string& buffer) const {
            if (buffer.size()<MagicTag.size())
                return false;
            return strncmp(buffer.c_str(),MagicTag.c_str(),MagicTag.size()) == 0;
        }

        bool known(const LinOp& linop) const {
            return linop.storageType()==LinOp::FULL || linop.storageType()==LinOp::SPARSE;
        }

#ifdef USE_MATIO

        void read_full(mat_t* mat,LinOp& linop) const throw(std::exception) {
            matvar_t* matvar = Mat_VarReadNext(mat);
            while (matvar!=NULL && (matvar->rank!=2 || matvar->data_type!=MAT_T_DOUBLE || matvar->class_type!=MAT_C_DOUBLE))
                matvar = Mat_VarReadNext(mat);
            if (matvar==NULL)
                throw maths::BadContent(identity(),"2D full double Matrix");

            Matrix& m = dynamic_cast<Matrix&>(linop);

            m.nlin() = matvar->dims[0];
            m.ncol() = matvar->dims[1];
            double* data = m.data();
            data = static_cast<double*>(matvar->data);
            matvar->mem_conserve = 1;
            Mat_VarFree(matvar);
            Mat_Close(mat);
        }

        void read_sparse(mat_t* mat,LinOp& linop) const throw(std::exception) {
            SparseMatrix& m = dynamic_cast<SparseMatrix&>(linop);

            matvar_t* matvar = Mat_VarReadNext(mat);

            while( matvar!=NULL && (matvar->rank!=2 || matvar->data_type!=MAT_T_DOUBLE || matvar->class_type!=MAT_C_SPARSE) )
                matvar = Mat_VarReadNext(mat);
            if (matvar==NULL)
                throw maths::BadContent(identity(),"2D sparse double Matrix");

            std::cout << "Using variable with name : " << matvar->name << std::endl;

            m.nlin() = matvar->dims[0];
            m.ncol() = matvar->dims[1];

            sparse_t* sparse = static_cast<sparse_t*>(matvar->data);
            size_t _nz = sparse->nzmax;

            double *data = static_cast<double*>(sparse->data);

            size_t current_col = 0;
            for (size_t k = 0; k < _nz; ++k) {
                size_t i = sparse->ir[k];
                double val = data[k];
                if (k < (size_t) sparse->jc[sparse->njc-1]) {
                    assert(current_col < (size_t) sparse->njc);
                    while((size_t) sparse->jc[current_col + 1] <= k) current_col++; // look for the last idx of jc such that jc[idx+1] > k
                    size_t j = current_col;
                    m(i,j) = val;
                }
            }
            matvar->mem_conserve = 1;
            Mat_VarFree(matvar);
            Mat_Close(mat);
        }

        void write_full(mat_t* mat,const LinOp& linop) const {
            const Matrix& m = dynamic_cast<const Matrix&>(linop);
            matvar_t* matvar;
            int dims[2] = { m.nlin(), m.ncol() };
            matvar = Mat_VarCreate("linop",MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,m.data(),0);
            Mat_VarWrite(mat,matvar,COMPRESSION_ZLIB);
            Mat_VarFree(matvar);
            Mat_Close(mat);
        }

        void write_sparse(mat_t* mat,const LinOp& linop) const {
            const SparseMatrix& m = dynamic_cast<const SparseMatrix&>(linop);

            SparseMatrix::Tank tank_inverted;
            for (SparseMatrix::const_iterator it=m.begin();it!=m.end();++it) {
                size_t i = it->first.first;
                size_t j = it->first.second;
                double val = it->second;
                tank_inverted[std::pair<size_t,size_t>(j,i)] = val;
            }

            double*      t  = new double[linop.size()];
            mat_int32_t* ir = new mat_int32_t[linop.size()];
            mat_int32_t* jc = new mat_int32_t[linop.ncol()+1];

            size_t cnt = 0;
            size_t current_col = (size_t)-1;

            for (SparseMatrix::const_iterator it = tank_inverted.begin(); it != tank_inverted.end(); ++it) {
                const size_t i   = it->first.second;
                const size_t j   = it->first.first;
                const double val = it->second;
                if (current_col != j) {
                    for (size_t k = current_col+1; k <= j; ++k)
                        jc[k] = cnt;
                    current_col = j;
                }
                ir[cnt] = i;
                t[cnt] = val;
                cnt++;
            }

            for (size_t k = current_col+1; k <= linop.ncol(); ++k)
                jc[k] = linop.size();

            int dims[2] = { linop.nlin(), linop.ncol() };
            matvar_t *matvar;
            sparse_t  sparse = {0,};
            sparse.nzmax = linop.size();
            sparse.nir   = linop.size();
            sparse.ir    = ir;
            sparse.njc   = linop.ncol()+1;
            sparse.jc    = jc;
            sparse.ndata = linop.size();
            sparse.data  = t;

            matvar = Mat_VarCreate("matrix",MAT_C_SPARSE,MAT_T_DOUBLE,2,dims,&sparse,MEM_CONSERVE);
            Mat_VarWrite(mat,matvar,COMPRESSION_ZLIB);
            Mat_VarFree(matvar);
            Mat_Close(mat);
        }

#endif

        void read(std::ifstream& is,LinOp& linop) const throw(std::exception) {
            #ifdef USE_MATIO
                if (is.is_open()) {
                    is.close();
                }
                mat_t* mat = Mat_Open(name().c_str(),MAT_ACC_RDONLY);
                if (mat) {
                    if (linop.storageType()==LinOp::SPARSE) {
                        read_sparse(mat,linop);
                    } else {
                        read_full(mat,linop);
                    }
                } else {
                    throw BadFileOpening(name(),maths::BadFileOpening::READ);
                }
            #else
                std::cerr << "You have to compile OpenMEEG with MATIO to load matlab files" << std::endl;
            #endif
        }

        void write(std::ofstream& os, const LinOp& linop) const throw(std::exception) {
            #ifdef USE_MATIO
                if (os.is_open()) {
                    std::cout << "closing stream" << std::endl;
                    os.close();
                    std::remove(name().c_str());
                }
                mat_t* mat = Mat_Open(name().c_str(),MAT_ACC_RDWR);
                if (mat) {
                    if (linop.storageType()==LinOp::SPARSE) {
                        write_sparse(mat,linop);
                    } else {
                        write_full(mat,linop);
                    }
                } else {
                    throw BadFileOpening(name(),maths::BadFileOpening::WRITE);
                }
            #else
                std::cerr << "You have to compile OpenMEEG with MATIO to save matlab files" << std::endl;
            #endif
        }

    private:

        MatlabIO(): MathsIOBase(MathsIO::ios()) { }
        ~MatlabIO() {};

        static Suffixes init() {
            Suffixes suffs;
            suffs.push_back("mat");
            return suffs;
        }

        static const MatlabIO    prototype;
        static const std::string MagicTag;
        static const Suffixes    suffs;
        static const std::string Identity;
    };
}

#endif  //! MATLABIO_H

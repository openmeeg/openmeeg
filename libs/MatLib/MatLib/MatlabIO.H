// FILE: $Id$

/*
Project Name : OpenMEEG

author            : $Author$
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

#include "MatrixIO.H"
#include "sparse_matrice_dcl.h"
// #include "matrice_dcl.h"

namespace Maths {

    struct MatlabIO: public MathsIOBase {

        const std::string& identity() const { return Identity; }
        const Suffixes&    suffixes() const { return suffs;    }

        bool identify(const std::string& buffer) const {
            if (buffer.size()<MagicTag.size())
                return false;
            return strncmp(buffer.c_str(),MagicTag.c_str(),MagicTag.size()) == 0;
        }

        bool known(const MatrixBase& matrix) const {
            return matrix.storageType()==MatrixBase::FULL || matrix.storageType()==MatrixBase::SPARSE;
        }

        void read_full(mat_t* mat,MatrixBase& matrix) const {
            matvar_t* matvar = Mat_VarReadNext(mat);
            while (matvar!=NULL && (matvar->rank!=2 || matvar->data_type!=MAT_T_DOUBLE || matvar->class_type!=MAT_C_DOUBLE))
                matvar = Mat_VarReadNext(mat);
            if (matvar==NULL)
                throw std::string("There is no 2D full double matrix in this Matlab file");

            matrice* M = dynamic_cast<matrice*>(&matrix);

            M->nlin() = matvar->dims[0];
            M->ncol() = matvar->dims[1];
            double* data = M->data();
            data = static_cast<double*>(matvar->data);
            matvar->mem_conserve = 1;
            Mat_VarFree(matvar);
            Mat_Close(mat);
        }

        void read_sparse(mat_t* mat,MatrixBase& matrix) const {
            sparse_matrice* M = dynamic_cast<sparse_matrice*>(&matrix);

            matvar_t* matvar = Mat_VarReadNext(mat);

            while( matvar!=NULL && (matvar->rank!=2 || matvar->data_type!=MAT_T_DOUBLE || matvar->class_type!=MAT_C_SPARSE) )
                matvar = Mat_VarReadNext(mat);
            if (matvar==NULL) {
                throw std::string("There is no 2D sparse double matrix in this Matlab file");
            }

            std::cout << "Using variable with name : " << matvar->name << std::endl;

            M->nlin() = matvar->dims[0];
            M->ncol() = matvar->dims[1];

            sparse_t* sparse = static_cast<sparse_t*>(matvar->data);
            size_t _nz = sparse->nzmax;

            double *data = static_cast<double*>(sparse->data);

            size_t current_col = 0;
            for(size_t k = 0; k < _nz; ++k)
            {
                size_t i = sparse->ir[k];
                double val = data[k];
                if(k < (size_t) sparse->jc[sparse->njc-1])
                {
                    assert(current_col < (size_t) sparse->njc);
                    while((size_t) sparse->jc[current_col + 1] <= k) // look for the last idx of jc such that jc[idx+1] > k
                    {
                        current_col++;
                    }
                    size_t j = current_col;
                    M->operator()(i,j) = val;
                }
            }
            matvar->mem_conserve = 1;
            Mat_VarFree(matvar);
            Mat_Close(mat);
        }

        void write_full(mat_t* mat,const MatrixBase& matrix) const {
            matrice* M = dynamic_cast<matrice*>(const_cast<MatrixBase*>(&matrix));
            matvar_t* matvar;
            int dims[2] = { M->nlin(), M->ncol() };
            matvar = Mat_VarCreate("matrix",MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,M->data(),0);
            Mat_VarWrite(mat,matvar,COMPRESSION_ZLIB);
            Mat_VarFree(matvar);
            Mat_Close(mat);
        }

        void write_sparse(mat_t* mat,const MatrixBase& matrix) const {
            sparse_matrice* M = dynamic_cast<sparse_matrice*>(const_cast<MatrixBase*>(&matrix));
            sparse_matrice::const_iterator it;
            
            sparse_matrice::Tank tank_inverted;
            
            for(it = M->begin(); it != M->end(); ++it) {
                size_t i = it->first.first;
                size_t j = it->first.second;
                double val = it->second;
                tank_inverted[std::pair<size_t,size_t>(j,i)] = val;
            }
            
            double* t = (double*)malloc(sizeof(double)*matrix.size());
            mat_int32_t  *ir = (mat_int32_t*)malloc(sizeof(mat_int32_t)*matrix.size());
            mat_int32_t  *jc = (mat_int32_t*)malloc(sizeof(mat_int32_t)*(matrix.ncol()+1));
            
            size_t cnt = 0;
            long current_col = -1;
            
            for(it = tank_inverted.begin(); it != tank_inverted.end(); ++it) {
                size_t i = it->first.second;
                size_t j = it->first.first;
                double val = it->second;
                if(current_col != j) {
                    for(size_t k = current_col+1; k <= j; ++k) {
                        jc[k] = cnt;
                    }
                    current_col = j;
                }
                ir[cnt] = i;
                t[cnt] = val;
                cnt++;
            }
            
            for(size_t k = current_col+1; k <= matrix.ncol(); ++k) {
                jc[k] = matrix.size();
            }
            
            int dims[2] = { matrix.nlin(), matrix.ncol() };
            matvar_t *matvar;
            sparse_t  sparse = {0,};
            sparse.nzmax = matrix.size();
            sparse.nir   = matrix.size();
            sparse.ir    = ir;
            sparse.njc   = matrix.ncol()+1;
            sparse.jc    = jc;
            sparse.ndata = matrix.size();
            sparse.data  = t;
            
            matvar = Mat_VarCreate("matrix",MAT_C_SPARSE,MAT_T_DOUBLE,2,dims,&sparse,MEM_CONSERVE);
            Mat_VarWrite(mat,matvar,COMPRESSION_ZLIB);
            Mat_VarFree(matvar);
            Mat_Close(mat);
        }

        void read(std::ifstream& is,MatrixBase& matrix) const {
            #ifdef USE_MATIO
                if(is.is_open()) {
                    is.close();
                }
                mat_t* mat = Mat_Open(name().c_str(),MAT_ACC_RDONLY);
                if (mat) {
                    if(matrix.storageType()==MatrixBase::SPARSE) {
                        read_sparse(mat,matrix);
                    } else {
                        read_full(mat,matrix);
                    }
                } else {
                    throw std::string("Unable to open : ")+name();
                }
            #else
                std::cerr << "You have to compile OpenMEEG with MATIO to load matlab files" << std::endl;
            #endif
        }

        void write(std::ofstream& os, const MatrixBase& matrix) const {
            #ifdef USE_MATIO
                if(os.is_open()) {
                    std::cout << "closing stream" << std::endl;
                    os.close();
                    std::remove(name().c_str());
                }
                mat_t* mat = Mat_Open(name().c_str(),MAT_ACC_RDWR);
                if (mat) {
                    if(matrix.storageType()==MatrixBase::SPARSE) {
                        write_sparse(mat,matrix);
                    } else {
                        write_full(mat,matrix);
                    }
                } else {
                    throw std::string("Unable to open : ")+name();
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

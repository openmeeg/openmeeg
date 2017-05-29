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

#include <matio.h>

#include <Exceptions.H>
#include "MathsIO.H"
#include "sparse_matrix.h"
#include "matrix.h"

namespace OpenMEEG {
    namespace maths {

        namespace details {
            template <typename T> class helper { };

            template <>
            struct helper<Vector> {
                static const int dim = 1;
                static size_t* dims(const Vector& m) {
                    static size_t dims[1];
                    dims[0] = m.nlin();
                    return dims;
                }
                static matvar_t* set_type(matvar_t* matvar,LinOpInfo& linop) {
                    linop.dimension()   = 1;
                    linop.storageType() = LinOpInfo::FULL;
                    linop.nlin()        = matvar->dims[0];
                    return matvar;
                }
                static bool good_type(const matvar_t* matvar) {
                    return matvar->rank==1 && matvar->class_type==MAT_C_DOUBLE && matvar->data_type==MAT_T_DOUBLE;
                }
                static const char message[];
            };

            template <>
            struct helper<Matrix> {
                static const int dim = 2;
                static size_t* dims(const Matrix& m) {
                    static size_t dims[2];
                    dims[0] = m.nlin();
                    dims[1] = m.ncol();
                    return dims;
                }
                static matvar_t* set_type(matvar_t* matvar,LinOpInfo& linop) {
                    linop.dimension()   = 2;
                    linop.storageType() = LinOpInfo::FULL;
                    linop.nlin()        = matvar->dims[0];
                    linop.ncol()        = matvar->dims[1];
                    return matvar;
                }
                static bool good_type(const matvar_t* matvar) {
                    return matvar->rank==2 && matvar->class_type==MAT_C_DOUBLE && matvar->data_type==MAT_T_DOUBLE;
                }
                static const char message[];
            };

            template <>
            struct helper<SymMatrix> {
                static matvar_t* saved;

                ~helper() { Mat_VarFree(saved); }

                static const int dim = 1;
                static size_t* dims(const Matrix& m) {
                    static size_t dims[1];
                    dims[0] = (m.nlin()*(m.nlin()+1))/2;
                    return dims;
                }
                static matvar_t* set_type(matvar_t* matvar,LinOpInfo& linop) {
                    static const void* size = "size";
                    linop.dimension()   = 2;
                    linop.storageType() = LinOpInfo::SYMMETRIC;
                    const matvar_t* sizevar = Mat_VarGetStructField(matvar,const_cast<void*>(size),MAT_BY_NAME,0);
                    linop.nlin() = linop.ncol() = *static_cast<int*>(sizevar->data);
                    return adjust(matvar);
                }
                static bool good_type(const matvar_t* matvar) {
                    return matvar->rank==2 && matvar->class_type==MAT_C_STRUCT && matvar->data_type==MAT_T_STRUCT;
                }
                static matvar_t* adjust(matvar_t* structvar) {
                    static const void* data = "data";
                    matvar_t* matvar = Mat_VarGetStructField(structvar,const_cast<void*>(data),MAT_BY_NAME,0);
                    saved = structvar;
                    return matvar;
                }
                static const char message[];
            };

            template <>
            struct helper<SparseMatrix> {
                static const int dim = 2;
                static matvar_t* set_type(matvar_t* matvar,LinOpInfo& linop) {
                    linop.dimension()   = 2;
                    linop.storageType() = LinOpInfo::SPARSE;
                    linop.nlin()        = matvar->dims[0];
                    linop.ncol()        = matvar->dims[1];
                    return matvar;
                }
                static bool good_type(const matvar_t* matvar) {
                    return matvar->rank==2 && matvar->data_type==MAT_T_DOUBLE && matvar->class_type==MAT_C_SPARSE;
                }
                static const char message[];
            };
        }

        struct OPENMEEGMATHS_EXPORT MatlabIO: public MathsIOBase {

            const std::string& identity() const { return Identity; }
            const Suffixes&    suffixes() const { return suffs;    }

            bool identify(const std::string& buffer) const {
                if (buffer.size()<MagicTag.size())
                    return false;
                return strncmp(buffer.c_str(),MagicTag.c_str(),MagicTag.size()) == 0;
            }

            bool known(const LinOp& linop) const {
                const unsigned st = linop.storageType();
                return st==LinOpInfo::FULL || st==LinOpInfo::SYMMETRIC || st==LinOpInfo::SPARSE;
            }

            LinOpInfo info(std::ifstream& is) const {
                if (is.is_open())
                    is.close();

                LinOpInfo linop;
                mat_t* mat = Mat_Open(name().c_str(),MAT_ACC_RDONLY);
                if (mat) {
                    try { read_header<Vector>(mat,linop);       return finish(mat,linop); } catch(...) { Mat_Rewind(mat); }
                    try { read_header<Matrix>(mat,linop);       return finish(mat,linop); } catch(...) { Mat_Rewind(mat); }
                    try { read_header<SymMatrix>(mat,linop);    return finish(mat,linop); } catch(...) { Mat_Rewind(mat); }
                    try { read_header<SparseMatrix>(mat,linop); return finish(mat,linop); } catch(...) { }
                    Mat_Close(mat);
                    throw ImpossibleObjectIdentification(name());
                } else {
                    throw BadFileOpening(name(),maths::BadFileOpening::READ);
                }

                return linop;
            }

            void read(std::ifstream& is,LinOp& linop) const {
                if (is.is_open()) {
                    is.close();
                }
                mat_t* mat = Mat_Open(name().c_str(),MAT_ACC_RDONLY);
                if (!mat)
                    throw BadFileOpening(name(),maths::BadFileOpening::READ);

                if (linop.dimension()==1)
                    read<Vector>(mat,linop);
                else {
                    switch (linop.storageType()) {
                        case LinOp::SPARSE:
                            read_sparse(mat,linop);
                            break;
                        case LinOp::SYMMETRIC:
                            read<SymMatrix>(mat,linop);
                            break;
                        case LinOp::FULL:
                            read<Matrix>(mat,linop);
                            break;
                        default: {
                            Mat_Close(mat);
                            throw BadStorageType(name());
                        }
                    }
                }

                Mat_Close(mat);
            }

            void write(std::ofstream& os, const LinOp& linop) const {
                if (os.is_open()) {
                    os.close();
                    std::remove(name().c_str());
                }
                mat_t* mat = Mat_CreateVer(name().c_str(),NULL,MAT_FT_MAT73);
                if (!mat)
                    throw BadFileOpening(name(),maths::BadFileOpening::WRITE);

                if (linop.dimension()==1)
                    write<Vector>(mat,linop);
                else {
                    switch (linop.storageType()) {
                        case LinOp::SPARSE:
                            write_sparse(mat,linop);
                            break;
                        case LinOp::SYMMETRIC:
                            write_symmetric(mat,linop);
                            break;
                        case LinOp::FULL:
                            write<Matrix>(mat,linop);
                            break;
                        default:
                            throw BadStorageType(name());
                    }
                }

                Mat_Close(mat);
            }

        private:

            static void matiologfunc(const int log_level,char* message) {

                //  Ugly copy of values in matio/io.c, but those are not exported...

                enum { LOG_LEVEL_ERROR=1, LOG_LEVEL_CRITICAL=2, LOG_LEVEL_WARNING=4, LOG_LEVEL_MESSAGE=8, LOG_LEVEL_DEBUG=16 };

                if (log_level & (LOG_LEVEL_CRITICAL|LOG_LEVEL_ERROR))
                    throw MatioError(message);

                std::string error;
                if (log_level==LOG_LEVEL_WARNING)
                    error = "Warning";
                else if (log_level==LOG_LEVEL_DEBUG)
                    error = "Debug";
                else if (log_level==LOG_LEVEL_MESSAGE)
                    error = "Message";
                else
                    error = "Unknown error";

                std::cerr << error << message << std::endl;
            }

            MatlabIO(): MathsIOBase(5) { Mat_LogInitFunc("OpenMEEG",matiologfunc); }

            ~MatlabIO() {};

            template <typename TYPE>
            matvar_t* read_header(mat_t* mat,LinOpInfo& linop) const {
                matvar_t* matvar = Mat_VarReadNext(mat);
                while (matvar!=NULL && !details::helper<TYPE>::good_type(matvar)) {
                    Mat_VarFree(matvar);
                    matvar = Mat_VarReadNext(mat);
                }
                if (matvar==NULL)
                    throw maths::BadContent(identity(),details::helper<TYPE>::message);
                return details::helper<TYPE>::set_type(matvar,linop);
            }

            template <typename TYPE>
            void read(mat_t* mat,LinOp& linop) const {
                matvar_t* matvar = read_header<TYPE>(mat,linop);
                TYPE& O = dynamic_cast<TYPE&>(linop);
                O.reference_data(static_cast<double*>(matvar->data));
                matvar->mem_conserve = 1;
                Mat_VarFree(matvar);
            }

            void read_sparse(mat_t* mat,LinOp& linop) const {
                matvar_t*     matvar = read_header<SparseMatrix>(mat,linop);
                SparseMatrix& m      = dynamic_cast<SparseMatrix&>(linop);
                std::cout << "Using variable with name : " << matvar->name << std::endl;

                mat_sparse_t* sparse = static_cast<mat_sparse_t*>(matvar->data);
                size_t _nz = sparse->nzmax;

                double *data = static_cast<double*>(sparse->data);

                size_t current_col = 0;
                for (size_t k = 0; k < _nz; ++k) {
                    size_t i = sparse->ir[k];
                    double val = data[k];
                    if (k < (size_t) sparse->jc[sparse->njc-1]) {
                        om_assert(current_col < (size_t) sparse->njc);
                        while((size_t) sparse->jc[current_col + 1] <= k) current_col++; // look for the last idx of jc such that jc[idx+1] > k
                        size_t j = current_col;
                        m(i,j) = val;
                    }
                }
                matvar->mem_conserve = 1;
                Mat_VarFree(matvar);
            }

            template <typename TYPE>
            void write(mat_t* mat,const LinOp& linop) const {
                const TYPE& m = dynamic_cast<const TYPE&>(linop);
                size_t* dims = details::helper<TYPE>::dims(m);
                matvar_t* matvar = Mat_VarCreate("linop",MAT_C_DOUBLE,MAT_T_DOUBLE,details::helper<TYPE>::dim,dims,m.data(),MAT_F_DONT_COPY_DATA);
                Mat_VarWrite(mat,matvar,MAT_COMPRESSION_ZLIB);
                Mat_VarFree(matvar);
            }

            void write_symmetric(mat_t* mat,const LinOp& linop) const {
                const SymMatrix& m  = dynamic_cast<const SymMatrix&>(linop);
                size_t dims[2] = { m.size(),1 };
                size_t dims1[2] = { 1, 1 };
                size_t size[1]  = { m.nlin() };
                matvar_t** matvar = new matvar_t*[3];
                matvar[0] = Mat_VarCreate("size",MAT_C_UINT32,MAT_T_UINT32,2,dims1,size,0);
                matvar[1] = Mat_VarCreate("data",MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,m.data(),MAT_F_DONT_COPY_DATA);
                matvar[2] = NULL;
                matvar_t* structvar = Mat_VarCreate("symmatrix",MAT_C_STRUCT,MAT_T_STRUCT,2,dims1,matvar,0);
                Mat_VarWrite(mat,structvar,MAT_COMPRESSION_ZLIB);
                delete[] matvar;
                Mat_VarFree(structvar);
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

                const int sz = linop.size();
                const int num = linop.ncol()+1;

                double* t  = new double[sz];
                int*    ir = new int[sz];
                int*    jc = new int[num];

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
                    jc[k] = sz;

                size_t dims[2] = { linop.nlin(), linop.ncol() };
                matvar_t *matvar;
                mat_sparse_t  sparse = { sz, ir, sz, jc, num, sz, t};

                matvar = Mat_VarCreate("matrix",MAT_C_SPARSE,MAT_T_DOUBLE,2,dims,&sparse,MAT_F_DONT_COPY_DATA);
                Mat_VarWrite(mat,matvar,MAT_COMPRESSION_ZLIB);
                Mat_VarFree(matvar);
                delete[] t;
                delete[] ir;
                delete[] jc;
            }

            static LinOpInfo finish(mat_t* mat,LinOpInfo linop) {
                Mat_Close(mat);
                return linop;
            }

            static Suffixes init() {
                Suffixes suffxes;
                suffxes.push_back("mat");
                return suffxes;
            }

            static const MatlabIO    prototype;
            static const std::string MagicTag;
            static const Suffixes    suffs;
            static const std::string Identity;
        };
    }
}

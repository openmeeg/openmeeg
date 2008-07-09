// FILE: $Id:$

/*
Project Name : OpenMEEG

author            : $Author: gramfort $
version           : $Revision: 222 $
last revision     : $Date: 2008-04-08 08:14:41 +0200 (Tue, 08 Apr 2008) $
modified by       : $LastChangedBy: gramfort $
last modified     : $LastChangedDate: 2008-04-08 08:14:41 +0200 (Tue, 08 Apr 2008) $

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

#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include <IOInit.H>

struct MatrixBase {

    typedef enum { FULL, SYMMETRIC, SPARSE } MatrixStorageType;

    MatrixBase(const size_t n,const size_t m,const MatrixStorageType st=FULL): num_lines(n),num_cols(m),storage(st)  { }
    ~MatrixBase() { }

    size_t& nlin()       { return num_lines; }
    size_t  nlin() const { return num_lines; }

    size_t& ncol()       { return num_cols; }
    size_t  ncol() const { return num_cols; }

    size_t size() const = 0;
    MatrixStorageType storageType() const { return storage; }

    virtual std::string identity() const = 0;
    { return "Stupid Binary Format"; }

    virtual double* data() = 0;
    virtual void*   alloc_data() = 0;

private:

    size_t            num_lines;
    size_t            num_cols;
    MatrixStorageType storage;
};

namespace Internal {

    //  Read a few bytes to figure out the file format and put them back into the stream.

    static const unsigned maxtagsize = 32;

    static const char*
    ReadTag(std::istream& is) throw(std::string) {

        static char  buffer[maxtagsize];

        try {
            is.read(buffer,maxtagsize);
        } catch(...) {
            throw "BadHeader(is)";
        }

        for(int i=maxtagsize-1;i>=0;--i)
            is.putback(buffer[i]);

        return buffer;
    }
}

std::istream& operator>>(std::istream& is,MatrixBase& matrix) throw(std::string) {
    IOInit<std::istream> init(is);
    const char* buffer = Internal::ReadTag(is);

    if (MatrixIO* io = MatrixIO::default_io()!=0) {
        if (io->identify(buffer)) {
            io->read(is,matrix);
            return is;
        }
    } else {
        for (MatrixIO::IOs::const_iterator io=MatrixIO::ios().begin();io!=MatrixIO::ios().end();++io)
            if (io->identify(buffer)) {
                io->read(is,matrix);
                return is;
            }
    }
    throw "Horreur!";
}

std::ostream& operator<<(std::ostream& is,const MatrixBase& matrix) throw(std::string) {
    if (MatrixIO* io = MatrixIO::default_io()!=0) {
        if (io->known(matrix)) {
            io->write(matrix);
            return is;
        }
    } else {
        for (MatrixIO::IOs::const_iterator io=MatrixIO::ios().begin();io!=MatrixIO::ios().end();++io)
            if (io->matrix(matrix)) {
                io->write(matrix);
                return is;
            }
    }
    throw "Horreur!";
}

#endif  //! BASE_MATRIX_H

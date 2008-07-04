// FILE: $Id:$

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

#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include "MatLibConfig.h"

struct MatrixBase {

    typedef enum { FULL, SYMMETRIC, SPARSE } MatrixStorageType;
    typedef enum { ONE, TWO } Dimension;

    MatrixBase() { }
    MatrixBase(const size_t m,const size_t n,const MatrixStorageType st, const Dimension d):
        num_lines(m),num_cols(n),storage(st),dim(d)  { }

    virtual ~MatrixBase() {};

    size_t& nlin()       { return num_lines; }
    size_t  nlin() const { return num_lines; }

    size_t& ncol()       { return num_cols; }
    size_t  ncol() const { return num_cols; }

    virtual size_t size() const = 0;
    MatrixStorageType storageType() const { return storage; }
    Dimension dimension() const { return dim; }
    virtual void   info() const = 0;

protected:

    size_t            num_lines;
    size_t            num_cols;
    MatrixStorageType storage;
    Dimension         dim;

};

#endif  //! BASE_MATRIX_H

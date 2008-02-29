/* FILE: $Id$ */

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

#ifndef H_genericMatrix

#define H_genericMatrix
#include <iostream>
#include "LinOp.h"

template<class T> class TgenericMatrix : public LinOp
{
public:
    virtual size_t nlin() const=0;
    virtual size_t ncol() const=0;
    virtual T operator()(size_t i,size_t j) const=0;
    virtual T& operator()(size_t i,size_t j)=0;
    virtual void saveTxt( const char *filename ) const=0;
    virtual void saveBin( const char *filename ) const=0;
    virtual void loadTxt( const char *filename )=0;
    virtual void loadBin( const char *filename )=0;
    virtual void write(std::ostream& f) const =0;
    virtual void read(std::istream& f) =0;
    virtual std::ostream& operator>>(std::ostream& f) const =0;

// #ifdef USE_MATIO
//     virtual void saveMat( const char *filename ) const=0;
//     virtual void loadMat( const char *filename )=0;
// #endif
    
};

typedef TgenericMatrix<double> genericMatrix;


#endif


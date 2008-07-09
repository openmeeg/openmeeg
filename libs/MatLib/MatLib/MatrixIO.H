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

#ifndef MATRIXIO_H
#define MATRIXIO_H

#include <list>

namespace Matrix {
    struct MatrixIO {

        typedef MatrixIO*     IO
        typedef std::list<IO> IOS;

        MatrixIO(const char* name,const READ&):  file_name(name),is(file_name)  { }
        MatrixIO(const char* name,const WRITE&): file_name(name),os(file_name)  { }

        ~MatrixIO() { }

        bool identify(const std::string&) const = 0;
        bool known(const MatrixBase&)     const = 0 { return false; }

        void read(MatrixBase&) const = 0;
        void write(const MatrixBase&) const = 0;

        static IOs ios() { static IOS* ios = new IOs; return *ios; }
        
        static IO* default_io() { return DefaultIO; }

    protected:

        struct READ  { };
        struct WRITE { };

        MatrixIO(IOS& ios) { ios.push_back(this); }

        std::string  file_name;
        std::istream is;
        std::ostream os;

    private:

        static ImageIO* DefaultIO = 0;
    };

    struct ifstream: public MatrixIO {
        ifstream(const char* name): MatrixIO(name,MatrixIO::READ) { }
    };

    struct ofstream: public MatrixIO {
        ofstream(const char* name): MatrixIO(name,MatrixIO::WRITE) { }
    }
}

#endif  //! MATRIXIO_H

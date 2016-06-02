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

#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <string>
#include <cstring>

#include "linop.h"

#include <Exceptions.H>

#ifdef WIN32
#pragma warning( disable : 4290)    //MSVC warning C4290
#endif

namespace OpenMEEG {
    namespace maths {

        class OPENMEEGMATHS_EXPORT MathsIOBase;

        //  Quite similar to ImageIO, how to fuse this ?

        struct OPENMEEGMATHS_EXPORT MathsIO {

            typedef MathsIOBase* IO;
            typedef std::set<IO> IOs;

        private:

            static IO DefaultIO;

        public:

            MathsIO(const unsigned pr): priority(pr) { };
            MathsIO(const char* file): file_name(file),priority(0) { }
            virtual ~MathsIO() {};

            static IOs& ios() { static IOs* ios = new IOs; return *ios; }

            static IO default_io() { return DefaultIO; }

            static bool permanent;

            const std::string& name() const { return file_name; }

            void setName(const std::string& n) { file_name = n; }

            //  This is very similar with ImageIO, how to commonize it ?

            //  Handle the default io.

            static IO GetCurrentFormat() {
                IO tmp = DefaultIO;
                if (!permanent)
                    DefaultIO = 0;
                return tmp;
            }

            static void SetCurrentFormat(IO io,const bool perm=false) {
                DefaultIO = io;
                permanent = perm;
            }

            static const IO& format(const std::string&);
            static const IO& format_from_suffix(const std::string&);

            static void SetCurrentFormat(const std::string& fmt,const bool perm) {
                if (fmt=="default") {
                    SetCurrentFormat(0,perm);
                    return;
                }
                SetCurrentFormat(format(fmt),perm);
            }

            static void SetCurrentFormatFromSuffix(const std::string& name,const bool perm) {
                SetCurrentFormat(format_from_suffix(name),perm);
            }

            operator bool() const {
                std::ifstream ifs(file_name.c_str());
                return static_cast<bool>(ifs);
            }

        protected:

            std::string    file_name;
            const unsigned priority;
        };

        class OPENMEEGMATHS_EXPORT MathsIOBase: public MathsIO {

            typedef MathsIO base;

        public:

            typedef base::IOs              IOs;
            typedef std::list<std::string> Suffixes;

            virtual const std::string& identity() const = 0;
            virtual const Suffixes&    suffixes() const = 0;

            virtual bool identify(const std::string&) const = 0;
            virtual bool known(const LinOp&) const = 0;

            virtual LinOpInfo info(std::ifstream&) const = 0;

            virtual void read(std::ifstream&,LinOp&) const = 0;
            virtual void write(std::ofstream&,const LinOp&) const = 0;

            virtual bool known_suffix(const char* suffix)  const throw() {
                const Suffixes& suffs = suffixes();
                for (Suffixes::const_iterator i=suffs.begin();i!=suffs.end();++i) {
                    if (!strcmp(suffix,i->c_str()))
                        return true;
                }
                return false;
            }

            bool operator<(const MathsIOBase& io) const { return priority<io.priority; }

        protected:

            MathsIOBase(const unsigned pr): MathsIO(pr) { base::ios().insert(this); }
            ~MathsIOBase() {};
        };

        typedef MathsIO ifstream;
        typedef MathsIO ofstream;

        OPENMEEGMATHS_EXPORT maths::ifstream& operator>>(maths::ifstream&,LinOp&);
        OPENMEEGMATHS_EXPORT maths::ofstream& operator<<(maths::ofstream&,const LinOp&);

        // The manip format() used to specify explicitely a format.
        // Almost similar to Images::format. How to fuse those.

        struct OPENMEEGMATHS_EXPORT format {

            typedef enum { FromFormatName, FromSuffix } FormatSpec;

            format(const char* id,       const FormatSpec sp=FromFormatName,const bool perm=false): identity(id),type(sp),permanent(perm) { }
            format(const std::string& id,const FormatSpec sp=FromFormatName,const bool perm=false): identity(id),type(sp),permanent(perm) { }

            friend maths::ifstream& operator>>(maths::ifstream&,const format&);
            friend maths::ofstream& operator<<(maths::ofstream&,const format&);

        private:

            void set() const {
                switch (type) {
                    case FromFormatName:
                        MathsIO::SetCurrentFormat(identity,permanent);
                        break;
                    case FromSuffix:
                        MathsIO::SetCurrentFormatFromSuffix(identity,permanent);
                        break;
                }
            }

            const std::string identity;
            const FormatSpec  type;
            const bool        permanent;
        };

        inline maths::ifstream& operator>>(maths::ifstream &is,const format& f) { f.set(); return is; }
        inline maths::ofstream& operator<<(maths::ofstream &os,const format& f) { f.set(); return os; }

        OPENMEEGMATHS_EXPORT LinOpInfo info(const char* name);
    }
}

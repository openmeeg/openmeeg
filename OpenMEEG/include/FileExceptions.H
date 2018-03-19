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
#include <stdexcept>  //  Needs class 'std::runtime_error'

//  Exception classes defined in this file (indentation shows derivation hierarchy) :

//  in namespace std::io_except :
//  io_error
//    read_error
//    write_error
//    in namespace std::io_except::file_except :
//    file_error
//      file_open_error
//        input_file_open_error
//        output_file_open_error

//  Define macro LANGUAGE with value 0 for english messages, or 1 for french messages.
//  Default is english. LANGUAGE is undefined at end of this file if it wasn't defined
//  before.

#ifndef LANGUAGE
#define LANGUAGE 0
#else
#define LANGUAGE_WAS_DEF
#endif

#if LANGUAGE < 0 || LANGUAGE > 1
#error : Invalid value for macro LANGUAGE
#endif

#if LANGUAGE == 0
#define UNEXP_EOF        string("Unexpected end of file")
#define FORMAT_ERR       string("Format error")
#define WRITE_ERR        string("Write error")
#define OPEN_ERR         string("Error opening file")
#define I_OPEN_ERR       string("Error opening input file")
#define O_OPEN_ERR       string("Error opening output file")
#define I_OPEN_ERR_CAUSE string("File read protected or doesn't exist")
#define O_OPEN_ERR_CAUSE string("File write protected")
#endif

#if LANGUAGE == 1
#define UNEXP_EOF        string("Fin de fichier inattendue")
#define FORMAT_ERR       string("Erreur de format")
#define WRITE_ERR        string("Erreur d'écriture")
#define OPEN_ERR         string("Erreur d'ouverture du fichier")
#define I_OPEN_ERR       string("Erreur d'ouverture en lecture du fichier")
#define O_OPEN_ERR       string("Erreur d'ouverture en écriture du fichier")
#define I_OPEN_ERR_CAUSE string("Fichier protégé en lecture ou inexistant")
#define O_OPEN_ERR_CAUSE string("Fichier protégé en écriture")
#endif

namespace std {
    namespace io_except {

        //  General input/output error. Requires a string describing the error.

        class io_error: public std::runtime_error {
        public:
            io_error(const std::string& what_arg):std::runtime_error(what_arg) { }
            ~io_error() throw() { }
        };

        //  Read error. Requires an istream object.

        class read_error: public io_error {
        public:
            read_error(const std::istream& is):io_error(is.eof() ? UNEXP_EOF : FORMAT_ERR) { }
        };

        //  Write error.

        class write_error: public io_error {
        public:
            write_error():io_error(WRITE_ERR) { }
        };

        namespace file_except {

            //  General input/output error. Requires a string describing the error.

            class file_error: public io_error {
            public:
                file_error(const std::string& what_arg):io_error(what_arg) { }
            };

            //  Error in opening a file.
            //  Optionnal arguments: the file name and the file open mode (in/out).

            class file_open_error: public file_error {
            public:

                file_open_error():file_error(std::string(OPEN_ERR)) { }
                file_open_error(const std::string& file_name):file_error(std::string(OPEN_ERR + " " + file_name)) { }

                file_open_error(const std::ios::openmode mode):
                    file_error(std::string((mode==ios::in) ? I_OPEN_ERR+"\n"+I_OPEN_ERR_CAUSE : O_OPEN_ERR+"\n"+O_OPEN_ERR_CAUSE)) { }

                file_open_error(const std::string& file_name,const std::ios::openmode mode):
                    file_error(std::string((mode==ios::in) ? I_OPEN_ERR+" "+file_name+"\n"+I_OPEN_ERR_CAUSE :
                                                        O_OPEN_ERR+" "+file_name+"\n"+O_OPEN_ERR_CAUSE)) { }
            };

            //  Error in opening an input file. Optionnal argument : the file name.

            class input_file_open_error: public file_open_error {
            public:

                input_file_open_error():file_open_error(ios::in) { }
                input_file_open_error(const std::string& file_name):file_open_error(file_name,ios::in) { }
            };

            //  Error in opening an output file. Optionnal argument: the file name.

            class output_file_open_error: public file_open_error {
            public:

                output_file_open_error():file_open_error(ios::out) { }
                output_file_open_error(const std::string& file_name):file_open_error(file_name,ios::out) { }
            };
        }
    }
}

//  Undefine local macros.
#undef UNEXP_EOF
#undef FORMAT_ERR
#undef WRITE_ERR
#undef OPEN_ERR
#undef I_OPEN_ERR
#undef O_OPEN_ERR
#undef I_OPEN_ERR_CAUSE
#undef O_OPEN_ERR_CAUSE

#ifndef LANGUAGE_WAS_DEF
#undef LANGUAGE
#else
#undef LANGUAGE_WAS_DEF
#endif

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

#include <string>
#include <sstream>
#include <iostream>
#include <exception>

namespace OpenMEEG {

    //TODO  remove all unused exception
#if 0 
    struct ReturnCode {
        static const unsigned code = 
    };
#endif

    typedef enum { UNEXPECTED = 128, IO_EXCPT, OPEN_ERROR,
                   BAD_FMT, NO_SUFFIX, NON_MATCH_FMT, BAD_HDR, BAD_DATA, WRONG_FILE_FMT, BAD_DIM, UNKN_DIM, BAD_SIZE_SPEC, UNKN_PIX, UNKN_PIX_TYPE,
                   UNKN_FILE_FMT, UNKN_FILE_SUFFIX, UNKN_NAMED_FILE_FMT, NON_MATCH_NAMED_FILE_FMT, NO_FILE_FMT,
                   BAD_PLGIN_LIST, BAD_PLGIN_FILE, BAD_PLGIN, ALREADY_KN_TAG, NON_EXISTING_DOMAIN,
                   NO_IMG_ARG, DIFF_IMG } ExceptionCode;

    class Exception: public std::exception {
    public:

        Exception(const std::string& w): whatstring(std::string("Mesh reader::Exception: ")+w) { }

        virtual ~Exception() throw() { }

        virtual const char*   what () const throw() { return whatstring.c_str(); }
        virtual ExceptionCode code()  const throw() = 0;

    private:

        const std::string whatstring;
    };
    
    struct UnexpectedException: public Exception {

        UnexpectedException(const std::string& func,const std::string& file,const unsigned line): Exception(message(func,file,line)) { }

        ExceptionCode code() const throw() { return UNEXPECTED; }

    private:

        static std::string message(const std::string& func,const std::string& file,const unsigned line) {
            std::ostringstream ost;
            ost << "Unexpected error in " << func << " at " << file << ':' << line << '.';
            return ost.str();
        }
    };

    struct IOException: public Exception {

        IOException(const std::string& str): Exception(str) { }

        template <typename CharT,typename Traits>
        IOException(std::basic_ios<CharT,Traits>& ios,const std::string& str): Exception(str) { ios.setstate(std::ios::failbit); }

        ExceptionCode code() const throw() { return IO_EXCPT; }
    };

    struct OpenError: public IOException {
        OpenError(const std::string& fmtname): IOException(message(fmtname)) { }

        ExceptionCode code() const throw() { return OPEN_ERROR; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Cannot open file "+fmtname+" !"); }   
    };

    struct BadFormat: public IOException {

        BadFormat(const std::string& fmtname):                  IOException(message(fmtname))    { }
        BadFormat(std::istream& is,const std::string& fmtname): IOException(is,message(fmtname)) { }

        ExceptionCode code() const throw() { return BAD_FMT; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Unable to read the input as a "+fmtname+" image file."); }
    };

    struct BadInterface: public Exception {

        BadInterface(const std::string& str): Exception(message(str))    { }

        ExceptionCode code() const throw() { return BAD_FMT; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Unable to find interface " + fmtname + "."); }
    };

    struct BadDomain: public Exception {

        BadDomain(const std::string& str): Exception(message(str))    { }

        ExceptionCode code() const throw() { return BAD_FMT; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Unable to find domain " + fmtname + "."); }
    };

    struct BadInterfaceIndex: public Exception {

        BadInterfaceIndex(const std::string& str): Exception(message(str))    { }

        ExceptionCode code() const throw() { return BAD_FMT; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Unable to read the input as a "+fmtname+" image file."); }
    };

    struct NoSuffix: public IOException {
        NoSuffix(const std::string& name): IOException(std::string("No identifiable suffix in name ")+name) { }

        ExceptionCode code() const throw() { return NO_SUFFIX; }
    };

    struct NonMatchingFormat: public IOException {

        NonMatchingFormat(const std::string& fmtname):                  IOException(message(fmtname))    { }
        NonMatchingFormat(std::ostream& os,const std::string& fmtname): IOException(os,message(fmtname)) { }

        ExceptionCode code() const throw() { return NON_MATCH_FMT; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Unable to save the image in the ")+fmtname+" format."; }
    };

    struct BadHeader: public IOException {

        BadHeader(const std::string& fmtname=""):                  IOException(message(fmtname))    { }
        BadHeader(std::istream& is,const std::string& fmtname=""): IOException(is,message(fmtname)) { }

        ExceptionCode code() const throw() { return BAD_HDR; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Bad ")+fmtname+" file header."; }
    };

    struct BadData: public IOException {

        BadData(const std::string& fmtname):                  IOException(message(fmtname))    { }
        BadData(std::istream& is,const std::string& fmtname): IOException(is,message(fmtname)) { }

        ExceptionCode code() const throw() { return BAD_DATA; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Bad ")+fmtname+" file data."; }
    };

    struct WrongFileFormat: public IOException {

        WrongFileFormat(const std::string& file): IOException(std::string("Wrong format for file ")+file+".") { }

        ExceptionCode code() const throw() { return WRONG_FILE_FMT; }
    };

    struct UnknownFileFormat: public IOException {

        UnknownFileFormat(const std::string& fmt): IOException(std::string("Unknown ")+fmt+" format.") { }

        template <typename CharT,typename Traits>
        UnknownFileFormat(std::basic_ios<CharT,Traits>& ios): IOException(ios,std::string("Unknown file format.")) { }

        ExceptionCode code() const throw() { return UNKN_FILE_FMT; }
    };

    struct UnknownFileSuffix: public IOException {
        UnknownFileSuffix(const std::string& fmt): IOException(std::string("Unknown ")+fmt+" suffix.") { }
        ExceptionCode code() const throw() { return UNKN_FILE_SUFFIX; }
    };

    struct UnknownNamedFileFormat: public Exception {
        UnknownNamedFileFormat(const std::string& name): Exception(std::string("Unknown format for file "+name+".")) { }
        ExceptionCode code() const throw() { return UNKN_NAMED_FILE_FMT; }
    };

    struct NonMatchingNamedFileFormat: public Exception {
        NonMatchingNamedFileFormat(const std::string& name,const std::string& fmt): Exception(std::string("Cannot save the image "+name+" as a "+fmt+" file.")) { }

        ExceptionCode code() const throw() { return NON_MATCH_NAMED_FILE_FMT; }
    };

    struct NoMatchingFileFormat: public IOException {
        NoMatchingFileFormat(std::ostream& os): IOException(os,std::string("No IO is able to write this type of image.")) { }

        ExceptionCode code() const throw() { return NO_FILE_FMT; }
    };

    template <typename T>
    struct NonExistingDomain: public std::exception {
        NonExistingDomain(const std::string& dom,const T& id): domain(dom), interface(id) { }

        virtual ~NonExistingDomain() throw() { }

        const char* what() const throw() {
            std::ostringstream oss;
            oss << "Definition of domain " << domain << ": interface " << interface << " does not exist";
            return oss.str().c_str();
        }

        ExceptionCode code() const throw() { return NON_EXISTING_DOMAIN; }

        const std::string domain;
        const T&          interface;
    };

    //  Revisit (TODO).

    struct NoImageArgument: public Exception {
        NoImageArgument(): Exception("No image argument.") { }

        ExceptionCode code() const throw() { return NO_IMG_ARG; }
    };

    struct DifferentImages: public Exception {
        DifferentImages(): Exception("Different image types.") { }

        ExceptionCode code() const throw() { return DIFF_IMG; }
    };
}

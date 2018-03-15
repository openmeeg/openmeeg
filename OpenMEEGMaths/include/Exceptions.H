#pragma once

//  Largely derivate work from Images/src/Exceptions.H (merge if possible).

#include <string>
#include <sstream>
#include <iostream>
#include <exception>

#include "OpenMEEGMaths_Export.h"

namespace OpenMEEG {
    namespace maths {

        typedef enum { UNEXPECTED = 128, IO_EXCPT,
                       BAD_FILE, BAD_FILE_OPEN, BAD_CONTENT, NO_SUFFIX, BAD_HDR, BAD_DATA, BAD_VECT, UNKN_DIM, BAD_SYMM_MAT,
                       BAD_STORAGE_TYPE, NO_IO, MATIO_ERROR, UNKN_FILE_FMT, UNKN_FILE_SUFFIX, NO_FILE_FMT, UNKN_NAMED_FILE_FMT,
                       IMPOSSIBLE_IDENTIFICATION } ExceptionCode;


        class OPENMEEGMATHS_EXPORT Exception: public std::exception {
        public:

            Exception(const std::string& w): whatstring(std::string("Exception: ")+w) { }

            virtual ~Exception() throw() { }

            virtual const char*   what () const throw() { return whatstring.c_str(); }
            virtual ExceptionCode code()  const throw() = 0;

        private:

            std::string whatstring;
        };
        
        struct OPENMEEGMATHS_EXPORT UnexpectedException: public Exception {

            UnexpectedException(const std::string& func,const std::string& file,const unsigned line): Exception(message(func,file,line)) { }

            ExceptionCode code() const throw() { return UNEXPECTED; }

        private:

            static std::string message(const std::string& func,const std::string& file,const unsigned line) {
                std::ostringstream ost;
                ost << "Unexpected error in " << func << " at " << file << ':' << line << '.';
                return ost.str();
            }
        };

        struct OPENMEEGMATHS_EXPORT IOException: public Exception {

            IOException(const std::string& str): Exception(str) { }

            template <typename CharT,typename Traits>
            IOException(std::basic_ios<CharT,Traits>& ios,const std::string& str): Exception(str) { ios.setstate(std::ios::failbit); }

            ExceptionCode code() const throw() { return IO_EXCPT; }
        };

        struct OPENMEEGMATHS_EXPORT BadFile: public IOException {
     
            BadFile(const std::string& fmtname):                  IOException(message(fmtname))    { }
            BadFile(std::istream& is,const std::string& fmtname): IOException(is,message(fmtname)) { }

            ExceptionCode code() const throw() { return BAD_FILE; }

        private:

            static std::string message(const std::string& fmtname) { return std::string("Unable to read the input as a "+fmtname+" image file."); }
        };

        struct OPENMEEGMATHS_EXPORT BadFileOpening: public IOException {
     
            typedef enum { READ, WRITE } Mode;

            BadFileOpening(const std::string& file,const Mode& mode):                  IOException(message(file,mode))    { }
            BadFileOpening(std::istream& is,const std::string& file,const Mode& mode): IOException(is,message(file,mode)) { }

            ExceptionCode code() const throw() { return BAD_FILE_OPEN; }

        private:

            static std::string message(const std::string& file,const Mode& mode) {
                return std::string("Unable to open the file ")+file+" for "+((mode==READ) ? "reading" : "writing");
            }
        };

        struct OPENMEEGMATHS_EXPORT BadContent: public IOException {

            BadContent(const std::string& fmtname,const std::string& object):                  IOException(message(fmtname,object))    { }
            BadContent(std::istream& is,const std::string& fmtname,const std::string& object): IOException(is,message(fmtname,object)) { }

            ExceptionCode code() const throw() { return BAD_CONTENT; }

        private:

            static std::string message(const std::string& fmtname,const std::string& object) {
                return std::string("This "+fmtname+" file does not contains a "+object+" as expected.");
            }
        };

        struct OPENMEEGMATHS_EXPORT NoSuffix: public IOException {
            NoSuffix(const std::string& name): IOException(std::string("No identifiable suffix in name ")+name) { }

            ExceptionCode code() const throw() { return NO_SUFFIX; }
        };

        struct OPENMEEGMATHS_EXPORT BadHeader: public IOException {

            BadHeader():                 IOException(message())    { }
            BadHeader(std::istream& is): IOException(is,message()) { }

            ExceptionCode code() const throw() { return BAD_HDR; }

        private:

            static std::string message() { return std::string("Bad file header."); }
        };

        struct OPENMEEGMATHS_EXPORT ImpossibleObjectIdentification: public IOException {

            ImpossibleObjectIdentification(const std::string& file):                  IOException(message(file))    { }
            ImpossibleObjectIdentification(std::istream& is,const std::string& file): IOException(is,message(file)) { }

            ExceptionCode code() const throw() { return IMPOSSIBLE_IDENTIFICATION; }

        private:

            static std::string message(const std::string& name) { return std::string("Impossible to identify the object in this file: ")+name; }
        };

        struct OPENMEEGMATHS_EXPORT BadStorageType: public IOException {

            BadStorageType(const std::string& file):                  IOException(message(file))    { }
            BadStorageType(std::istream& is,const std::string& file): IOException(is,message(file)) { }

            ExceptionCode code() const throw() { return BAD_STORAGE_TYPE; }

        private:

            static std::string message(const std::string& file) { return std::string("Bad storage type in file ")+file+"."; }
        };

        struct OPENMEEGMATHS_EXPORT BadData: public IOException {

            BadData(const std::string& fmtname):                  IOException(message(fmtname))    { }
            BadData(std::istream& is,const std::string& fmtname): IOException(is,message(fmtname)) { }

            ExceptionCode code() const throw() { return BAD_DATA; }

        private:

            static std::string message(const std::string& fmtname) { return std::string("Bad ")+fmtname+" file data."; }
        };

        struct OPENMEEGMATHS_EXPORT BadVector: public IOException {

            BadVector(const std::size_t dim): IOException(message(dim)) { }

            template <typename CharT,typename Traits>
            BadVector(std::basic_ios<CharT,Traits>& ios,const std::size_t dim): IOException(ios,message(dim)) { }

            ExceptionCode code() const throw() { return BAD_VECT; }

        private:

            static std::string message(const std::size_t dim) {
                std::ostringstream ost;
                ost << "Bad file (expected a vector, got a matrix with " << dim << " columns).";
                return ost.str();
            }
        };

        struct OPENMEEGMATHS_EXPORT BadSymmMatrix: public IOException {

            BadSymmMatrix(const std::size_t dim1,const std::size_t dim2): IOException(message(dim1,dim2)) { }

            template <typename CharT,typename Traits>
            BadSymmMatrix(std::basic_ios<CharT,Traits>& ios,const std::size_t dim1,const std::size_t dim2): IOException(ios,message(dim1,dim2)) { }

            ExceptionCode code() const throw() { return BAD_SYMM_MAT; }

        private:

            static std::string message(const std::size_t dim1,const std::size_t dim2) {
                std::ostringstream ost;
                ost << "Symmetric matrix is expected to be square (got an " << dim1 << 'x' << dim2 << " matrix instead).";
                return ost.str();
            }
        };

        struct OPENMEEGMATHS_EXPORT NoIO: public IOException {

            typedef enum { READ, WRITE } Mode;

            NoIO(const std::string& file,const Mode& mode): IOException(message(file,mode)) { }

            template <typename CharT,typename Traits>
            NoIO(std::basic_ios<CharT,Traits>& ios,const std::string& file,const Mode& mode): IOException(ios,message(file,mode)) { }

            ExceptionCode code() const throw() { return NO_IO; }

        private:

            static std::string message(const std::string& file,const Mode& mode) {
                return std::string("Unable to find ")+((mode==READ) ? "reader" : "writer")+" for file "+file+".";
            }
        };

        struct OPENMEEGMATHS_EXPORT MatioError: public IOException {
            MatioError(const std::string& err): IOException(err) { }
            ExceptionCode code() const throw() { return MATIO_ERROR; }
        };

        struct OPENMEEGMATHS_EXPORT UnknownFileFormat: public IOException {

            UnknownFileFormat(const std::string& fmt): IOException(std::string("Unknown ")+fmt+" format.") { }

            template <typename CharT,typename Traits>
            UnknownFileFormat(std::basic_ios<CharT,Traits>& ios): IOException(ios,std::string("Unknown file format.")) { }

            ExceptionCode code() const throw() { return UNKN_FILE_FMT; }
        };

        struct OPENMEEGMATHS_EXPORT UnknownFileSuffix: public IOException {
            UnknownFileSuffix(const std::string& fmt): IOException(std::string("Unknown ")+fmt+" suffix.") { }
            ExceptionCode code() const throw() { return UNKN_FILE_SUFFIX; }
        };

        struct OPENMEEGMATHS_EXPORT UnknownNamedFileFormat: public Exception {
            UnknownNamedFileFormat(const std::string& name): Exception(std::string("Unknown format for file "+name+".")) { }
            ExceptionCode code() const throw() { return UNKN_NAMED_FILE_FMT; }
        };
    }
}

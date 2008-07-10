#ifndef MATLIB_EXCEPTIONS_H
#define MATLIB_EXCEPTIONS_H

//  Largely derivate work from Images/src/Exceptions.H (merge if possible).

#include <string>
#include <sstream>
#include <iostream>
#include <exception>

namespace maths {

    typedef enum { UNEXPECTED = 128, IO_EXCPT,
                   BAD_FMT, NO_SUFFIX, NON_MATCH_FMT, BAD_HDR, BAD_DATA, BAD_DIM, UNKN_DIM, BAD_SIZE_SPEC, UNKN_PIX, UNKN_PIX_TYPE,
                   UNKN_FILE_FMT, UNKN_FILE_SUFFIX, UNKN_NAMED_FILE_FMT, NON_MATCH_NAMED_FILE_FMT, NO_FILE_FMT,
                   BAD_PLGIN_LIST, BAD_PLGIN_FILE, BAD_PLGIN, ALREADY_KN_TAG } ExceptionCode;

    class Exception: public std::exception {
    public:

        Exception(const std::string& w): whatstring(std::string("Images::Exception: ")+w) { }

        virtual ~Exception() throw() { }

        virtual const char*   what () const throw() { return whatstring.c_str(); }
        virtual ExceptionCode code()  const throw() = 0;

    private:

        std::string whatstring;
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

    struct BadFormat: public IOException {

        BadFormat(const std::string& fmtname):                  IOException(message(fmtname))    { }
        BadFormat(std::istream& is,const std::string& fmtname): IOException(is,message(fmtname)) { }

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

    struct BadDimension: public IOException {

        BadDimension(const std::size_t dim): IOException(message(dim)) { }

        template <typename CharT,typename Traits>
        BadDimension(std::basic_ios<CharT,Traits>& ios,const std::size_t dim): IOException(ios,message(dim)) { }

        ExceptionCode code() const throw() { return BAD_DIM; }

    private:

        static std::string message(const std::size_t dim) {
            std::ostringstream ost;
            ost << "Incoming image is not of dimension " << dim << '.';
            return ost.str();
        }
    };

    struct BadSizeSpec: public IOException {

        BadSizeSpec(const std::size_t dim): IOException(message(dim)) { }

        template <typename CharT,typename Traits>
        BadSizeSpec(std::basic_ios<CharT,Traits>& ios,const std::size_t dim): IOException(ios,message(dim)) { }

        ExceptionCode code() const throw() { return BAD_SIZE_SPEC; }

    private:

        static std::string message(const std::size_t dim) {
            std::ostringstream ost;
            ost << "Size specification is not of dimension " << dim << '.';
            return ost.str();
        }
    };

    struct UnknownDimension: public Exception {

        UnknownDimension(const std::size_t dim): Exception(message(dim)) { }

        ExceptionCode code() const throw() { return UNKN_DIM; }

    private:

        static std::string message(const std::size_t dim) {
            std::ostringstream ost;
            ost << "Image is of unknown dimension " << dim << '.';
            return ost.str();
        }
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
}

#endif // !MATLIB_EXCEPTIONS_H

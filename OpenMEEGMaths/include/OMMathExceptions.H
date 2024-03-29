#pragma once

//  Largely derivate work from Images/src/Exceptions.H (merge if possible).

#include <string>
#include <sstream>
#include <iostream>
#include <exception>

#include "range.h"
#include "OpenMEEGMaths_Export.h"

#if defined(_MSC_VER)
    #pragma warning(disable : 4275) /* non dll-interface class X used as base for dll-interface struct Y */
#endif


namespace OpenMEEG {
    namespace maths {

        typedef enum { UNEXPECTED = 128, IO_EXCPT,
                       BAD_FILE, BAD_FILE_OPEN, BAD_CONTENT, NO_SUFFIX, BAD_HDR, BAD_DATA, BAD_VECT, UNKN_DIM, BAD_SYMM_MAT,
                       BAD_STORAGE_TYPE, NO_IO, NON_EXISTING_RANGE, OVERLAPPING_RANGES, NON_EXISTING_BLOCK, MATIO_ERROR,
                       UNKN_FILE_FMT, UNKN_FILE_SUFFIX, NO_FILE_FMT, UNKN_NAMED_FILE_FMT, IMPOSSIBLE_IDENTIFICATION,
                       BAD_SPARSE_OPERATION, BAD_LINEAR_ALGEBRA
                     } ExceptionCode;

        class Exception: public std::exception {
        public:

            Exception(const std::string& w): whatstring(std::string("Exception: ")+w) { }

            virtual ~Exception() noexcept { }

            virtual const char*   what () const noexcept { return whatstring.c_str(); }
            virtual ExceptionCode code()  const noexcept = 0;

        private:

            std::string whatstring;
        };

        struct UnexpectedException: public Exception {

            UnexpectedException(const std::string& func,const std::string& file,const unsigned line): Exception(message(func,file,line)) { }

            ExceptionCode code() const noexcept { return UNEXPECTED; }

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

            ExceptionCode code() const noexcept { return IO_EXCPT; }
        };

        struct BadFile: public IOException {

            BadFile(const std::string& fmtname):                  IOException(message(fmtname))    { }
            BadFile(std::istream& is,const std::string& fmtname): IOException(is,message(fmtname)) { }

            ExceptionCode code() const noexcept { return BAD_FILE; }

        private:

            static std::string message(const std::string& fmtname) { return std::string("Unable to read the input as a "+fmtname+" image file."); }
        };

        struct BadFileOpening: public IOException {

            typedef enum { READ, WRITE } Mode;

            BadFileOpening(const std::string& file,const Mode& mode):                  IOException(message(file,mode))    { }
            BadFileOpening(std::istream& is,const std::string& file,const Mode& mode): IOException(is,message(file,mode)) { }

            ExceptionCode code() const noexcept { return BAD_FILE_OPEN; }

        private:

            static std::string message(const std::string& file,const Mode& mode) {
                return std::string("Unable to open the file ")+file+" for "+((mode==READ) ? "reading" : "writing");
            }
        };

        struct BadContent: public IOException {

            BadContent(const std::string& fmtname,const std::string& object):                  IOException(message(fmtname,object))    { }
            BadContent(std::istream& is,const std::string& fmtname,const std::string& object): IOException(is,message(fmtname,object)) { }

            ExceptionCode code() const noexcept { return BAD_CONTENT; }

        private:

            static std::string message(const std::string& fmtname,const std::string& object) {
                return std::string("This "+fmtname+" file does not contains a "+object+" as expected.");
            }
        };

        struct NoSuffix: public IOException {
            NoSuffix(const std::string& name): IOException(std::string("No identifiable suffix in name ")+name) { }

            ExceptionCode code() const noexcept { return NO_SUFFIX; }
        };

        struct BadHeader: public IOException {

            BadHeader():                 IOException(message())    { }
            BadHeader(std::istream& is): IOException(is,message()) { }

            ExceptionCode code() const noexcept { return BAD_HDR; }

        private:

            static std::string message() { return std::string("Bad file header."); }
        };

        struct ImpossibleObjectIdentification: public IOException {

            ImpossibleObjectIdentification(const std::string& file):                  IOException(message(file))    { }
            ImpossibleObjectIdentification(std::istream& is,const std::string& file): IOException(is,message(file)) { }

            ExceptionCode code() const noexcept { return IMPOSSIBLE_IDENTIFICATION; }

        private:

            static std::string message(const std::string& name) { return std::string("Impossible to identify the object in this file: ")+name; }
        };

        struct BadStorageType: public IOException {

            BadStorageType(const std::string& file):                  IOException(message(file))    { }
            BadStorageType(std::istream& is,const std::string& file): IOException(is,message(file)) { }

            ExceptionCode code() const noexcept { return BAD_STORAGE_TYPE; }

        private:

            static std::string message(const std::string& file) { return std::string("Bad storage type in file ")+file+"."; }
        };

        struct BadData: public IOException {

            BadData(const std::string& fmtname):                  IOException(message(fmtname))    { }
            BadData(std::istream& is,const std::string& fmtname): IOException(is,message(fmtname)) { }

            ExceptionCode code() const noexcept { return BAD_DATA; }

        private:

            static std::string message(const std::string& fmtname) { return std::string("Bad ")+fmtname+" file data."; }
        };

        struct BadVector: public IOException {

            BadVector(const std::size_t dim): IOException(message(dim)) { }

            template <typename CharT,typename Traits>
            BadVector(std::basic_ios<CharT,Traits>& ios,const std::size_t dim): IOException(ios,message(dim)) { }

            ExceptionCode code() const noexcept { return BAD_VECT; }

        private:

            static std::string message(const std::size_t dim) {
                std::ostringstream ost;
                ost << "Bad file (expected a vector, got a matrix with " << dim << " columns).";
                return ost.str();
            }
        };

        struct BadSymmMatrix: public IOException {

            BadSymmMatrix(const std::size_t dim1,const std::size_t dim2): IOException(message(dim1,dim2)) { }

            template <typename CharT,typename Traits>
            BadSymmMatrix(std::basic_ios<CharT,Traits>& ios,const std::size_t dim1,const std::size_t dim2): IOException(ios,message(dim1,dim2)) { }

            ExceptionCode code() const noexcept { return BAD_SYMM_MAT; }

        private:

            static std::string message(const std::size_t dim1,const std::size_t dim2) {
                std::ostringstream ost;
                ost << "Symmetric matrix is expected to be square (got an " << dim1 << 'x' << dim2 << " matrix instead).";
                return ost.str();
            }
        };

        struct NoIO: public IOException {

            typedef enum { READ, WRITE } Mode;

            NoIO(const std::string& file,const Mode& mode): IOException(message(file,mode)) { }

            template <typename CharT,typename Traits>
            NoIO(std::basic_ios<CharT,Traits>& ios,const std::string& file,const Mode& mode): IOException(ios,message(file,mode)) { }

            ExceptionCode code() const noexcept { return NO_IO; }

        private:

            static std::string message(const std::string& file,const Mode& mode) {
                return std::string("Unable to find ")+((mode==READ) ? "reader" : "writer")+" for file "+file+".";
            }
        };

        struct OPENMEEGMATHS_EXPORT NonExistingRange: public Exception {

            NonExistingRange(const Range& range): Exception(message(range)) { }

            ExceptionCode code() const noexcept { return NON_EXISTING_RANGE; }

        private:

            static std::string message(const Range& range) {
                std::ostringstream ost;
                ost << "Attempt to access a non-existing range " << range << '.';
                return ost.str();
            }
        };

        struct OPENMEEGMATHS_EXPORT OverlappingRanges: public Exception {

            OverlappingRanges(const Range& r1,const Range& r2): Exception(message(r1,r2)) { }

            ExceptionCode code() const noexcept { return OVERLAPPING_RANGES; }

        private:

            static std::string message(const Range& r1,const Range& r2) {
                std::ostringstream ost;
                ost << "Ranges " << r1 << " and " << r2 << " are overlapping blocks in a block matrix.";
                return ost.str();
            }
        };

        struct OPENMEEGMATHS_EXPORT NonExistingBlock: public Exception {

            NonExistingBlock(const unsigned ind): Exception(message(ind)) { }

            ExceptionCode code() const noexcept { return NON_EXISTING_BLOCK; }

        private:

            static std::string message(const unsigned ind) {
                std::ostringstream ost;
                ost << "Attempt to access a non-existing block at index " << ind << '.';
                return ost.str();
            }
        };

        struct MatioError: public IOException {
            MatioError(const std::string& err): IOException(err) { }
            ExceptionCode code() const noexcept { return MATIO_ERROR; }
        };

        struct UnknownFileFormat: public IOException {

            UnknownFileFormat(const std::string& fmt): IOException(std::string("Unknown ")+fmt+" format.") { }

            template <typename CharT,typename Traits>
            UnknownFileFormat(std::basic_ios<CharT,Traits>& ios): IOException(ios,std::string("Unknown file format.")) { }

            ExceptionCode code() const noexcept { return UNKN_FILE_FMT; }
        };

        struct UnknownFileSuffix: public IOException {
            UnknownFileSuffix(const std::string& fmt): IOException(std::string("Unknown ")+fmt+" suffix.") { }
            ExceptionCode code() const noexcept { return UNKN_FILE_SUFFIX; }
        };

        struct UnknownNamedFileFormat: public Exception {
            UnknownNamedFileFormat(const std::string& name): Exception(std::string("Unknown format for file "+name+".")) { }
            ExceptionCode code() const noexcept { return UNKN_NAMED_FILE_FMT; }
        };

        struct BadSparseOperation: public Exception {
            BadSparseOperation(const std::string& what): Exception(what) { }
            ExceptionCode code() const noexcept { return BAD_SPARSE_OPERATION; }
        };

        struct LinearAlgebraError: public Exception {
            LinearAlgebraError(const std::string& what): Exception(what) { }
            ExceptionCode code() const noexcept { return BAD_LINEAR_ALGEBRA; }
        };
    }
}

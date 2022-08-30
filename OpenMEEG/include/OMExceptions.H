// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <mutex>
#include <exception>

// TODO We should rename this file Exceptions.H probably, but this collides
// with OpenMEEGMaths Exceptions.H. This file and that file also contain many
// shared error types. This should be refactored at some point.

namespace OpenMEEG {

    //TODO  remove all unused exception
#if 0
    struct ReturnCode {
        static const unsigned code =
    };
#endif

    typedef enum { UNEXPECTED = 128, IO_EXCPT, OPEN_ERROR,
                   BAD_FMT, NO_SUFFIX, NON_MATCH_FMT, BAD_HDR, BAD_DATA, BAD_CONTENT, WRONG_FILE_FMT, BAD_DIM, UNKN_DIM, BAD_SIZE_SPEC, UNKN_PIX, UNKN_PIX_TYPE,
                   UNKN_FILE_FMT, UNKN_FILE_SUFFIX, UNKN_NAMED_FILE_FMT, NON_MATCH_NAMED_FILE_FMT, NO_FILE_FMT,
                   BAD_PLGIN_LIST, BAD_PLGIN_FILE, BAD_PLGIN, ALREADY_KN_TAG, NON_EXISTING_DOMAIN, UNKNOWN_VERTEX,
                   SOURCE_MESH_OVERLAPS_GEOMETRY, NO_IMG_ARG, DIFF_IMG, BAD_VTK, BAD_SENSOR, BAD_DIPOLE, BAD_GENERIC } ExceptionCode;

    class Exception: public std::exception {
    public:

        Exception(const std::string& w): whatstring(w) { }

        virtual ~Exception() noexcept { }

        virtual const char*   what() const noexcept { return whatstring.c_str(); }
        virtual ExceptionCode code() const noexcept = 0;

        std::string& what() { return whatstring; }

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

    struct OpenError: public IOException {

        OpenError(const std::string& filename): IOException(message(filename)) { }
        OpenError(const std::string& filename,const std::string& fmtname): IOException(message(filename,fmtname)) { }

        ExceptionCode code() const noexcept { return OPEN_ERROR; }

    private:

        static std::string message(const std::string& filename) { return std::string("Cannot open file "+filename+" !"); }

        static std::string message(const std::string& filename,const std::string& fmtname) {
            return std::string("Error opening "+fmtname+" file: "+filename+" for reading.");
        }
    };

    struct BadFormat: public IOException {

        BadFormat(const std::string& fmtname):                  IOException(message(fmtname))    { }
        BadFormat(std::istream& is,const std::string& fmtname): IOException(is,message(fmtname)) { }

        ExceptionCode code() const noexcept { return BAD_FMT; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Unable to read the input as a "+fmtname+" image file."); }
    };

    struct BadInterface: public Exception {

        BadInterface(const std::string& str): Exception(message(str))    { }

        ExceptionCode code() const noexcept { return BAD_FMT; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Unable to find interface " + fmtname + "."); }
    };

    struct BadDomain: public Exception {

        BadDomain(const std::string& str): Exception(message(str))    { }

        ExceptionCode code() const noexcept { return BAD_FMT; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Unable to find domain " + fmtname + "."); }
    };

    struct BadInterfaceIndex: public Exception {

        BadInterfaceIndex(const std::string& str): Exception(message(str))    { }

        ExceptionCode code() const noexcept { return BAD_FMT; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Unable to read the input as a "+fmtname+" image file."); }
    };

    struct NoSuffix: public IOException {
        NoSuffix(const std::string& name): IOException(std::string("No identifiable suffix in name ")+name) { }

        ExceptionCode code() const noexcept { return NO_SUFFIX; }
    };

    struct NonMatchingFormat: public IOException {

        NonMatchingFormat(const std::string& fmtname):                  IOException(message(fmtname))    { }
        NonMatchingFormat(std::ostream& os,const std::string& fmtname): IOException(os,message(fmtname)) { }

        ExceptionCode code() const noexcept { return NON_MATCH_FMT; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Unable to save the image in the ")+fmtname+" format."; }
    };

    struct BadHeader: public IOException {

        BadHeader(const std::string& fmtname=""):                  IOException(message(fmtname))    { }
        BadHeader(std::istream& is,const std::string& fmtname=""): IOException(is,message(fmtname)) { }

        ExceptionCode code() const noexcept { return BAD_HDR; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Bad ")+fmtname+" file header."; }
    };

    struct BadData: public IOException {

        BadData(const std::string& fmtname):                  IOException(message(fmtname))    { }
        BadData(std::istream& is,const std::string& fmtname): IOException(is,message(fmtname)) { }

        ExceptionCode code() const noexcept { return BAD_DATA; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Bad ")+fmtname+" file data."; }
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

    struct WrongFileFormat: public IOException {

        WrongFileFormat(const std::string& file): IOException(std::string("Wrong format for file ")+file+".") { }

        ExceptionCode code() const noexcept { return WRONG_FILE_FMT; }
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

    struct NonMatchingNamedFileFormat: public Exception {
        NonMatchingNamedFileFormat(const std::string& name,const std::string& fmt): Exception(std::string("Cannot save the image "+name+" as a "+fmt+" file.")) { }

        ExceptionCode code() const noexcept { return NON_MATCH_NAMED_FILE_FMT; }
    };

    struct NoMatchingFileFormat: public IOException {
        NoMatchingFileFormat(std::ostream& os): IOException(os,std::string("No IO is able to write this type of image.")) { }

        ExceptionCode code() const noexcept { return NO_FILE_FMT; }
    };

    template <typename T>
    struct NonExistingDomain: public std::exception {
        NonExistingDomain(const std::string& dom,const T& id): domain(dom), interface(id) { }

        virtual ~NonExistingDomain() noexcept { }

        const char* what() const noexcept {
            std::ostringstream oss;
            oss << "Definition of domain " << domain << ": interface " << interface << " does not exist";
            std::string temp = oss.str();
            const char* out = temp.c_str();
            return out;
        }

        ExceptionCode code() const noexcept { return NON_EXISTING_DOMAIN; }

        const std::string domain;
        const T&          interface;
    };

    struct UnknownVertex: public Exception {
        UnknownVertex(const std::string& msg): Exception(message(msg)) { }

        ExceptionCode code() const noexcept { return UNKNOWN_VERTEX; }

    private:

        static std::string message(const std::string& msg) { return std::string("Unknown vertex in triangle: " + msg); }
    };

    struct OverlappingSourceMesh: public Exception {
        OverlappingSourceMesh(): Exception("Source mesh overlaps the geometry") { }

        ExceptionCode code() const noexcept { return SOURCE_MESH_OVERLAPS_GEOMETRY; }
    };

    //  Revisit (TODO).

    struct NoImageArgument: public Exception {
        NoImageArgument(): Exception("No image argument.") { }

        ExceptionCode code() const noexcept { return NO_IMG_ARG; }
    };

    struct DifferentImages: public Exception {
        DifferentImages(): Exception("Different image types.") { }

        ExceptionCode code() const noexcept { return DIFF_IMG; }
    };

    // VTK

    struct VTKError: public Exception {

        VTKError(const std::string& str): Exception(message(str)) { }
        ExceptionCode code() const noexcept { return BAD_VTK; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("VTK error: " + fmtname); }
    };

    // Sensors

    struct SensorError: public Exception {

        SensorError(const std::string& str): Exception(message(str)) { }
        ExceptionCode code() const noexcept { return BAD_SENSOR; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Sensor error: " + fmtname); }
    };

    // Dipoles

    struct DipoleError: public Exception {

        DipoleError(const std::string& str): Exception(message(str)) { }
        ExceptionCode code() const noexcept { return BAD_DIPOLE; }

    private:

        static std::string message(const std::string& fmtname) { return std::string("Dipoles file format error: " + fmtname); }
    };

    // Danielsson

    struct GenericError: public Exception {

        GenericError(const std::string& str): Exception(str) { }

        ExceptionCode code() const noexcept { return BAD_GENERIC; }
    };

    // OpenMP
    // https://stackoverflow.com/questions/11828539/elegant-exceptionhandling-in-openmp
    class ThreadException {
        std::exception_ptr Ptr;
        std::mutex         Lock;
    public:
        ThreadException(): Ptr(nullptr) {}
        void Rethrow(){
            if(this->Ptr) std::rethrow_exception(this->Ptr);
        }
        void CaptureException() {
            std::unique_lock<std::mutex> guard(this->Lock);
            this->Ptr = std::current_exception();
        }
        template <typename Function, typename... Parameters>
        void Run(Function f, Parameters... params)
        {
            try
            {
                f(params...);
            }
            catch (...)
            {
                CaptureException();
            }
        }
    };
}

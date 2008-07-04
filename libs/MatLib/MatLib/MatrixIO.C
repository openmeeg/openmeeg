#include "MatrixIO.H"

namespace Maths {

    namespace Internal {

        //  Read a few bytes to figure out the file format and put them back into the stream.
        static const unsigned maxtagsize = 32;

        static const char*
        ReadTag(std::istream& is) throw(std::string) {

            static char buffer[maxtagsize];

            try {
                is.read(buffer,maxtagsize);
            } catch(...) {
                throw std::string("BadHeader(is)");
            }

            for(int i=maxtagsize-1;i>=0;--i)
                is.putback(buffer[i]);

            return buffer;
        }
    }

    MathsIO::IO MathsIO::DefaultIO = 0;
    bool MathsIO::permanent = false;

    const MathsIO::IO& MathsIO::format(const std::string& fmt) throw(UnknownFileFormat) {
        for (IOs::const_iterator i=ios().begin();i!=ios().end();++i)
            if (fmt==(*i)->identity())
                return *i;
        throw UnknownFileFormat(std::string("Unknown file format : ")+fmt);
    }

    const MathsIO::IO& MathsIO::format_from_suffix(const std::string& name) throw(NoSuffix,UnknownFileSuffix) {
        const std::string::size_type pos = name.find_last_of(".");
        if (pos==std::string::npos)
            throw NoSuffix(std::string("Could find suffix for :")+name);

        const std::string suffix = name.substr(pos+1);
        for(IOs::const_iterator i=ios().begin();i!=ios().end();++i) {
            if ((*i)->known_suffix(suffix.c_str()))
                return *i;
        }
        throw UnknownFileSuffix(std::string("Unkown suffix :")+suffix);
    }

    Maths::ifstream& operator>>(Maths::ifstream& mio,MatrixBase& matrix) throw(std::string) {
        std::ifstream is(mio.name().c_str());
        if(is.fail()) {
            throw std::string("Unable to open : ")+mio.name();
        }

        const char* buffer = Internal::ReadTag(is);

        if (Maths::MathsIO::IO io = Maths::MathsIO::default_io()) {
            if (io->identify(std::string(buffer))) {
                io->setName(mio.name());
                io->read(is,matrix);
                return mio;
            }
        } else {
            for (Maths::MathsIO::IOs::const_iterator io=Maths::MathsIO::ios().begin();io!=Maths::MathsIO::ios().end();++io) {
                if ((*io)->identify(std::string(buffer))) {
                    (*io)->setName(mio.name());
                    (*io)->read(is,matrix);
                    return mio;
                }
            }
        }
        throw std::string("Unable to find proper reader.");
    }

    Maths::ofstream& operator<<(Maths::ofstream& mio,const MatrixBase& matrix) throw(std::string) {
        std::ofstream os(mio.name().c_str());
        if(os.fail()) {
            throw std::string("Unable to open : ")+mio.name();
        }

        if (Maths::MathsIO::IO io = Maths::MathsIO::default_io()) {
            if (io->known(matrix)) {
                io->setName(mio.name());
                io->write(os,matrix);
                return mio;
            }
        } else {
            for (Maths::MathsIO::IOs::const_iterator io=Maths::MathsIO::ios().begin();io!=Maths::MathsIO::ios().end();++io) {
                if ((*io)->known(matrix)) {
                    (*io)->setName(mio.name());
                    (*io)->write(os,matrix);
                    return mio;
                }
            }
        }
        throw std::string("Unable to find proper writer.");
    }


}

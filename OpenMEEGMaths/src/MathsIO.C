#include "MathsIO.H"

namespace OpenMEEG {

    namespace maths {

        namespace Internal {

            //  Read a few bytes to figure out the file format and put them back into the stream.
            static const unsigned maxtagsize = 32;

            static const char*
            ReadTag(std::istream& is) {

                static char buffer[maxtagsize];

                try {
                    is.read(buffer,maxtagsize);
                } catch(...) {
                    throw BadHeader();
                }

                for(int i=maxtagsize-1;i>=0;--i)
                    is.putback(buffer[i]);

                return buffer;
            }
        }

        MathsIO::IO MathsIO::DefaultIO = 0;
        bool MathsIO::permanent = false;

        const MathsIO::IO& MathsIO::format(const std::string& fmt) {
            for (IOs::const_iterator i=ios().begin();i!=ios().end();++i) {
                if (fmt==(*i)->identity())
                    return *i;
            }
            throw UnknownFileFormat(fmt);
        }

        const MathsIO::IO& MathsIO::format_from_suffix(const std::string& name) {
            const std::string::size_type pos = name.find_last_of(".");
            if (pos==std::string::npos)
                throw NoSuffix(name);

            const std::string suffix = name.substr(pos+1);
            for (IOs::const_iterator i=ios().begin();i!=ios().end();++i)
                if ((*i)->known_suffix(suffix.c_str()))
                    return *i;
            throw UnknownFileSuffix(suffix);
        }

        maths::ifstream& operator>>(maths::ifstream& mio,LinOp& linop) {
            std::ifstream is(mio.name().c_str(),std::ios::binary);
            if (is.fail())
                throw BadFileOpening(mio.name(),BadFileOpening::READ);

            const char* buffer = Internal::ReadTag(is);

            if (maths::MathsIO::IO dio = maths::MathsIO::GetCurrentFormat()) {
                if (dio->identify(std::string(buffer))) {
                    dio->setName(mio.name());
                    dio->read(is,linop);
                    linop.default_io() = dio;
                    return mio;
                }
            } else {
                for (maths::MathsIO::IOs::const_iterator io=maths::MathsIO::ios().begin();io!=maths::MathsIO::ios().end();++io) {
                    if ((*io)->identify(std::string(buffer))) {
                        (*io)->setName(mio.name());
                        (*io)->read(is,linop);
                        linop.default_io() = *io;
                        return mio;
                    }
                }
            }
            throw NoIO(mio.name(),NoIO::READ);
        }

        maths::ofstream& operator<<(maths::ofstream& mio,const LinOp& linop) {

            std::ofstream os(mio.name().c_str(),std::ios::binary);
            if (os.fail())
                throw BadFileOpening(mio.name(),BadFileOpening::WRITE);

            if (maths::MathsIO::IO dio = maths::MathsIO::GetCurrentFormat()) {
                if (dio->known(linop)) {
                    dio->setName(mio.name());
                    dio->write(os,linop);
                    return mio;
                }
            } else {
                for (maths::MathsIO::IOs::const_iterator io=maths::MathsIO::ios().begin();io!=maths::MathsIO::ios().end();++io) {
                    if ((*io)->known(linop)) {
                        (*io)->setName(mio.name());
                        (*io)->write(os,linop);
                        return mio;
                    }
                }
            }
            throw NoIO(mio.name(),NoIO::WRITE);
        }

        LinOpInfo info(const char* name) {
            std::ifstream is(name,std::ios::binary);
            if(is.fail())
                throw BadFileOpening(name,BadFileOpening::READ);

            const char* buffer = Internal::ReadTag(is);

            if (maths::MathsIO::IO dio = maths::MathsIO::default_io()) {
                if (dio->identify(std::string(buffer))) {
                    dio->setName(name);
                    return dio->info(is);
                }
            } else {
                for (maths::MathsIO::IOs::const_iterator io=maths::MathsIO::ios().begin();io!=maths::MathsIO::ios().end();++io) {
                    if ((*io)->identify(std::string(buffer))) {
                        (*io)->setName(name);
                        return (*io)->info(is);
                    }
                }
            }
            throw NoIO(name,NoIO::READ);
        }
    }
}

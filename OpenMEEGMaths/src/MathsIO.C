#include "MathsIO.H"

namespace OpenMEEG {

    namespace maths {

        namespace Internal {

            //  Read a few bytes to figure out the file format and put them back into the stream.

            static const unsigned maxtagsize = 32; // Reserve space for adding a null char at the end.

            static std::string
            ReadTag(std::istream& is) {

                char buffer[maxtagsize+1]; // Not static: shared would corrupt concurrent readers.

                try {
                    is.read(buffer,maxtagsize);

                    for (int i=maxtagsize-1; i>=0; --i)
                        is.putback(buffer[i]);

                    buffer[maxtagsize] = '\0'; // Add an end of string.

                } catch(...) {
                    throw BadHeader();
                }

                return std::string(buffer);
            }
        }

        thread_local MathsIO::IO MathsIO::DefaultIO = 0;
        thread_local bool MathsIO::permanent = false;

        const MathsIO::IO& MathsIO::format(const std::string& fmt) {
            for (const IO& io : ios())
                if (fmt==io->identity())
                    return io;
            throw UnknownFileFormat(fmt);
        }

        const MathsIO::IO& MathsIO::format_from_suffix(const std::string& name) {
            const std::string::size_type pos = name.find_last_of(".");
            if (pos==std::string::npos)
                throw NoSuffix(name);

            const std::string suffix = name.substr(pos+1);
            for (const IO& io : ios())
                if (io->known_suffix(suffix.c_str()))
                    return io;
            throw UnknownFileSuffix(suffix);
        }

        maths::ifstream& operator>>(maths::ifstream& mio,LinOp& linop) {
            std::ifstream is(mio.name().c_str(),std::ios::binary);
            if (is.fail())
                throw BadFileOpening(mio.name(),BadFileOpening::READ);

            const std::string& buffer = Internal::ReadTag(is);

            if (maths::MathsIO::IO dio = maths::MathsIO::GetCurrentFormat()) {
                if (dio->identify(buffer)) {
                    dio->read(mio.name(),is,linop);
                    linop.default_io() = dio;
                    return mio;
                }
            } else {
                for (const maths::MathsIO::IO& io : MathsIO::ios())
                    if (io->identify(buffer)) {
                        io->read(mio.name(),is,linop);
                        linop.default_io() = io;
                        return mio;
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
                    dio->write(mio.name(),os,linop);
                    return mio;
                }
            } else {
                for (const maths::MathsIO::IO& io : MathsIO::ios())
                    if (io->known(linop)) {
                        io->write(mio.name(),os,linop);
                        return mio;
                    }
            }
            throw NoIO(mio.name(),NoIO::WRITE);
        }

        LinOpInfo info(const char* name) {
            std::ifstream is(name,std::ios::binary);
            if (is.fail())
                throw BadFileOpening(name,BadFileOpening::READ);

            const std::string& buffer = Internal::ReadTag(is);

            if (maths::MathsIO::IO dio = maths::MathsIO::default_io()) {
                if (dio->identify(buffer)) {
                    return dio->info(name,is);
                }
            } else {
                for (const maths::MathsIO::IO& io : MathsIO::ios())
                    if (io->identify(buffer)) {
                        return io->info(name,is);
                    }
            }
            throw NoIO(name,NoIO::READ);
        }
    }
}

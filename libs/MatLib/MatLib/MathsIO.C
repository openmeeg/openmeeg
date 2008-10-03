#include "MathsIO.H"

//  An ugly hack to insert TrivialBinIO at the end of the IO list. 
//  Remove when TrivialBinIo is deleted.
#include <TrivialBinIO.H>

namespace maths {

    namespace Internal {

        //  Read a few bytes to figure out the file format and put them back into the stream.
        static const unsigned maxtagsize = 32;

        static const char*
        ReadTag(std::istream& is) throw(BadHeader) {

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

    //  An ugly hack to insert TrivialBinIO at the end of the IO list. 
    //  Remove when TrivialBinIo is deleted.

    void MathsIOBase::InsertTrivialBinIO() {
        static bool TrivialBinIOInserted = false;

        if (!TrivialBinIOInserted) {
            maths::TrivialBinIO* prototype = new maths::TrivialBinIO();
            maths::MathsIO::ios().push_back(prototype);
            TrivialBinIOInserted = true;
        }
    }

    const MathsIO::IO& MathsIO::format(const std::string& fmt) throw(UnknownFileFormat) {
        //  An ugly hack to insert TrivialBinIO at the end of the IO list. 
        //  Remove when TrivialBinIo is deleted.

        maths::MathsIOBase::InsertTrivialBinIO();

        for (IOs::const_iterator i=ios().begin();i!=ios().end();++i) {
            if (fmt==(*i)->identity())
                return *i;
        }

        throw UnknownFileFormat(fmt);
    }

    const MathsIO::IO& MathsIO::format_from_suffix(const std::string& name) throw(NoSuffix,UnknownFileSuffix) {
        //  An ugly hack to insert TrivialBinIO at the end of the IO list. 
        //  Remove when TrivialBinIo is deleted.

        maths::MathsIOBase::InsertTrivialBinIO();
        
        const std::string::size_type pos = name.find_last_of(".");
        if (pos==std::string::npos)
            throw NoSuffix(name);

        const std::string suffix = name.substr(pos+1);
        for(IOs::const_iterator i=ios().begin();i!=ios().end();++i) {
            if ((*i)->known_suffix(suffix.c_str()))
                return *i;
        }
        throw UnknownFileSuffix(suffix);
    }

    maths::ifstream& operator>>(maths::ifstream& mio,LinOp& linop) throw(BadFileOpening,NoIO) {

        //  An ugly hack to insert TrivialBinIO at the end of the IO list. 
        //  Remove when TrivialBinIo is deleted.

        maths::MathsIOBase::InsertTrivialBinIO();

        std::ifstream is(mio.name().c_str(),std::ios::binary);
        if(is.fail())
            throw BadFileOpening(mio.name(),BadFileOpening::READ);

        const char* buffer = Internal::ReadTag(is);

        if (maths::MathsIO::IO io = maths::MathsIO::default_io()) {
            if (io->identify(std::string(buffer))) {
                io->setName(mio.name());
                io->read(is,linop);
                return mio;
            }
        } else {
            for (maths::MathsIO::IOs::const_iterator io=maths::MathsIO::ios().begin();io!=maths::MathsIO::ios().end();++io) {
                if ((*io)->identify(std::string(buffer))) {
                    (*io)->setName(mio.name());
                    (*io)->read(is,linop);
                    return mio;
                }
            }
        }
        throw NoIO(mio.name(),NoIO::READ);
    }

    maths::ofstream& operator<<(maths::ofstream& mio,const LinOp& linop) throw(BadFileOpening,NoIO) {

        //  An ugly hack to insert TrivialBinIO at the end of the IO list. 
        //  Remove when TrivialBinIo is deleted.

        maths::MathsIOBase::InsertTrivialBinIO();

        std::ofstream os(mio.name().c_str(),std::ios::binary);
        if(os.fail())
            throw BadFileOpening(mio.name(),BadFileOpening::WRITE);

        if (maths::MathsIO::IO io = maths::MathsIO::default_io()) {
            if (io->known(linop)) {
                io->setName(mio.name());
                io->write(os,linop);
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


}
